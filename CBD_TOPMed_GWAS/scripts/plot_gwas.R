#!/usr/bin/env Rscript
# =============================================================================
# Manhattan + QQ plots, genomic inflation factor and top-hit table
#
# Rscript plot_gwas.R --sumstats sumstats.tsv.gz --outdir DIR
#        [--freq-cases freq_cases.afreq --freq-controls freq_controls.afreq]
#        [--suggestive 1e-5] [--title "CBD GWAS"]
#
# sumstats: columns #CHROM POS MarkerID ... P (output of 11_post_gwas.sh)
# Requires: data.table
# =============================================================================
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) {
    if (is.null(default)) stop("missing required argument ", flag, call. = FALSE)
    return(default)
  }
  args[i + 1]
}
ss_file    <- get_arg("--sumstats")
outdir     <- get_arg("--outdir")
fq_case    <- get_arg("--freq-cases", "")
fq_ctrl    <- get_arg("--freq-controls", "")
suggestive <- as.numeric(get_arg("--suggestive", "1e-5"))
title      <- get_arg("--title", "GWAS")
GWS <- 5e-8

ss <- if (grepl("\\.gz$", ss_file)) fread(cmd = paste("zcat", shQuote(ss_file))) else fread(ss_file)
setnames(ss, "#CHROM", "CHROM")
ss <- ss[!is.na(P) & P > 0]

# ---- genomic inflation --------------------------------------------------------
chisq  <- qchisq(ss$P, df = 1, lower.tail = FALSE)
lambda <- median(chisq) / qchisq(0.5, df = 1)
writeLines(sprintf("lambda_GC\t%.4f\nn_variants\t%d\nn_gws_p<5e-8\t%d\nn_suggestive_p<%g\t%d",
                   lambda, nrow(ss), sum(ss$P < GWS), suggestive, sum(ss$P < suggestive)),
           file.path(outdir, "lambda.txt"))
cat(sprintf("lambda_GC = %.4f over %d variants\n", lambda, nrow(ss)))

# ---- Manhattan ------------------------------------------------------------------
setorder(ss, CHROM, POS)
chr_len <- ss[, .(len = max(POS)), by = CHROM][order(CHROM)]
chr_len[, offset := cumsum(as.numeric(shift(len, fill = 0)))]
ss[chr_len, cum_pos := POS + i.offset, on = "CHROM"]
axis_pos <- ss[, .(mid = (min(cum_pos) + max(cum_pos)) / 2), by = CHROM]

# thin the uninformative bulk (P > 0.01) to keep the PNG light
set.seed(1)
plot_dt <- ss[P < 0.01 | runif(.N) < 0.1]
plot_dt[, logp := -log10(P)]
ymax <- max(8, ceiling(max(plot_dt$logp)) + 0.5)

png(file.path(outdir, "manhattan.png"), width = 3000, height = 1300, res = 220)
par(mar = c(4.5, 4.5, 3, 1))
plot(plot_dt$cum_pos, plot_dt$logp, pch = 20, cex = 0.45,
     col = ifelse(plot_dt$CHROM %% 2 == 0, "#3B6FB6", "#7FA7D9"),
     xaxt = "n", xlab = "Chromosome", ylab = expression(-log[10](italic(P))),
     ylim = c(0, ymax), xaxs = "i", bty = "l", main = title)
axis(1, at = axis_pos$mid, labels = axis_pos$CHROM, cex.axis = 0.75, tick = FALSE)
abline(h = -log10(GWS), col = "#C0392B", lty = 2)
abline(h = -log10(suggestive), col = "#7F8C8D", lty = 3)
invisible(dev.off())

# ---- QQ ---------------------------------------------------------------------------
obs <- sort(-log10(ss$P), decreasing = TRUE)
exp <- -log10(ppoints(length(obs)))
keep <- obs > 2 | runif(length(obs)) < 0.05
png(file.path(outdir, "qq.png"), width = 1400, height = 1400, res = 220)
par(mar = c(4.5, 4.5, 3, 1))
plot(exp[keep], obs[keep], pch = 20, cex = 0.5, col = "#3B6FB6", bty = "l",
     xlab = expression(Expected ~ -log[10](italic(P))),
     ylab = expression(Observed ~ -log[10](italic(P))),
     main = sprintf("%s\nlambda = %.3f", title, lambda))
abline(0, 1, col = "#C0392B")
invisible(dev.off())

# ---- top hits ---------------------------------------------------------------------
top <- ss[P < suggestive][order(P)]
read_freq <- function(f, label) {
  fq <- fread(f)
  setnames(fq, c("ID", "ALT_FREQS", "OBS_CT"), c("MarkerID", paste0("ALT_FREQ_", label), paste0("N_ALLELES_", label)))
  fq[, c("MarkerID", paste0("ALT_FREQ_", label), paste0("N_ALLELES_", label)), with = FALSE]
}
if (nrow(top) > 0 && nzchar(fq_case) && nzchar(fq_ctrl)) {
  top <- merge(top, read_freq(fq_case, "cases"), by = "MarkerID", all.x = TRUE, sort = FALSE)
  top <- merge(top, read_freq(fq_ctrl, "controls"), by = "MarkerID", all.x = TRUE, sort = FALSE)
  setorder(top, P)
}
top[, c("cum_pos") := NULL]
fwrite(top, file.path(outdir, "top_hits.tsv"), sep = "\t")
cat(sprintf("top hits (P < %g): %d -> %s\n", suggestive, nrow(top), file.path(outdir, "top_hits.tsv")))
