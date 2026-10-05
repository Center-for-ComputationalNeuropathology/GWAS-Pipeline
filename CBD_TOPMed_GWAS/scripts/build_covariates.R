#!/usr/bin/env Rscript
# =============================================================================
# Merge sample metadata with PLINK2 principal components -> GWAS covariate file
#
# Rscript build_covariates.R --metadata meta.txt --eigenvec pca.eigenvec \
#         --n-pcs 10 --out covariates.txt [--plot pca_plots.pdf]
#
# metadata : header with FID, IID, PHENO (1 case / 0 control) and any other
#            covariates (Sex, Age, CHIP, ...). Existing PC columns are replaced.
# Output   : FID IID <metadata covariates> PC1..PCn PHENO  (tab-separated)
#            Only samples present in BOTH files are written; the rest are
#            listed in <out>.missing_samples.txt.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) {
    if (is.null(default)) stop("missing required argument ", flag, call. = FALSE)
    return(default)
  }
  args[i + 1]
}
meta_file <- get_arg("--metadata")
eig_file  <- get_arg("--eigenvec")
n_pcs     <- as.integer(get_arg("--n-pcs", "10"))
out_file  <- get_arg("--out")
plot_file <- get_arg("--plot", "")

# ---- read inputs --------------------------------------------------------------
meta <- read.table(meta_file, header = TRUE, stringsAsFactors = FALSE,
                   colClasses = c(FID = "character", IID = "character"))
for (col in c("FID", "IID", "PHENO")) {
  if (!col %in% names(meta)) stop("metadata lacks column ", col, call. = FALSE)
}
meta <- meta[, !grepl("^PC[0-9]+$", names(meta)), drop = FALSE]

# PLINK2 writes "#FID IID PC1 ..." (or "#IID PC1 ..." when FIDs are absent)
eig <- read.table(eig_file, header = TRUE, comment.char = "", check.names = FALSE,
                  stringsAsFactors = FALSE)
names(eig) <- sub("^#", "", names(eig))
pc_cols <- paste0("PC", seq_len(n_pcs))
eig <- eig[, c("IID", pc_cols)]
eig$IID <- as.character(eig$IID)

# ---- merge (by IID, as FIDs differ between sources in some cohorts) -----------
merged <- merge(meta, eig, by = "IID", sort = FALSE)
missing_ids <- setdiff(meta$IID, eig$IID)

cat(sprintf("metadata samples : %d\n", nrow(meta)))
cat(sprintf("PCA samples      : %d\n", nrow(eig)))
cat(sprintf("merged           : %d\n", nrow(merged)))
if (length(missing_ids) > 0) {
  miss_file <- paste0(out_file, ".missing_samples.txt")
  writeLines(missing_ids, miss_file)
  cat(sprintf("WARNING: %d metadata samples have no genotypes/PCs -> %s\n",
              length(missing_ids), miss_file))
}

other <- setdiff(names(meta), c("FID", "IID", "PHENO"))
out <- merged[, c("FID", "IID", other, pc_cols, "PHENO")]
write.table(out, out_file, quote = FALSE, row.names = FALSE, sep = "\t")

cat("\nPHENO (1 = case, 0 = control):\n")
print(table(out$PHENO, useNA = "ifany"))
cat(sprintf("\nwritten: %s\n", out_file))

# ---- PCA plots ------------------------------------------------------------------
if (nzchar(plot_file)) {
  pdf(plot_file, width = 11, height = 5.5)
  par(mfrow = c(1, 2))
  colour_by <- function(v, title) {
    f <- factor(ifelse(is.na(v), "NA", as.character(v)))
    pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
    cols <- pal[(as.integer(f) - 1) %% length(pal) + 1]
    for (pair in list(c(1, 2), c(2, 3))) {
      x <- out[[pc_cols[pair[1]]]]; y <- out[[pc_cols[pair[2]]]]
      plot(x, y, col = cols, pch = 20, cex = 0.6,
           xlab = pc_cols[pair[1]], ylab = pc_cols[pair[2]], main = paste("coloured by", title))
      legend("topright", legend = levels(f), col = pal[seq_along(levels(f))], pch = 20, bty = "n", cex = 0.8)
    }
  }
  out_ph <- out$PHENO
  colour_by(ifelse(out_ph == 1, "case", ifelse(out_ph == 0, "control", NA)), "phenotype")
  if ("CHIP" %in% names(out)) colour_by(out$CHIP, "CHIP")
  invisible(dev.off())
  cat(sprintf("plots  : %s\n", plot_file))
}
