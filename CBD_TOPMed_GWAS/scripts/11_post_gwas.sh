#!/usr/bin/env bash
# =============================================================================
# Step 11 - Post-GWAS: summary statistics, plots, allele frequencies
#
# Input : WORK_DIR/gwas_<RUN_TAG>/GWAS_<RUN_TAG>_firth.PHENO.glm.firth
# Output: WORK_DIR/gwas_<RUN_TAG>/
#   CBD_GWAS_<RUN_TAG>_hg38_locuszoom.tsv.gz(.tbi)  upload-ready for my.locuszoom.org
#   manhattan.png, qq.png, lambda.txt               genomic inflation (lambda GC)
#   top_hits.tsv                                    P < SUGGESTIVE_P, with case/control freqs
#   freq_cases.afreq, freq_controls.afreq
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_PLINK2}" "${MOD_HTSLIB}" "${MOD_R}"

banner "Step 11 | post-GWAS | ${RUN_TAG}"
RESULT=${GWAS_PREFIX}.PHENO.glm.firth
require_file "${RESULT}"
cd "${GWAS_DIR}"

# -- 1. LocusZoom-ready summary statistics (keeps only ERRCODE == "." rows) ---------
LZ=CBD_GWAS_${RUN_TAG}_hg38_locuszoom.tsv
awk 'BEGIN {OFS = "\t"}
     NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i
               print "#CHROM", "POS", "MarkerID", "REF", "ALT", "A1", "A1_FREQ", "OR", "SE", "P"; next }
     $c["ERRCODE"] == "." && $c["P"] != "NA" { print $c["#CHROM"], $c["POS"], $c["ID"], $c["REF"], $c["ALT"], $c["A1"],
                            $c["A1_FREQ"], $c["OR"], $c["LOG(OR)_SE"], $c["P"] }' \
  "${RESULT}" > "${LZ}"
bgzip -f "${LZ}"
tabix -f -s 1 -b 2 -e 2 -c '#' "${LZ}.gz"
log "summary stats: ${GWAS_DIR}/${LZ}.gz  ($(zcat "${LZ}.gz" | tail -n +2 | wc -l) variants with P)"

# -- 2. allele frequencies in cases and controls --------------------------------------
make_pheno_lists "${GWAS_DIR}"
for GRP in cases controls; do
  plink2 --bfile "${QC_PREFIX}" --keep "${GRP}.txt" --freq \
    --out "freq_${GRP}" --threads "${THREADS}" > /dev/null
done

# -- 3. Manhattan / QQ / lambda / top hits ---------------------------------------------
Rscript "${PIPELINE_ROOT}/scripts/plot_gwas.R" \
  --sumstats "${LZ}.gz" \
  --freq-cases freq_cases.afreq \
  --freq-controls freq_controls.afreq \
  --suggestive "${SUGGESTIVE_P}" \
  --title "CBD GWAS (${RUN_TAG})" \
  --outdir "${GWAS_DIR}"

log "DONE"
