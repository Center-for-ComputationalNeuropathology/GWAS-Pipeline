#!/usr/bin/env bash
# =============================================================================
# Step 04 - Keep well-imputed, non-rare variants (per cohort, per chromosome)
#
# LSF array: one task per (cohort x chromosome) = n_cohorts * 22 tasks
# Manual   : scripts/04_filter_imputed.sh <task-index>
#
# Filter : INFO/R2 >= IMPUTE_R2  and  INFO/MAF >= IMPUTE_MAF   (minimac4 fields)
# Input  : IMPUTED_DIR/<cohort>/chr<N>.dose.vcf.gz
# Output : WORK_DIR/filtered_vcf/<cohort>_chr<N>_filtered.vcf.gz(.csi)
#
# Every cohort is filtered with the SAME thresholds so that imputation quality
# is comparable between case and control batches.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_BCFTOOLS}" "${MOD_HTSLIB}"

read -r COHORT CHR <<< "$(cohort_chr_from_index "$(task_index "${1:-}")")"
IN=${IMPUTED_DIR}/${COHORT}/chr${CHR}.dose.vcf.gz
OUT=${FILTERED_VCF_DIR}/${COHORT}_chr${CHR}_filtered.vcf.gz

banner "Step 04 | filter R2>=${IMPUTE_R2} MAF>=${IMPUTE_MAF} | ${COHORT} chr${CHR}"
require_file "${IN}"
mkdir -p "${FILTERED_VCF_DIR}"

bcftools view -i "MAF>=${IMPUTE_MAF} & R2>=${IMPUTE_R2}" -Oz -o "${OUT}" "${IN}"
bcftools index -f "${OUT}"

log "kept $(bcftools index -n "${OUT}") variants -> ${OUT}"
