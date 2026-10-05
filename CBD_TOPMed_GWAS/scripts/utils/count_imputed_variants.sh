#!/usr/bin/env bash
# =============================================================================
# Variant counts per cohort before and after the step-04 imputation filter.
#
#   scripts/utils/count_imputed_variants.sh > imputation_filter_summary.tsv
#
# Fast when the VCFs are indexed; un-indexed TOPMed dose files are streamed
# (submit with bsub if counting all cohorts).
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"
load_modules "${MOD_BCFTOOLS}" >/dev/null 2>&1

# count records: use the index when present, otherwise stream the file
n_records() {
  bcftools index -n "$1" 2>/dev/null || bcftools view -H "$1" | wc -l
}

printf "cohort\tbefore_filter\tafter_filter\tpct_retained\n"
for COHORT in $(cohort_names); do
  BEFORE=0; AFTER=0
  for CHR in $(seq 1 22); do
    IN=${IMPUTED_DIR}/${COHORT}/chr${CHR}.dose.vcf.gz
    OUT=${FILTERED_VCF_DIR}/${COHORT}_chr${CHR}_filtered.vcf.gz
    [[ -f "${IN}" ]]  && BEFORE=$(( BEFORE + $(n_records "${IN}") ))
    [[ -f "${OUT}" ]] && AFTER=$(( AFTER + $(n_records "${OUT}") ))
  done
  awk -v c="${COHORT}" -v b="${BEFORE}" -v a="${AFTER}" \
    'BEGIN { printf "%s\t%d\t%d\t%.1f\n", c, b, a, (b > 0 ? 100 * a / b : 0) }'
done
