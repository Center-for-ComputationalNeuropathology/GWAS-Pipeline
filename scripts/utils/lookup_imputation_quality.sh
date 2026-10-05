#!/usr/bin/env bash
# =============================================================================
# Imputation quality (MAF, R2, typed/imputed) of chosen variants in every cohort.
# Useful to check that a top hit is well imputed in BOTH case and control batches.
#
#   scripts/utils/lookup_imputation_quality.sh snps.tsv > snp_quality.tsv
#
# snps.tsv: tab-separated, no header:  <label> <chr> <pos_hg38>
#   rs242559    17  45948522
#   rs13147207  4   52058612
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

[[ $# -eq 1 ]] || die "usage: $0 snps.tsv"
require_file "$1"

printf "label\tcohort\tchr\tpos\tref\talt\tMAF\tR2\tstatus\n"
while IFS=$'\t' read -r LABEL CHR POS; do
  [[ -z "${LABEL}" || "${LABEL}" == \#* ]] && continue
  for COHORT in $(cohort_names); do
    # chr<N>.info.gz is the sites-only TOPMed file (no index needed; streamed,
    # stops at the first position past the target; "|| true" absorbs SIGPIPE)
    INFO=${IMPUTED_DIR}/${COHORT}/chr${CHR}.info.gz
    HIT=""
    [[ -f "${INFO}" ]] && HIT=$(zcat "${INFO}" | awk -F'\t' -v p="${POS}" '
        !/^#/ && $2 == p {
          maf = r2 = "NA"; typed = ($8 ~ /(^|;)TYPED(;|$)/) ? "typed" : "imputed"
          n = split($8, kv, ";")
          for (i = 1; i <= n; i++) { split(kv[i], x, "="); if (x[1] == "MAF") maf = x[2]; if (x[1] == "R2") r2 = x[2] }
          print $4"\t"$5"\t"maf"\t"r2"\t"typed
        }
        !/^#/ && $2 > p { exit }' || true)
    if [[ -z "${HIT}" ]]; then
      printf "%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tnot_present\n" "${LABEL}" "${COHORT}" "${CHR}" "${POS}"
      continue
    fi
    while IFS=$'\t' read -r REF ALT MAF R2 STATUS; do
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${LABEL}" "${COHORT}" "${CHR}" "${POS}" "${REF}" "${ALT}" "${MAF}" "${R2}" "${STATUS}"
    done <<< "${HIT}"
  done
done < "$1"
