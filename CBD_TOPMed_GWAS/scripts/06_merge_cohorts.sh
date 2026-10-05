#!/usr/bin/env bash
# =============================================================================
# Step 06 - Merge all cohorts into one PLINK dataset (union of variants)
#
# Input : WORK_DIR/per_cohort/<cohort>_bed.{bed,bim,fam}   (all cohorts)
# Output: WORK_DIR/merged/AllCohorts_merged_raw.{bed,bim,fam}
#         WORK_DIR/merged/duplicate_ids.txt   IIDs present in >1 cohort
#
# Notes
#   * Variants missing from a cohort become missing genotypes for that cohort's
#     samples; they are removed later by the call-rate filters in step 08.
#   * Variants with >2 alleles across cohorts (PLINK writes them to
#     *-merge.missnp) are excluded from every cohort and the merge is retried.
#   * Samples with the same FID/IID in two cohorts are MERGED into one sample
#     (PLINK sets discordant calls to missing). This is reported, not hidden.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_PLINK19}"

banner "Step 06 | merge cohorts"
mkdir -p "${MERGED_DIR}"

mapfile -t COHORTS < <(cohort_names)
for C in "${COHORTS[@]}"; do require_bfile "${PER_COHORT_DIR}/${C}_bed"; done

BASE=${PER_COHORT_DIR}/${COHORTS[0]}_bed
MERGE_LIST=${MERGED_DIR}/merge_list.txt
printf "%s\n" "${COHORTS[@]:1}" | sed "s|^|${PER_COHORT_DIR}/|; s|$|_bed|" > "${MERGE_LIST}"

# -- duplicate sample report ----------------------------------------------------
for C in "${COHORTS[@]}"; do
  awk -v c="${C}" '{print $2"\t"c}' "${PER_COHORT_DIR}/${C}_bed.fam"
done | sort -k1,1 \
  | awk '{n[$1]++; c[$1]=c[$1] (c[$1]?",":"") $2} END {for (i in n) if (n[i]>1) print i"\t"c[i]}' \
  | sort > "${MERGED_DIR}/duplicate_ids.txt"
if [[ -s "${MERGED_DIR}/duplicate_ids.txt" ]]; then
  warn "$(count_lines "${MERGED_DIR}/duplicate_ids.txt") IIDs occur in more than one cohort; they will be merged into one sample."
  cut -f2 "${MERGED_DIR}/duplicate_ids.txt" | sort | uniq -c | sed 's/^/          /'
fi

# -- merge (retry once without multi-allelic conflicts) ----------------------------
set +e
plink --bfile "${BASE}" --merge-list "${MERGE_LIST}" \
  --make-bed --out "${MERGED_PREFIX}" --threads "${THREADS}" > "${MERGED_DIR}/first_merge_attempt.log" 2>&1
STATUS=$?
set -e

if [[ ${STATUS} -ne 0 ]]; then
  MISSNP=${MERGED_PREFIX}-merge.missnp
  [[ -s "${MISSNP}" ]] || die "merge failed for a reason other than allele conflicts - see ${MERGED_DIR}/first_merge_attempt.log"
  warn "$(count_lines "${MISSNP}") variants with >2 alleles across cohorts - excluding and retrying"

  : > "${MERGED_DIR}/merge_list_excl.txt"
  for C in "${COHORTS[@]}"; do
    plink --bfile "${PER_COHORT_DIR}/${C}_bed" --exclude "${MISSNP}" \
      --make-bed --out "${MERGED_DIR}/${C}_excl" --threads "${THREADS}" > /dev/null
    [[ "${C}" == "${COHORTS[0]}" ]] || echo "${MERGED_DIR}/${C}_excl" >> "${MERGED_DIR}/merge_list_excl.txt"
  done
  plink --bfile "${MERGED_DIR}/${COHORTS[0]}_excl" --merge-list "${MERGED_DIR}/merge_list_excl.txt" \
    --make-bed --out "${MERGED_PREFIX}" --threads "${THREADS}"
  rm -f "${MERGED_DIR}"/*_excl.{bed,bim,fam,log,nosex}
fi

require_bfile "${MERGED_PREFIX}"
log "merged: $(count_lines "${MERGED_PREFIX}.fam") samples, $(count_lines "${MERGED_PREFIX}.bim") variants"
