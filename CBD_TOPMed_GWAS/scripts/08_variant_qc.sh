#!/usr/bin/env bash
# =============================================================================
# Step 08 - Variant QC on the merged dataset
#
# Input : WORK_DIR/merged/AllCohorts_merged_raw, SAMPLE_METADATA
# Output: WORK_DIR/qc_<RUN_TAG>/AllCohorts_qc.{bed,bim,fam}
#         WORK_DIR/qc_<RUN_TAG>/qc_summary.tsv     variants remaining after each filter
#
# Filters, in order:
#   a) call rate:  --geno QC_GENO across all samples
#        or, if QC_SEPARATE_GENO is set, call rate >= 1-QC_SEPARATE_GENO in
#        cases AND in controls separately (stricter for unbalanced batches)
#   b) --maf QC_MAF  and  --hwe QC_HWE
#   c) differential missingness (if QC_DIFFMISS is set): remove variants whose
#      missing rate exceeds QC_DIFFMISS in cases OR in controls. This removes
#      variants that are only well covered in some batches and would otherwise
#      create case/control artefacts.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_PLINK2}"

banner "Step 08 | variant QC | ${RUN_TAG}"
require_bfile "${MERGED_PREFIX}"
mkdir -p "${QC_DIR}"
cd "${QC_DIR}"

make_pheno_lists "${QC_DIR}"
SUMMARY=${QC_DIR}/qc_summary.tsv
printf "step\tsamples\tvariants\n" > "${SUMMARY}"
record() { printf "%s\t%s\t%s\n" "$1" "$(count_lines "$2.fam")" "$(count_lines "$2.bim")" >> "${SUMMARY}"; }
record "merged_raw" "${MERGED_PREFIX}"

REMOVE=()
if [[ "${APPLY_SAMPLE_QC}" == "true" ]]; then
  require_file "${SAMPLE_QC_DIR}/samples_to_remove.txt"
  REMOVE=(--remove "${SAMPLE_QC_DIR}/samples_to_remove.txt")
fi

# -- a) call rate --------------------------------------------------------------------
if [[ -n "${QC_SEPARATE_GENO}" ]]; then
  for GRP in cases controls; do
    plink2 --bfile "${MERGED_PREFIX}" --keep "${GRP}.txt" --geno "${QC_SEPARATE_GENO}" \
      --make-just-bim --out "${GRP}_callrate_pass" --threads "${THREADS}"
    cut -f2 "${GRP}_callrate_pass.bim" | sort > "${GRP}_callrate_pass.ids"
  done
  comm -12 cases_callrate_pass.ids controls_callrate_pass.ids > callrate_pass_both.txt
  log "variants passing call rate in cases AND controls: $(count_lines callrate_pass_both.txt)"
  CALLRATE=(--extract callrate_pass_both.txt)
else
  CALLRATE=(--geno "${QC_GENO}")
fi

plink2 --bfile "${MERGED_PREFIX}" "${REMOVE[@]}" "${CALLRATE[@]}" \
  --make-bed --out step_a_callrate --threads "${THREADS}"
record "a_callrate" step_a_callrate

# -- b) MAF + HWE ------------------------------------------------------------------
plink2 --bfile step_a_callrate --maf "${QC_MAF}" --hwe "${QC_HWE}" \
  --make-bed --out step_b_maf_hwe --threads "${THREADS}"
record "b_maf_hwe" step_b_maf_hwe

# -- c) differential missingness -------------------------------------------------------
if [[ -n "${QC_DIFFMISS}" ]]; then
  for GRP in cases controls; do
    plink2 --bfile step_b_maf_hwe --keep "${GRP}.txt" --missing variant-only \
      --out "${GRP}_missing" --threads "${THREADS}"
    awk -v t="${QC_DIFFMISS}" 'NR > 1 && $NF > t {print $2}' "${GRP}_missing.vmiss" > "fail_diffmiss_${GRP}.txt"
    log "variants with missingness > ${QC_DIFFMISS} in ${GRP}: $(count_lines "fail_diffmiss_${GRP}.txt")"
  done
  sort -u fail_diffmiss_cases.txt fail_diffmiss_controls.txt > exclude_diffmiss.txt
  plink2 --bfile step_b_maf_hwe --exclude exclude_diffmiss.txt \
    --make-bed --out "${QC_PREFIX}" --threads "${THREADS}"
  record "c_diffmiss" "${QC_PREFIX}"
else
  for ext in bed bim fam; do mv "step_b_maf_hwe.${ext}" "${QC_PREFIX}.${ext}"; done
fi

rm -f step_a_callrate.{bed,bim,fam} step_b_maf_hwe.{bed,bim,fam}
log "QC summary:"
column -t "${SUMMARY}"
