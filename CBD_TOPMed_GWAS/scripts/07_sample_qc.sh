#!/usr/bin/env bash
# =============================================================================
# Step 07 - Sample QC report: missingness and relatedness
#
# Input : WORK_DIR/merged/AllCohorts_merged_raw
# Output: WORK_DIR/sample_qc/
#           sample_missingness.smiss      per-sample missing rate
#           fail_mind.txt                 samples with missingness > QC_MIND
#           related_pairs.genome          pairs with PI_HAT > IBD_PIHAT
#           samples_to_remove.txt         (only used if APPLY_SAMPLE_QC=true)
#
# Why removal is opt-in: cohorts genotyped on sparser arrays carry more
# missingness after the union merge. In the CBD data all 69 GSA-array cases
# exceed --mind 0.1, so a blanket --mind would delete a whole case batch.
# Batch-driven missingness is handled at the VARIANT level instead (step 08,
# QC_DIFFMISS). Always review these reports before the GWAS.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_PLINK19}" "${MOD_PLINK2}"

banner "Step 07 | sample QC report"
require_bfile "${MERGED_PREFIX}"
mkdir -p "${SAMPLE_QC_DIR}"
O=${SAMPLE_QC_DIR}

# -- 1. sample missingness ---------------------------------------------------------
plink2 --bfile "${MERGED_PREFIX}" --missing sample-only \
  --out "${O}/sample_missingness" --threads "${THREADS}"
awk -v t="${QC_MIND}" 'NR > 1 && $NF > t {print $1"\t"$2"\t"$NF}' \
  "${O}/sample_missingness.smiss" > "${O}/fail_mind.txt"
log "samples with missingness > ${QC_MIND}: $(count_lines "${O}/fail_mind.txt")"

# -- 2. relatedness on common, LD-pruned variants ----------------------------------------
plink2 --bfile "${MERGED_PREFIX}" --geno 0.02 --maf 0.05 --hwe 1e-6 \
  --indep-pairwise ${PRUNE_PARAMS} --out "${O}/relatedness_prune" --threads "${THREADS}"
plink --bfile "${MERGED_PREFIX}" --extract "${O}/relatedness_prune.prune.in" \
  --genome --min "${IBD_PIHAT}" --out "${O}/related_pairs" --threads "${THREADS}"
N_PAIRS=$(( $(count_lines "${O}/related_pairs.genome") - 1 ))
log "pairs with PI_HAT > ${IBD_PIHAT}: ${N_PAIRS}  (PI_HAT ~1 = duplicate/MZ twin, ~0.5 = 1st degree)"

# -- 3. removal list ---------------------------------------------------------------
# One sample per related pair is dropped (the second of the pair) plus all
# --mind failures.
{
  cut -f1,2 "${O}/fail_mind.txt"
  awk 'NR > 1 {print $3"\t"$4}' "${O}/related_pairs.genome"
} | sort -u > "${O}/samples_to_remove.txt"

if [[ "${APPLY_SAMPLE_QC}" == "true" ]]; then
  log "APPLY_SAMPLE_QC=true -> $(count_lines "${O}/samples_to_remove.txt") samples will be removed in step 08"
else
  log "APPLY_SAMPLE_QC=false -> report only; no samples removed"
fi
