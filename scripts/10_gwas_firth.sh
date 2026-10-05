#!/usr/bin/env bash
# =============================================================================
# Step 10 - Case/control association: PLINK2 Firth logistic regression
#
# Input : WORK_DIR/qc_<RUN_TAG>/AllCohorts_qc, WORK_DIR/qc_<RUN_TAG>/covariates.txt
# Output: WORK_DIR/gwas_<RUN_TAG>/GWAS_<RUN_TAG>_firth.PHENO.glm.firth
#
# Model : logit(P(case)) = b0 + b1*genotype + covariates (GWAS_COVARS)
#   * Firth-penalised regression for every variant: unbiased with unbalanced
#     designs (here ~1 case : 22 controls) and with rare alleles.
#   * --1                           PHENO coded 1 = case, 0 = control
#   * --covar-variance-standardize  rescales covariates (avoids convergence problems)
#   * Samples with PHENO = NA or absent from the covariate file are excluded.
#   * Rows with ERRCODE != "." (FIRTH_CONVERGE_FAIL, UNFINISHED) are unreliable
#     and are dropped from the summary statistics in step 11.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_PLINK2}"

banner "Step 10 | Firth GWAS | ${RUN_TAG}"
require_bfile "${QC_PREFIX}"
require_file "${COVAR_FILE}"
mkdir -p "${GWAS_DIR}"

log "covariates: ${GWAS_COVARS}"
plink2 \
  --bfile "${QC_PREFIX}" \
  --pheno "${COVAR_FILE}" --pheno-name PHENO --1 \
  --covar "${COVAR_FILE}" --covar-name "${GWAS_COVARS}" \
  --covar-variance-standardize \
  --glm firth hide-covar \
  --out "${GWAS_PREFIX}" \
  --threads "${THREADS}"

RESULT=${GWAS_PREFIX}.PHENO.glm.firth
require_file "${RESULT}"
grep -E "binary phenotype loaded|covariates loaded" "${GWAS_PREFIX}.log" || true
log "tested variants: $(( $(count_lines "${RESULT}") - 1 ))  -> ${RESULT}"
