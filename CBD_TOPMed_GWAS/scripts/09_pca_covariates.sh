#!/usr/bin/env bash
# =============================================================================
# Step 09 - Principal components and covariate file
#
# Input : WORK_DIR/qc_<RUN_TAG>/AllCohorts_qc, SAMPLE_METADATA
# Output: WORK_DIR/qc_<RUN_TAG>/pca.eigenvec / pca.eigenval
#         WORK_DIR/qc_<RUN_TAG>/covariates.txt   FID IID <metadata cols> PC1..PCn PHENO
#         WORK_DIR/qc_<RUN_TAG>/pca_plots.pdf
#
# PCs are computed on the final QC'd variant set (LD-pruned) so that they
# capture both ancestry and residual array/batch structure.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_PLINK19}" "${MOD_PLINK2}" "${MOD_R}"

banner "Step 09 | PCA + covariates | ${RUN_TAG}"
require_bfile "${QC_PREFIX}"
require_file "${SAMPLE_METADATA}"
cd "${QC_DIR}"

# -- LD pruning ---------------------------------------------------------------------
plink --bfile "${QC_PREFIX}" --indep-pairwise ${PRUNE_PARAMS} \
  --out pca_prune --threads "${THREADS}"
plink --bfile "${QC_PREFIX}" --extract pca_prune.prune.in \
  --make-bed --out pruned_for_pca --threads "${THREADS}"
log "pruned variant set: $(count_lines pruned_for_pca.bim)"

# -- PCA ----------------------------------------------------------------------------
plink2 --bfile pruned_for_pca --pca "${N_PCS}" --out pca --threads "${THREADS}"

# -- covariates ---------------------------------------------------------------------
Rscript "${PIPELINE_ROOT}/scripts/build_covariates.R" \
  --metadata "${SAMPLE_METADATA}" \
  --eigenvec pca.eigenvec \
  --n-pcs "${N_PCS}" \
  --out "${COVAR_FILE}" \
  --plot pca_plots.pdf

rm -f pruned_for_pca.{bed,bim,fam}
log "covariates: ${COVAR_FILE}"
