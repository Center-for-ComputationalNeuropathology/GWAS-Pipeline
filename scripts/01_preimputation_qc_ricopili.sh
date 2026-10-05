#!/usr/bin/env bash
# =============================================================================
# Step 01 - Pre-imputation QC with RICOPILI (preimp_dir)
#
# Runs per cohort on the raw PLINK files. RICOPILI submits its own cluster jobs,
# so run this script from a LOGIN node (not through bsub).
#
# Usage:
#   scripts/01_preimputation_qc_ricopili.sh <cohort>        # one cohort
#   scripts/01_preimputation_qc_ricopili.sh all             # every cohort in the table
#
# Input : RAW_DIR/<raw_dir>/<raw_bfile>.{bed,bim,fam}  (hg19, case/control status in .fam)
# Output: RAW_DIR/<raw_dir>/qc/<qc_prefix>.{bed,bim,fam} + <qc_prefix>.pdf QC report
#
# RICOPILI thresholds (defaults, recorded in qc/*.meta):
#   sample call rate >= 0.98 (--mind 0.02), SNP call rate >= 0.98 (--geno 0.02)
#   |F_het| < 0.2, case/control missingness difference < 0.02,
#   HWE p > 1e-6 in controls / 1e-10 in cases, sex check.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[[ $# -eq 1 ]] || die "usage: $0 <cohort|all>"
if [[ "$1" == "all" ]]; then mapfile -t COHORTS < <(cohort_names); else COHORTS=("$1"); fi

# shellcheck disable=SC1091
source activate "${RICOPILI_CONDA_ENV}" 2>/dev/null || conda activate "${RICOPILI_CONDA_ENV}"
command -v preimp_dir >/dev/null || die "preimp_dir not on PATH - is RICOPILI installed in ${RICOPILI_CONDA_ENV}?"

for COHORT in "${COHORTS[@]}"; do
  RAW=$(cohort_field "${COHORT}" 2)
  BFILE=$(cohort_field "${COHORT}" 3)
  COHORT_DIR=${RAW_DIR}/${RAW}

  banner "Step 01 | RICOPILI preimp_dir | ${COHORT}"
  require_bfile "${COHORT_DIR}/${BFILE}"
  log "samples: $(count_lines "${COHORT_DIR}/${BFILE}.fam")  variants: $(count_lines "${COHORT_DIR}/${BFILE}.bim")"

  cd "${COHORT_DIR}"
  # On the first run preimp_dir writes <dis>.names; set STUDYNAME there (5 chars,
  # e.g. gsa1) and rerun. The studyname becomes part of the qc_prefix.
  preimp_dir \
    --dis "${RICOPILI_DIS}" \
    --popname "${RICOPILI_POP}" \
    --outname "${RICOPILI_OUT}" \
    "${BFILE}.bed"

  log "submitted. Monitor with 'bjobs -w'. When finished check:"
  log "  ${COHORT_DIR}/qc/$(cohort_field "${COHORT}" 4).{bed,bim,fam,pdf}"
done
