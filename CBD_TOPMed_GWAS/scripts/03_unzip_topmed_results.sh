#!/usr/bin/env bash
# =============================================================================
# Step 03 - Extract the password-protected TOPMed result archives
#
# LSF array: one task per (cohort x chromosome) = n_cohorts * 22 tasks
# Manual   : scripts/03_unzip_topmed_results.sh <task-index>
#
# Input : IMPUTED_DIR/<cohort>/chr_<N>.zip           (downloaded from TOPMed)
#         config/topmed_passwords.tsv                 (<cohort>\t<password>)
# Output: IMPUTED_DIR/<cohort>/chr<N>.dose.vcf.gz, chr<N>.info.gz, ...
#
# Passwords are read inside the job from a git-ignored file - they never appear
# on the bsub command line, in `bjobs` output or in the repository.
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_P7ZIP}"

read -r COHORT CHR <<< "$(cohort_chr_from_index "$(task_index "${1:-}")")"
DIR=${IMPUTED_DIR}/${COHORT}
ZIP=${DIR}/chr_${CHR}.zip

banner "Step 03 | unzip TOPMed results | ${COHORT} chr${CHR}"
require_file "${ZIP}" "${TOPMED_PASSWORDS}"

if [[ -s "${DIR}/chr${CHR}.dose.vcf.gz" ]]; then
  log "already extracted - skipping"
  exit 0
fi

PASSWORD=$(awk -F'\t' -v c="${COHORT}" '$1==c {print $2}' "${TOPMED_PASSWORDS}")
[[ -n "${PASSWORD}" ]] || die "no password for ${COHORT} in ${TOPMED_PASSWORDS}"

7zz x -p"${PASSWORD}" -o"${DIR}" -y "${ZIP}" >/dev/null
require_file "${DIR}/chr${CHR}.dose.vcf.gz"
log "extracted ${DIR}/chr${CHR}.dose.vcf.gz"
