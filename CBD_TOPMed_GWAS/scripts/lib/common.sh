#!/usr/bin/env bash
# =============================================================================
# Shared helpers - sourced by every pipeline step.
#   * loads config/config.sh (override with CONFIG_FILE=/path/to/config.sh)
#   * logging / error handling
#   * cohort-table lookups
#   * LSF array-index helpers (steps also run interactively without LSF)
# =============================================================================
set -euo pipefail

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PIPELINE_ROOT

CONFIG_FILE="${CONFIG_FILE:-${PIPELINE_ROOT}/config/config.sh}"
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "ERROR: config not found: ${CONFIG_FILE}" >&2
  echo "       cp ${PIPELINE_ROOT}/config/config.example.sh ${PIPELINE_ROOT}/config/config.sh" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "${CONFIG_FILE}"

# ---------------------------------------------------------------- logging ----
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $*" >&2; }
die()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; exit 1; }

banner() {
  echo "=============================================================================="
  echo " $*"
  echo " host: $(hostname)   start: $(date)"
  echo "=============================================================================="
}

require_file() {
  local f
  for f in "$@"; do
    [[ -s "${f}" ]] || die "required input missing or empty: ${f}"
  done
}

require_bfile() {
  local p
  for p in "$@"; do require_file "${p}.bed" "${p}.bim" "${p}.fam"; done
}

# Load environment modules when available (Lmod/Environment Modules).
load_modules() {
  if type module >/dev/null 2>&1; then
    module load "$@"
  else
    warn "'module' not available - assuming these tools are on PATH: $*"
  fi
}

count_lines() { wc -l < "$1" | tr -d ' '; }

# ---------------------------------------------------------- cohort table ----
# Columns: 1 cohort | 2 raw_dir | 3 raw_bfile | 4 qc_prefix | 5 description
cohort_rows()  { grep -v -E '^[[:space:]]*(#|$)' "${COHORT_TABLE}"; }
cohort_names() { cohort_rows | cut -f1; }
n_cohorts()    { cohort_rows | wc -l | tr -d ' '; }

# cohort_field <cohort> <column-number>
cohort_field() {
  cohort_rows | awk -F'\t' -v c="$1" -v k="$2" '$1==c {print $k; found=1} END {exit !found}' \
    || die "cohort '$1' not found in ${COHORT_TABLE}"
}

# cohort_by_index <1-based index>
cohort_by_index() {
  local name
  name=$(cohort_names | sed -n "${1}p")
  [[ -n "${name}" ]] || die "no cohort at index $1 (table has $(n_cohorts) rows)"
  echo "${name}"
}

# ---------------------------------------------------------- array helpers ----
# Index of this task: LSF array index, else first CLI argument.
task_index() {
  local idx="${LSB_JOBINDEX:-${1:-}}"
  [[ -n "${idx}" && "${idx}" != "0" ]] || die "no task index: run as an LSF array job or pass the index as argument"
  echo "${idx}"
}

# Map a 1-based index over (cohort x chr1..22) to "<cohort> <chr>".
cohort_chr_from_index() {
  local idx=$(( $1 - 1 ))
  echo "$(cohort_by_index $(( idx / 22 + 1 ))) $(( idx % 22 + 1 ))"
}

# --------------------------------------------------------- standard paths ----
FILTERED_VCF_DIR=${WORK_DIR}/filtered_vcf          # step 04
PER_COHORT_DIR=${WORK_DIR}/per_cohort              # step 05
MERGED_DIR=${WORK_DIR}/merged                      # step 06
SAMPLE_QC_DIR=${WORK_DIR}/sample_qc                # step 07
QC_DIR=${WORK_DIR}/qc_${RUN_TAG}                   # steps 08-09
GWAS_DIR=${WORK_DIR}/gwas_${RUN_TAG}               # steps 10-11

MERGED_PREFIX=${MERGED_DIR}/AllCohorts_merged_raw
QC_PREFIX=${QC_DIR}/AllCohorts_qc
COVAR_FILE=${QC_DIR}/covariates.txt
GWAS_PREFIX=${GWAS_DIR}/GWAS_${RUN_TAG}_firth

mkdir -p "${LOG_DIR}"

# Case / control keep-lists (FID IID) derived from SAMPLE_METADATA.
# make_pheno_lists <outdir>  ->  <outdir>/cases.txt, <outdir>/controls.txt
make_pheno_lists() {
  local out="$1"
  require_file "${SAMPLE_METADATA}"
  awk -v out="${out}" '
    NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i
              if (!("FID" in col) || !("IID" in col) || !("PHENO" in col)) {
                print "SAMPLE_METADATA needs FID, IID and PHENO columns" > "/dev/stderr"; exit 1 }
              next }
    $col["PHENO"] == 1 { print $col["FID"], $col["IID"] > (out "/cases.txt") }
    $col["PHENO"] == 0 { print $col["FID"], $col["IID"] > (out "/controls.txt") }
  ' "${SAMPLE_METADATA}"
  log "phenotype lists: $(count_lines "${out}/cases.txt") cases, $(count_lines "${out}/controls.txt") controls"
}
