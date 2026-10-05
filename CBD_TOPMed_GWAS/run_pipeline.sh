#!/usr/bin/env bash
# =============================================================================
# Submit pipeline steps to IBM LSF, chained with job dependencies.
#
#   ./run_pipeline.sh 02                 # prepare TOPMed VCFs (one array job)
#   ./run_pipeline.sh 03-11              # everything after the TOPMed download
#   ./run_pipeline.sh 08 09 10 11        # re-run QC -> GWAS with new thresholds
#   ./run_pipeline.sh --dry-run 03-11    # print the bsub commands only
#
# Step 01 (RICOPILI) submits its own jobs - run it directly on a login node.
# The TOPMed upload between steps 02 and 03 is manual (see README).
# =============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT}/scripts/lib/common.sh"

DRY_RUN=false
STEPS=()
for a in "$@"; do
  case "${a}" in
    --dry-run) DRY_RUN=true ;;
    -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
    *-*) for ((s = 10#${a%-*}; s <= 10#${a#*-}; s++)); do STEPS+=("$(printf '%02d' "${s}")"); done ;;
    *)   STEPS+=("$(printf '%02d' "$((10#${a}))")") ;;
  esac
done
[[ ${#STEPS[@]} -gt 0 ]] || { sed -n '2,13p' "$0"; exit 1; }

N=$(n_cohorts)
NCHR=$(( N * 22 ))

# step -> "script | array size (0 = single job) | cores | mem MB per core | walltime"
declare -A STEP_SPEC=(
  [02]="02_prepare_topmed_vcfs.sh|${N}|4|6000|6:00"
  [03]="03_unzip_topmed_results.sh|${NCHR}|1|4000|10:00"
  [04]="04_filter_imputed.sh|${NCHR}|1|8000|4:00"
  [05]="05_convert_cohort.sh|${N}|8|4000|24:00"
  [06]="06_merge_cohorts.sh|0|4|16000|24:00"
  [07]="07_sample_qc.sh|0|8|8000|24:00"
  [08]="08_variant_qc.sh|0|8|8000|24:00"
  [09]="09_pca_covariates.sh|0|8|8000|12:00"
  [10]="10_gwas_firth.sh|0|8|8000|48:00"
  [11]="11_post_gwas.sh|0|8|4000|4:00"
)

PREV=""
for STEP in "${STEPS[@]}"; do
  [[ "${STEP}" == "01" ]] && die "step 01 is run directly: scripts/01_preimputation_qc_ricopili.sh <cohort|all>"
  [[ -n "${STEP_SPEC[${STEP}]:-}" ]] || die "unknown step ${STEP}"
  IFS='|' read -r SCRIPT ARRAY CORES MEM WALL <<< "${STEP_SPEC[${STEP}]}"

  NAME="cbdgwas_${STEP}_${RUN_TAG}"
  OUT="${LOG_DIR}/step${STEP}_%J"
  [[ "${ARRAY}" -gt 0 ]] && { NAME="${NAME}[1-${ARRAY}]"; OUT="${OUT}_%I"; }

  CMD=(bsub -P "${LSF_PROJECT}" -q "${LSF_QUEUE}" -J "${NAME}"
       -n "${CORES}" -R "span[hosts=1]" -R "rusage[mem=${MEM}]" -W "${WALL}"
       -oo "${OUT}.out" -eo "${OUT}.err")
  [[ -n "${PREV}" ]] && CMD+=(-w "done(${PREV})")
  CMD+=(env CONFIG_FILE="${CONFIG_FILE}" bash "${ROOT}/scripts/${SCRIPT}")

  if ${DRY_RUN}; then
    printf '%q ' "${CMD[@]}"; echo
    PREV="<job-${STEP}>"
  else
    JOB=$("${CMD[@]}" | sed -n 's/^Job <\([0-9]*\)>.*/\1/p')
    [[ -n "${JOB}" ]] || die "submission of step ${STEP} failed"
    log "step ${STEP} (${SCRIPT}) -> job ${JOB}${PREV:+  (after ${PREV})}"
    PREV=${JOB}
  fi
done

${DRY_RUN} || log "monitor with: bjobs -w | grep cbdgwas ; logs in ${LOG_DIR}"
