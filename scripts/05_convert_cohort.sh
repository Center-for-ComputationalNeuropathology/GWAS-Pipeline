#!/usr/bin/env bash
# =============================================================================
# Step 05 - Per cohort: concatenate chromosomes and convert to PLINK
#
# LSF array: one task per cohort
# Manual   : scripts/05_convert_cohort.sh <cohort-index>
#
# Input : WORK_DIR/filtered_vcf/<cohort>_chr{1..22}_filtered.vcf.gz
# Output: WORK_DIR/per_cohort/<cohort>.{pgen,pvar,psam}   dosages (DS field)
#         WORK_DIR/per_cohort/<cohort>_bed.{bed,bim,fam}  hard calls, used for merging
#
# Variant IDs are rewritten to CHR:POS:REF:ALT so the same variant has the same
# ID in every cohort (TOPMed IDs are not unique across multi-allelic sites).
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_BCFTOOLS}" "${MOD_HTSLIB}" "${MOD_PLINK2}"

COHORT=$(cohort_by_index "$(task_index "${1:-}")")
P=${PER_COHORT_DIR}/${COHORT}

banner "Step 05 | concatenate + convert | ${COHORT}"
mkdir -p "${PER_COHORT_DIR}"

# -- 1. concatenate chr1..22 in numeric order -----------------------------------
LIST=${P}_vcf_list.txt
: > "${LIST}"
for CHR in $(seq 1 22); do
  F=${FILTERED_VCF_DIR}/${COHORT}_chr${CHR}_filtered.vcf.gz
  require_file "${F}"
  echo "${F}" >> "${LIST}"
done

bcftools concat -f "${LIST}" -Oz -o "${P}_filtered.vcf.gz" --threads "${THREADS}"
bcftools index -f "${P}_filtered.vcf.gz"
log "concatenated: $(bcftools index -n "${P}_filtered.vcf.gz") variants"

# -- 2. VCF -> pgen (dosages) ----------------------------------------------------
plink2 --vcf "${P}_filtered.vcf.gz" dosage=DS \
  --make-pgen --out "${P}_tmp" --threads "${THREADS}"

# -- 3. unique CHR:POS:REF:ALT IDs ---------------------------------------------------
plink2 --pfile "${P}_tmp" \
  --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len "${MAX_ALLELE_LEN}" \
  --make-pgen --out "${P}" --threads "${THREADS}"

# -- 4. hard-call BED for the PLINK 1.9 merge -----------------------------------------
plink2 --pfile "${P}" --make-bed --out "${P}_bed" --threads "${THREADS}"

rm -f "${P}_tmp".{pgen,pvar,psam,log}
log "samples: $(count_lines "${P}_bed.fam")  variants: $(count_lines "${P}_bed.bim")"
