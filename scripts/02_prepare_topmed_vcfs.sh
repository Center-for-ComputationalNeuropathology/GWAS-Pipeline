#!/usr/bin/env bash
# =============================================================================
# Step 02 - Lift QC'd genotypes to hg38 and write per-chromosome VCFs for the
#           TOPMed Imputation Server.
#
# LSF array: one task per cohort  (run_pipeline.sh 02)
# Manual   : scripts/02_prepare_topmed_vcfs.sh <cohort-index>
#
# Input : RAW_DIR/<raw_dir>/qc/<qc_prefix>.{bed,bim,fam}       (hg19, RICOPILI-QC'd)
# Output: PREP_DIR/<cohort>/<cohort>_chr{1..22}.vcf.gz(.tbi)   <- upload these
#
# What happens:
#   1. autosomal SNPs -> UCSC BED -> liftOver hg19->hg38 (unmapped SNPs dropped)
#   2. rewrite .bim with hg38 CHR/POS, sort, drop duplicate positions
#   3. per chromosome: export VCF, sort, add "chr" prefix, split multi-allelics,
#      check REF against the hg38 FASTA and swap REF/ALT where needed (TOPMed
#      rejects files whose REF does not match the reference)
# =============================================================================
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"
load_modules "${MOD_LIFTOVER}" "${MOD_PLINK2}" "${MOD_BCFTOOLS}" "${MOD_HTSLIB}"

COHORT=$(cohort_by_index "$(task_index "${1:-}")")
RAW=$(cohort_field "${COHORT}" 2)
QCP=$(cohort_field "${COHORT}" 4)
IN=${RAW_DIR}/${RAW}/qc/${QCP}
OUT=${PREP_DIR}/${COHORT}
P=${OUT}/${COHORT}

banner "Step 02 | liftOver + TOPMed VCFs | ${COHORT}"
require_bfile "${IN}"
require_file "${LIFTOVER_CHAIN}" "${REF_FASTA_HG38}" "${PIPELINE_ROOT}/config/chr_rename.txt"
mkdir -p "${OUT}"

# -- 1. liftOver ----------------------------------------------------------------
awk '$1 ~ /^([1-9]|1[0-9]|2[0-2])$/ { print "chr"$1"\t"($4-1)"\t"$4"\t"$2 }' \
  "${IN}.bim" > "${P}_hg19.bed"
log "autosomal input variants: $(count_lines "${P}_hg19.bed")"

liftOver "${P}_hg19.bed" "${LIFTOVER_CHAIN}" "${P}_hg38.bed" "${P}_unmapped.bed"
log "lifted: $(count_lines "${P}_hg38.bed")  unmapped: $(grep -vc '^#' "${P}_unmapped.bed" || true)"

# SNP -> new chr/pos; drop SNPs that moved to non-autosomal contigs
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ { sub("chr", "", $1); print $4"\t"$1"\t"$3 }' \
  "${P}_hg38.bed" > "${P}_snp_lookup.txt"
cut -f1 "${P}_snp_lookup.txt" > "${P}_keep_snps.txt"

# -- 2. rewrite coordinates, sort, deduplicate ------------------------------------
plink2 --bfile "${IN}" --extract "${P}_keep_snps.txt" --make-bed --out "${P}_tmp" --threads 4

awk 'NR==FNR { chr[$1]=$2; pos[$1]=$3; next }
     { print chr[$2]"\t"$2"\t"$3"\t"pos[$2]"\t"$5"\t"$6 }' \
  "${P}_snp_lookup.txt" "${P}_tmp.bim" > "${P}_tmp_hg38.bim"
cp "${P}_tmp.bed" "${P}_tmp_hg38.bed"
cp "${P}_tmp.fam" "${P}_tmp_hg38.fam"

# sort-vars fixes "split chromosome" errors caused by SNPs that jumped chromosomes
plink2 --bfile "${P}_tmp_hg38" --make-pgen sort-vars --out "${P}_tmp_sorted" --threads 4
plink2 --pfile "${P}_tmp_sorted" --rm-dup exclude-mismatch --make-pgen \
  --out "${P}_hg38_dedup" --threads 4
log "after dedup: $(grep -vc '^#' "${P}_hg38_dedup.pvar") variants"

# -- 3. per-chromosome VCFs ---------------------------------------------------------
for CHR in $(seq 1 22); do
  plink2 --pfile "${P}_hg38_dedup" --chr "${CHR}" \
    --export vcf bgz id-paste=iid --out "${P}_chr${CHR}_raw" --threads 4 >/dev/null

  bcftools sort "${P}_chr${CHR}_raw.vcf.gz" -Ou -T "${OUT}/sort_tmp_${CHR}" \
  | bcftools annotate --rename-chrs "${PIPELINE_ROOT}/config/chr_rename.txt" -Ou \
  | bcftools norm -m -any --check-ref ws --fasta-ref "${REF_FASTA_HG38}" \
      -Oz -o "${P}_chr${CHR}.vcf.gz"
  tabix -f -p vcf "${P}_chr${CHR}.vcf.gz"

  rm -f "${P}_chr${CHR}_raw.vcf.gz" "${P}_chr${CHR}_raw.log"
  log "chr${CHR}: $(bcftools index -n "${P}_chr${CHR}.vcf.gz") variants"
done

rm -f "${P}"_tmp*.{bed,bim,fam,log,pgen,pvar,psam}

log "DONE. Upload ${OUT}/${COHORT}_chr{1..22}.vcf.gz to https://imputation.biodatacatalyst.nhlbi.nih.gov"
