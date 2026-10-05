# =============================================================================
# CBD GWAS pipeline - configuration
#
#   cp config/config.example.sh config/config.sh      # config.sh is git-ignored
#   cp config/cohorts.example.tsv config/cohorts.tsv
#
# Every step sources this file. Values below are the ones used for the
# Mount Sinai (Minerva) CBD analysis and serve as a worked example.
# =============================================================================

# -----------------------------------------------------------------------------
# Paths  (data lives OUTSIDE the repository - never commit genotypes or IDs)
# -----------------------------------------------------------------------------
PROJECT_DIR=/sc/arion/projects/tauomics/Shrishtee/CBD_GWAS_ricopili

COHORT_TABLE=${PIPELINE_ROOT}/config/cohorts.tsv     # one row per genotyping batch
RAW_DIR=${PROJECT_DIR}/cohorts                       # raw PLINK files + RICOPILI qc/ output
PREP_DIR=${PROJECT_DIR}/imputation_prep              # VCFs uploaded to TOPMed
IMPUTED_DIR=${PROJECT_DIR}/Imputed_files_TopMed      # TOPMed downloads (one sub-dir per cohort)
WORK_DIR=${PROJECT_DIR}/pipeline_output              # everything produced after imputation

# Sample metadata: tab/space-delimited with header containing at least
#   FID IID Sex PHENO   (+ any covariates you list in GWAS_COVARS, e.g. Age CHIP)
#   PHENO: 1 = case, 0 = control, NA = exclude.  Sex: 1 = male, 2 = female.
SAMPLE_METADATA=${PROJECT_DIR}/Imputed_files_TopMed/Filtered_merged_files/FINAL_covariates_regenie_v2.txt

# TOPMed result passwords: "<cohort>\t<password>" per line. Git-ignored. chmod 600.
TOPMED_PASSWORDS=${PIPELINE_ROOT}/config/topmed_passwords.tsv

# -----------------------------------------------------------------------------
# Reference files
# -----------------------------------------------------------------------------
LIFTOVER_CHAIN=${HOME}/liftover/hg19ToHg38.over.chain.gz   # UCSC chain, hg19 -> hg38
REF_FASTA_HG38=/sc/arion/projects/tauomics/Shrishtee/Reference_files/hg38.fa   # UCSC-style "chr1" names

# -----------------------------------------------------------------------------
# Cluster (IBM LSF) and software modules
# -----------------------------------------------------------------------------
LSF_PROJECT=acc_tauomics
LSF_QUEUE=premium
LOG_DIR=${WORK_DIR}/logs
THREADS=8

MOD_PLINK19=plink/1.90b6.21
MOD_PLINK2=plink2/v2.00a5.14
MOD_BCFTOOLS=bcftools/1.22
MOD_HTSLIB=htslib
MOD_LIFTOVER=liftover/24-Jan-2025
MOD_R=R/4.4.1
MOD_P7ZIP=p7zip

# RICOPILI (pre-imputation QC) conda environment
RICOPILI_CONDA_ENV=rp_env

# -----------------------------------------------------------------------------
# Step 01 - RICOPILI pre-imputation QC
#   Output prefix becomes <DIS>_<studyname>_<POP>_<OUT>-qc1  (e.g. cbd_gsa1_eur_sk-qc1)
# -----------------------------------------------------------------------------
RICOPILI_DIS=cbd
RICOPILI_POP=eur
RICOPILI_OUT=sk

# -----------------------------------------------------------------------------
# Step 04 - post-imputation filter (applied per cohort, per chromosome)
# -----------------------------------------------------------------------------
IMPUTE_R2=0.8          # minimac R2 (imputation quality)
IMPUTE_MAF=0.01        # within-cohort minor allele frequency

# Step 05 - variant IDs are rewritten to CHR:POS:REF:ALT; plink2 stops with an
#           error if an allele is longer than this (raise it and rerun)
MAX_ALLELE_LEN=100

# -----------------------------------------------------------------------------
# Step 07 - sample QC (report always produced; removal is opt-in)
# -----------------------------------------------------------------------------
APPLY_SAMPLE_QC=false  # true -> drop samples failing --mind and one of each related pair
QC_MIND=0.1            # sample missingness threshold
IBD_PIHAT=0.2          # report pairs with PI_HAT above this (0.2 ~ 2nd-degree relatives)

# -----------------------------------------------------------------------------
# Step 08 - variant QC on the merged dataset
# -----------------------------------------------------------------------------
QC_GENO=0.02           # max variant missingness across all samples
QC_SEPARATE_GENO=""    # e.g. 0.01 -> require call rate >= 99% in cases AND controls
                       #   separately (replaces QC_GENO). Leave empty to disable.
QC_MAF=0.02            # min MAF in the merged dataset
QC_HWE=1e-6            # HWE exact test p-value threshold
QC_DIFFMISS=0.05       # drop variants with missingness > this in cases OR controls.
                       #   Leave empty to disable.

# Tag used to name QC/PCA/GWAS output folders. Change it whenever you change the
# thresholds above so runs never overwrite each other.
RUN_TAG=geno02_maf02_hwe1e6_dm

# -----------------------------------------------------------------------------
# Step 09 - PCA
# -----------------------------------------------------------------------------
PRUNE_PARAMS="200 50 0.25"   # --indep-pairwise window step r2
N_PCS=10

# -----------------------------------------------------------------------------
# Step 10 - association (PLINK2 Firth logistic regression)
# -----------------------------------------------------------------------------
GWAS_COVARS="Sex,PC1,PC2,PC3,PC4,PC5"

# -----------------------------------------------------------------------------
# Step 11 - post-GWAS
# -----------------------------------------------------------------------------
SUGGESTIVE_P=1e-5
