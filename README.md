[![Built With](https://img.shields.io/badge/Built%20With-Nextflow-brightgreen.svg)](https://www.nextflow.io/)
[![Compatibility](https://img.shields.io/badge/Compatibility-Linux+%2F+OSX-blue.svg)]()
[![GitHub Issues](https://img.shields.io/github/issues/your-username/your-repo-name.svg)](https://github.com/your-username/your-repo-name/issues)
[![GitHub Open](https://img.shields.io/badge/open-1-yellow.svg)]()

# GWAS Pipeline

A comprehensive pipeline for performing Genome-Wide Association Studies (GWAS), from quality control through association analysis and visualization.

Date: February 2025

## Overview

This pipeline provides a standardized workflow for GWAS analysis, including quality control, imputation, ancestry prediction, and association testing. The pipeline is designed to work with genotype data and includes tools for handling population stratification and relatedness.

## 🆕 TOPMed-imputation pipeline (multi-array case/control): `CBD_TOPMed_GWAS/`

A configurable, LSF-ready pipeline used for the CBD GWAS: RICOPILI pre-imputation QC → hg38 liftover →
**TOPMed Imputation Server** → per-cohort filtering → multi-array merge → sample/variant QC
(incl. differential missingness) → PCA → PLINK2 Firth GWAS → LocusZoom-ready summary statistics,
Manhattan/QQ plots. See **[CBD_TOPMed_GWAS/README.md](CBD_TOPMed_GWAS/README.md)** for step-by-step
instructions with examples. The sections below describe the original (1000G / Minimac4) workflow.

## Pipeline Structure

1. Sample and SNP QC: Data Cleaning
2. Imputation
3. Ancestry prediction
4. GWAS
5. Visualize results on Locus Zoom
<!-- 6. Run case-case GWAS -->
<!-- 7. Calculate Polygenic Risk Scores / Pathway scores -->

## Installation / Prerequisites

* PLINK (v1.90b6.21 or newer)
* FlashPCA
* Minimac4
* vcftools
* bcftools
* R
* Somalier

## Setup

### Clone this repository
```bash
git clone https://github.com/Shrishtee-kandoi/GWAS-Pipeline.git
cd GWAS-Pipeline
```

### Schematic diagram of how GWAS works!
![image](https://github.com/Shrishtee-kandoi/GWAS_Pipeline_CraryLab/assets/98359418/4d515baa-2f33-4be3-ad51-fbf7ea45e7f2)

# 1. Quality Control (Sample and SNP level): QC_GWAS/

* Missingness filtering (SNPs and individuals)
* Sex discrepancy checks
* Minor allele frequency filtering
* Hardy-Weinberg equilibrium testing
* Heterozygosity checks
* Relatedness analysis
* Principal Component Analysis

## 2. Imputation: Imputation_scripts/

* Data preparation
* Chromosome separation
* VCF conversion
* Michigan Imputation Server or Minimac4

### === CONFIGURATION ===

```bash
INPUT_FILE="your_clean_data"     # Base name of your input PLINK files after QC
OUTPUT_DIR="imputation_ready"    # Directory to store processed files

CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"  # Chromosomes to process
```

### === LOAD MODULES ===

#### module load plink
#### module load vcftools
#### module load bcftools
#### module load minimac4

### === CREATE DIRECTORIES ===
mkdir -p ${OUTPUT_DIR}/chrom

### === PROCESS DATA ===
echo "Starting imputation preparation for ${INPUT_FILE}"

## 1. Generate frequency file (useful for QC and post-imputation filtering)
* Step 1: Generating frequency file

```bash
plink --bfile ${INPUT_FILE} \
      --keep-allele-order \
      --freq \
      --out ${OUTPUT_DIR}/${INPUT_FILE}_freq \
      --allow-no-sex

## 2. Split data by chromosome
* Step 2: Splitting data by chromosome

for chr in ${CHROMOSOMES}; do
    echo "  Processing chromosome ${chr}"
    plink --bfile ${INPUT_FILE} \
          --chr ${chr} \
          --make-bed \
          --out ${OUTPUT_DIR}/chrom/${INPUT_FILE}_chr${chr}
done
```

## 3. Convert to VCF format
* Step 3: Converting to VCF format

```bash
for chr in ${CHROMOSOMES}; do
    # Skip sex chromosomes if needed
    if [[ ${chr} == "X" ]] || [[ ${chr} == "Y" ]]; then
        echo "  Skipping chromosome ${chr} for this example (modify if needed)"
        continue
    fi
    
    echo "Converting chromosome ${chr} to VCF"
    plink --bfile ${OUTPUT_DIR}/chrom/${INPUT_FILE}_chr${chr} \
          --recode vcf \
          --chr ${chr} \
          --out ${OUTPUT_DIR}/chrom/${INPUT_FILE}_chr${chr} \
          --allow-no-sex
done
```

# 4. Sort and compress VCF files
* Step 4: Sorting and compressing VCF files

```bash
for chr in ${CHROMOSOMES}; do
    # Skip sex chromosomes if needed
    if [[ ${chr} == "X" ]] || [[ ${chr} == "Y" ]]; then
        continue
    fi
    
    echo "  Processing chromosome ${chr}"
    vcf_file="${OUTPUT_DIR}/chrom/${INPUT_FILE}_chr${chr}.vcf"
    
    if [ -f "$vcf_file" ]; then
        vcf-sort ${vcf_file} | bgzip -c > ${OUTPUT_DIR}/chrom/${INPUT_FILE}_chr${chr}.vcf.gz
    else
        echo "  Warning: VCF file for chromosome ${chr} not found"
    fi
done
```

# 5. Index VCF files
* Step 5: Indexing VCF files
```bash
for chr in ${CHROMOSOMES}; do
    # Skip sex chromosomes if needed
    if [[ ${chr} == "X" ]] || [[ ${chr} == "Y" ]]; then
        continue
    fi
    
    gz_file="${OUTPUT_DIR}/chrom/${INPUT_FILE}_chr${chr}.vcf.gz"
    
    if [ -f "$gz_file" ]; then
        echo "  Indexing chromosome ${chr}"
        tabix -p vcf ${gz_file}
    else
        echo "  Warning: Compressed VCF file for chromosome ${chr} not found"
    fi
done
```

### Next steps

a. Michigan Imputation Server

https://imputationserver.sph.umich.edu/index.html

* Build: hg19
* Reference Panel: 1000g-phase-3-v5 (hg19)
* Population: mixed
* Phasing: eagle
* Mode: imputation

b. Minimac4: https://github.com/statgen/Minimac4

* Reference panel: https://share.sph.umich.edu/minimac4/panels/
    * g1k_p3_msav_files_with_estimates.tar.gz (hg19)
    * Note: To use hg38, generate one on your own using a phased 1000g call set with `minimac4 --compress-reference`

```bash
ml minimac4
i=${1}

dir=/sc/arion/projects/tauomics/Shrishtee/PART_DATA_SOW_KAT_CHOP/KAT
msav_files=/sc/arion/projects/tauomics/Shrishtee/Imputation/g1k_p3_msav_files_with_estimates

minimac4 ${msav_files}/1000g_phase3_v5.chr$i.with_parameter_estimates.msav ${dir}/chrom/pre_impute_PART_KatQC_Hg19_separate.$i.vcf.gz -o ${dir}/Imputed/imputed_PART_KatQC_Hg19_separate.$i.vcf.gz
```
# 3. Ancestry prediction

This project uses [somalier](https://github.com/brentp/somalier) for relatedness analysis of sequencing data.

# 4. GWAS

Now, we take the Imputed files and run an association analysis using plink.
* Plink version: plink/1.90b6.21

```bash
i=${1}

ml plink

cd WORKING_DIRECTORY

plink \
--bfile PLINK_FILE$i \
--logistic \
--ci 0.95 \
--pheno Covariate_file \
--pheno-name status \
--covar Covariate_file \
--covar-name PC01,PC02,PC03,sex,express,gsa,ucla,mssm \
--keep-allele-order \
--allow-no-sex \
--adjust \
--out Output_file_chr$i

awk 'NR==1 || $5 == "ADD" {print}' Output_file_chr$i > Output_file_chr$i.ADD


cat  Output_file_chr*.ADD | head -n1 >  Output_file_all_chrom.ADD_gwas
cat Output_file_chr*.ADD | grep -v "SE" >> Output_file_all_chrom.ADD_gwas
cat Output_file_all_chrom.ADD_gwas | grep -v "NA" > Output_file_all_chrom.ADD_gwas_without_NA

```

#### Prep the association analysis results to upload to Locus Zoom
```Rscript Association_analysis/Locus_zoom_prep_and_Manhattan_plot.R ```



