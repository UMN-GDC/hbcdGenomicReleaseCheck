#!/bin/bash
#SBATCH --job-name=prep_impute
#SBATCH --time=00:10:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/prep_imputeX.out
#SBATCH --error=logs/prep_imputeX.err

# step 14: convert onlyQc chrX to VCF for TOPMed (separate due to X-chr handling)

module load plink
module load bcftools

set -euo pipefail

DATA_DIR=/home/ood-coffm049/hbcdData
PREFIX="${DATA_DIR}/onlyQc"
OUTDIR="/home/ood-coffm049/hbcdData/data"
mkdir -p $OUTDIR

echo "--- Processing chrX ---"

plink2 --bfile "$PREFIX" \
    --chr X \
    --sort-vars \
    --exclude-palindromic-snps \
    --set-invalid-haploid-missing \
    --output-chr chrM \
    --out ${OUTDIR}/temp_raw_chrX \
    --make-pgen \
    --snps-only 'just-acgt'

plink2 --pfile "${OUTDIR}/temp_raw_chrX" \
    --output-chr chrM \
    --export vcf-4.2 bgz \
    --out "${OUTDIR}/temp_raw_chrX"

bcftools +fixref "${OUTDIR}/temp_raw_chrX.vcf.gz" -Oz -o "${OUTDIR}/chrX.vcf.gz" -- \
    -f /scratch.global/GDC/referenceData/hg38.fa \
    -m flip

tabix -f -p vcf "${OUTDIR}/chrX.vcf.gz"

rm "${OUTDIR}/temp_raw_chrX.vcf.gz"
echo "--- Finished chrX ---"
