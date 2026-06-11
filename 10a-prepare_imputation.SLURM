#!/bin/bash
#SBATCH --job-name=prep_impute
#SBATCH --array=0-22
#SBATCH --time=00:10:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/prep_impute_%A_%a.out
#SBATCH --error=logs/prep_impute_%A_%a.err

# step 13: convert onlyQc PLINK fileset into per-chromosome VCF for TOPMed
# Input:  ${DATA_DIR}/onlyQc.{bed,bim,fam}  (output of step 2)
# Output: per-chromosome VCFs in OUTDIR

module load plink
module load bcftools

set -euo pipefail

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
chr="${CHROMS[$SLURM_ARRAY_TASK_ID]}"

#DATA_DIR=/projects/standard/basu_hbcd/shared/data
DATA_DIR=/home/ood-coffm049/hbcdData
PREFIX="${DATA_DIR}/onlyQc"
OUTDIR="/home/ood-coffm049/hbcdData/data"
mkdir -p $OUTDIR

echo "--- Processing chr${chr} (task ${SLURM_ARRAY_TASK_ID}) ---"

# exclude palindromic snps because they require strand flipping and allele matching
plink2 --bfile "$PREFIX" \
    --chr "$chr" \
    --sort-vars \
    --exclude-palindromic-snps \
    --output-chr chrM \
    --out ${OUTDIR}/temp_raw_chr${chr} \
    --make-pgen \
    --snps-only 'just-acgt'

plink2 --pfile "${OUTDIR}/temp_raw_chr${chr}" \
    --output-chr chrM \
    --export vcf-4.2 bgz \
    --out "${OUTDIR}/temp_raw_chr${chr}"

bcftools +fixref "${OUTDIR}/temp_raw_chr${chr}.vcf.gz" -Oz -o "${OUTDIR}/chr${chr}.vcf.gz" -- \
    -f /scratch.global/GDC/referenceData/hg38.fa \
    -m flip

tabix -f -p vcf "${OUTDIR}/chr${chr}.vcf.gz"

rm "${OUTDIR}/temp_raw_chr${chr}.vcf.gz"
echo "--- Finished chr${chr} ---"
