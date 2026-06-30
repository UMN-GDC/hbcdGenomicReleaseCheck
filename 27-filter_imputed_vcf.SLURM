#!/bin/bash
#SBATCH --job-name=filter_vcf
#SBATCH --array=0-23
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/filter_vcf_%A_%a.out
#SBATCH --error=logs/filter_vcf_%A_%a.err

# step 27: filter imputed VCFs to only include release subjects
# Uses keep_list.txt from step 25 as the sample whitelist.
# Run AFTER step 25 (keep_list.txt must exist).
#
# Env: HBCD_RELEASE  — release tag (default: br_21p3)
#       RELEASE_BASE — override release base dir (default: derived)

set -euo pipefail

module load bcftools

DATA_DIR=/projects/standard/basu_hbcd/shared/data
: "${HBCD_RELEASE:=br_21p3}"
: "${RELEASE_BASE:=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_${HBCD_RELEASE}}"

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr="${CHROMS[$SLURM_ARRAY_TASK_ID]}"

KEEP_LIST="${RELEASE_BASE}/keep_list.txt"
IMPUTED_DIRS=("${DATA_DIR}/imputed/c1" "${DATA_DIR}/imputed/c2" "${DATA_DIR}/imputed/c3" "${DATA_DIR}/imputed/cX")
OUTDIR="${RELEASE_BASE}/data/imputed"

mkdir -p "$OUTDIR" logs

# Find the VCF for this chromosome
VCF_FILE=$(find "${IMPUTED_DIRS[@]}" -name "chr${chr}.dose.vcf.gz" ! -name "*.tbi" 2>/dev/null | head -1)

if [ -z "$VCF_FILE" ]; then
  echo "Skipping chr${chr}: VCF not found"
  exit 0
fi

echo "--- Filtering chr${chr} ---"
echo "  Input : $VCF_FILE"
echo "  Output: ${OUTDIR}/chr${chr}.dose.vcf.gz"

# Extract IID column from keep_list (column 2) and use as sample whitelist
# bcftools view --samples-file expects one sample ID per line
awk '{print $2}' "$KEEP_LIST" > "${OUTDIR}/keep_samples_${chr}.txt"

bcftools view \
    "$VCF_FILE" \
    --samples-file "${OUTDIR}/keep_samples_${chr}.txt" \
    --force-samples \
    -Oz -o "${OUTDIR}/chr${chr}.dose.vcf.gz"

tabix -f -p vcf "${OUTDIR}/chr${chr}.dose.vcf.gz"

rm -f "${OUTDIR}/keep_samples_${chr}.txt"

n_before=$(bcftools query -l "$VCF_FILE" | wc -l)
n_after=$(bcftools query -l "${OUTDIR}/chr${chr}.dose.vcf.gz" | wc -l)
echo "  Samples: ${n_before} → ${n_after}"

# Also copy the .info.gz if it exists (no sample filtering needed, but keep paired)
INFO_FILE=$(find "${IMPUTED_DIRS[@]}" -name "chr${chr}.info.gz" 2>/dev/null | head -1)
if [ -n "$INFO_FILE" ]; then
  cp "$INFO_FILE" "${OUTDIR}/chr${chr}.info.gz"
  echo "  Copied .info.gz"
fi

echo "--- Finished chr${chr} ---"
