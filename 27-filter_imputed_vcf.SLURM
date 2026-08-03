#!/bin/bash
#SBATCH --job-name=filter_vcf
#SBATCH --array=0-23
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/filter_vcf_%A_%a.out
#SBATCH --error=logs/filter_vcf_%A_%a.err

# step 27: filter imputed VCFs to release subjects (samples already de-identified)
# VCF sample IDs are already release_candid-based (10-digit + C/M); just filter to release.
# Run AFTER steps 25-26 (keep_list.txt exists).
#
# Env: HBCD_RELEASE         — release tag (default: br_21p3)
#       HBCD_DATA_DIR        — source data dir (default: shared/data)
#       HBCD_IMPUTATION_DIR  — imputation data root (default: shared/hbcdSandboxData)
#       RELEASE_BASE         — override release base dir (default: derived from DATA_DIR)

set -euo pipefail

module load bcftools

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline

DATA_DIR="${HBCD_DATA_DIR:-/projects/standard/basu_hbcd/shared/data}"
: "${HBCD_RELEASE:=br_21p3}"
IMPUTATION_DIR="${HBCD_IMPUTATION_DIR:-/projects/standard/basu_hbcd/shared/hbcdSandboxData}"
IMPDIR="${IMPUTATION_DIR}/imputed"
: "${RELEASE_BASE:=$(dirname "$DATA_DIR")/HBCD_genomics_release_${HBCD_RELEASE}}"

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr="${CHROMS[$SLURM_ARRAY_TASK_ID]}"

KEEP_LIST="${RELEASE_BASE}/keep_list.txt"
IMPUTED_DIRS=("${IMPDIR}/c1" "${IMPDIR}/c2" "${IMPDIR}/c3")
OUTDIR="${RELEASE_BASE}/genotype_microarray/imputed"

mkdir -p "$OUTDIR" logs

# Find the VCF for this chromosome
VCF_FILE=$(find "${IMPUTED_DIRS[@]}" -name "chr${chr}.dose.vcf.gz" ! -name "*.tbi" 2>/dev/null | head -1)

if [ -z "$VCF_FILE" ]; then
  echo "Skipping chr${chr}: VCF not found"
  exit 0
fi

echo "--- Filtering & de-identifying chr${chr} ---"
echo "  Input : $VCF_FILE"
echo "  Output: ${OUTDIR}/chr${chr}_dose.vcf.gz"

# ------------------------------------------------------------------
# Step 1: find keep-list / VCF sample intersection
# ------------------------------------------------------------------

cat > "${OUTDIR}/.build_keep_${chr}.py" << 'PYEOF'
import sys
from pathlib import Path

keep_list = Path(sys.argv[1])
outdir = Path(sys.argv[2])
chr_tag = sys.argv[3]

# Load release IIDs from keep_list.txt (FID IID)
keep = set()
with open(keep_list) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            keep.add(parts[1])

# Read VCF sample names
with open(outdir / f"vcf_samples_{chr_tag}.txt") as f:
    vcf_samples = [line.strip() for line in f if line.strip()]

keep_lines = [s for s in vcf_samples if s in keep]

(outdir / f"keep_old_{chr_tag}.txt").write_text("\n".join(keep_lines) + "\n")

print(f"  VCF samples             : {len(vcf_samples)}")
print(f"  Release subjects kept   : {len(keep_lines)}")
print(f"  Excluded                : {len(vcf_samples) - len(keep_lines)}")
PYEOF

# Extract VCF sample names first
bcftools query -l "$VCF_FILE" > "${OUTDIR}/vcf_samples_${chr}.txt"

python3 "${OUTDIR}/.build_keep_${chr}.py" \
    "$KEEP_LIST" "$OUTDIR" "$chr"

n_keep=$(wc -l < "${OUTDIR}/keep_old_${chr}.txt" | tr -d ' ')

if [ "$n_keep" -eq 0 ]; then
  echo "  WARNING: No release subjects found in VCF for chr${chr} — writing empty VCF"
  bcftools view -h "$VCF_FILE" | bgzip -c > "${OUTDIR}/chr${chr}_dose.vcf.gz"
  tabix -f -p vcf "${OUTDIR}/chr${chr}_dose.vcf.gz"
  rm -f "${OUTDIR}/.build_keep_${chr}.py" \
        "${OUTDIR}/vcf_samples_${chr}.txt" \
        "${OUTDIR}/keep_old_${chr}.txt"
  echo "--- Finished chr${chr} (empty) ---"
  exit 0
fi

# ------------------------------------------------------------------
# Step 2: filter VCF to release subjects
# ------------------------------------------------------------------

bcftools view \
    "$VCF_FILE" \
    --samples-file "${OUTDIR}/keep_old_${chr}.txt" \
    --force-samples \
    -Oz -o "${OUTDIR}/chr${chr}_dose.vcf.gz"

# ------------------------------------------------------------------
# Step 3: index and cleanup
# ------------------------------------------------------------------

tabix -f -p vcf "${OUTDIR}/chr${chr}_dose.vcf.gz"

rm -f "${OUTDIR}/.build_keep_${chr}.py" \
      "${OUTDIR}/vcf_samples_${chr}.txt" \
      "${OUTDIR}/keep_old_${chr}.txt"

n_before=$(bcftools query -l "$VCF_FILE" | wc -l)
n_after=$(bcftools query -l "${OUTDIR}/chr${chr}_dose.vcf.gz" | wc -l)
echo "  Samples: ${n_before} → ${n_after}"

# Also copy the .info.gz if it exists (no sample filtering needed, but keep paired)
INFO_FILE=$(find "${IMPUTED_DIRS[@]}" -name "chr${chr}.info.gz" 2>/dev/null | head -1)
if [ -n "$INFO_FILE" ]; then
  cp "$INFO_FILE" "${OUTDIR}/chr${chr}.info.gz"
  echo "  Copied .info.gz"
fi

echo "--- Finished chr${chr} ---"
