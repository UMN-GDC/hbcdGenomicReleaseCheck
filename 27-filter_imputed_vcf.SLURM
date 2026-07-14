#!/bin/bash
#SBATCH --job-name=filter_vcf
#SBATCH --array=0-23
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/filter_vcf_%A_%a.out
#SBATCH --error=logs/filter_vcf_%A_%a.err

# step 27: filter imputed VCFs to release subjects and de-identify sample IDs
# Sample IDs are renamed from pscid-based to release_candid-based (10-digit + C/M).
# Run AFTER steps 25-26 (keep_list.txt exists).
#
# Env: HBCD_RELEASE  — release tag (default: br_21p3)
#       RELEASE_BASE — override release base dir (default: derived)

set -euo pipefail

module load bcftools

source ~/miniconda3/etc/profile.d/conda.sh
conda activate python

DATA_DIR=/home/ood-coffm049/hbcdData
IMPDIR="${DATA_DIR}/imputed"
: "${HBCD_RELEASE:=br_21p3}"
: "${RELEASE_BASE:=${DATA_DIR}/HBCD_genomics_release_${HBCD_RELEASE}}"

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
# Step 1: generate rename mapping and keep list using identifiers
# ------------------------------------------------------------------

cat > "${OUTDIR}/.build_mapping_${chr}.py" << 'PYEOF'
import sys
import pandas as pd
from pathlib import Path

data_dir = Path(sys.argv[1])
keep_list = Path(sys.argv[2])
outdir = Path(sys.argv[3])
chr_tag = sys.argv[4]

# Load identifiers crosswalk
ids = pd.read_csv(data_dir / "release_identifiers_20260628.csv")
ids = ids[ids["release_candid"] != "release_candid"]
ids["release_candid"] = pd.to_numeric(ids["release_candid"])
ids = ids.dropna(subset=["release_candid"])
ids["pscid"] = ids["pscid"].astype(str).str.strip()
pscid_to_rc = dict(zip(ids["pscid"], ids["release_candid"].astype(int).astype(str)))

# Load release IIDs
keep = pd.read_csv(keep_list, sep=r"\s+", header=None, names=["FID", "IID"])
release_iids = set(keep["IID"].astype(str))

# Read VCF sample names
with open(outdir / f"vcf_samples_{chr_tag}.txt") as f:
    vcf_samples = [line.strip() for line in f if line.strip()]

rename_lines = []
keep_lines = []
n_skip_no_rel = 0
n_skip_no_rc = 0
n_skip_not_release = 0

for s in vcf_samples:
    if len(s) < 2:
        n_skip_no_rel += 1
        continue
    pscid = s[:-1]
    rel = s[-1]
    if rel not in ("C", "M"):
        n_skip_no_rel += 1
        continue
    rc = pscid_to_rc.get(pscid)
    if rc is None:
        n_skip_no_rc += 1
        continue
    new_iid = f"{rc}{rel}"
    if new_iid not in release_iids:
        n_skip_not_release += 1
        continue
    rename_lines.append(f"{s}\t{new_iid}")
    keep_lines.append(s)

(outdir / f"rename_{chr_tag}.txt").write_text("\n".join(rename_lines) + "\n")
(outdir / f"keep_old_{chr_tag}.txt").write_text("\n".join(keep_lines) + "\n")

print(f"  VCF samples             : {len(vcf_samples)}")
print(f"  Release subjects kept   : {len(keep_lines)}")
print(f"  Skip — no C/M suffix    : {n_skip_no_rel}")
print(f"  Skip — unmapped pscid   : {n_skip_no_rc}")
print(f"  Skip — not in release   : {n_skip_not_release}")
PYEOF

# Extract VCF sample names first
bcftools query -l "$VCF_FILE" > "${OUTDIR}/vcf_samples_${chr}.txt"

python3 "${OUTDIR}/.build_mapping_${chr}.py" \
    "$DATA_DIR" "$KEEP_LIST" "$OUTDIR" "$chr"

n_keep=$(wc -l < "${OUTDIR}/keep_old_${chr}.txt" | tr -d ' ')

if [ "$n_keep" -eq 0 ]; then
  echo "  WARNING: No release subjects found in VCF for chr${chr} — writing empty VCF"
  # Write an empty (header-only) VCF so tabix doesn't fail downstream
  bcftools view -h "$VCF_FILE" | bgzip -c > "${OUTDIR}/chr${chr}_dose.vcf.gz"
  tabix -f -p vcf "${OUTDIR}/chr${chr}_dose.vcf.gz"
  rm -f "${OUTDIR}/.build_mapping_${chr}.py" \
        "${OUTDIR}/vcf_samples_${chr}.txt" \
        "${OUTDIR}/keep_old_${chr}.txt" \
        "${OUTDIR}/rename_${chr}.txt"
  echo "--- Finished chr${chr} (empty) ---"
  exit 0
fi

# ------------------------------------------------------------------
# Step 2: filter VCF to release subjects (by original sample name)
# ------------------------------------------------------------------

bcftools view \
    "$VCF_FILE" \
    --samples-file "${OUTDIR}/keep_old_${chr}.txt" \
    --force-samples \
    -Oz -o "${OUTDIR}/chr${chr}_dose_tmp.vcf.gz"

# ------------------------------------------------------------------
# Step 3: rename samples to de-identified IIDs
# ------------------------------------------------------------------

bcftools reheader \
    --samples "${OUTDIR}/rename_${chr}.txt" \
    -o "${OUTDIR}/chr${chr}_dose.vcf.gz" \
    "${OUTDIR}/chr${chr}_dose_tmp.vcf.gz"

# ------------------------------------------------------------------
# Step 4: index and cleanup
# ------------------------------------------------------------------

tabix -f -p vcf "${OUTDIR}/chr${chr}_dose.vcf.gz"

rm -f "${OUTDIR}/.build_mapping_${chr}.py" \
      "${OUTDIR}/vcf_samples_${chr}.txt" \
      "${OUTDIR}/keep_old_${chr}.txt" \
      "${OUTDIR}/rename_${chr}.txt" \
      "${OUTDIR}/chr${chr}_dose_tmp.vcf.gz"

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
