#!/bin/bash -l
#SBATCH --job-name=grm_gcta
#SBATCH --output=logs/grm_gcta_%j.out
#SBATCH --error=logs/grm_gcta_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH -p agsmall

# =============================================================================
# grm_gcta.sh
#
# Computes a Genetic Relatedness Matrix (GRM) using GCTA from
# hbcd_rsid_harmonized PLINK files. Output written to PCA directory.
#
# The GRM is computed in chromosome chunks (--make-grm-part) to manage
# memory, then merged into a single GRM.
# =============================================================================

set -euo pipefail

module load intel-oneapi-mkl/2023.1.0-intel-oneapi-mpi-2021.9.0-oneapi-2023.1.0-kanr5ou

GCTA=/projects/standard/gdc/shared/grm_tools/gcta64
PLINK_STEM=/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026/AE/rsid_updated/hbcd_rsid_harmonized
OUT_DIR=/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026/PCA
OUT_STEM="${OUT_DIR}/hbcd_gcta_grm"
N_PARTS=22   # one part per autosome

mkdir -p "${OUT_DIR}" logs

echo "[$(date)] Starting GCTA GRM computation"
echo "  Input : ${PLINK_STEM}"
echo "  Output: ${OUT_STEM}"

# ---------------------------------------------------------------------------
# Step 1 — Compute GRM in parts (one per chromosome chunk)
# This avoids loading all variants into memory at once
# ---------------------------------------------------------------------------
echo "[$(date)] Step 1: Computing GRM in ${N_PARTS} parts"

for PART in $(seq 1 ${N_PARTS}); do
    echo "  Computing part ${PART}/${N_PARTS}..."
    "${GCTA}" \
        --bfile   "${PLINK_STEM}" \
        --make-grm-part ${N_PARTS} ${PART} \
        --out     "${OUT_STEM}" \
        --thread-num "${SLURM_CPUS_PER_TASK:-10}"
done

# ---------------------------------------------------------------------------
# Step 2 — Merge all parts into a single GRM
# ---------------------------------------------------------------------------
echo "[$(date)] Step 2: Merging GRM parts"

# Build list of part files
> "${OUT_DIR}/grm_parts.txt"
for PART in $(seq 1 ${N_PARTS}); do
    echo "${OUT_STEM}.part_${N_PARTS}_${PART}" >> "${OUT_DIR}/grm_parts.txt"
done

"${GCTA}" \
    --mgrm    "${OUT_DIR}/grm_parts.txt" \
    --make-grm \
    --out     "${OUT_STEM}" \
    --thread-num "${SLURM_CPUS_PER_TASK:-10}"

# Clean up part files
rm -f "${OUT_STEM}".part_*.grm.bin \
      "${OUT_STEM}".part_*.grm.N.bin \
      "${OUT_STEM}".part_*.grm.id

echo ""
echo "===== GCTA GRM Summary ====="
echo "  Output GRM  : ${OUT_STEM}.grm.bin"
echo "  Sample IDs  : ${OUT_STEM}.grm.id"
echo "  N samples   : $(wc -l < ${OUT_STEM}.grm.id)"
echo "============================"
echo "[$(date)] Done."
