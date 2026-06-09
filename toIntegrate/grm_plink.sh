#!/bin/bash -l
#SBATCH --job-name=grm_plink
#SBATCH --output=logs/grm_plink_%j.out
#SBATCH --error=logs/grm_plink_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH -p agsmall

# =============================================================================
# grm_plink.sh
#
# Computes a Genetic Relatedness Matrix (GRM) using PLINK2 --make-rel
# from hbcd_rsid_harmonized PLINK files. Output written to PCA directory.
#
# PLINK2 produces three files:
#   .rel      — lower-triangular relatedness matrix (text)
#   .rel.id   — sample IDs corresponding to matrix rows/columns
# =============================================================================

set -euo pipefail

module load plink

PLINK_STEM=/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026/AE/rsid_updated/hbcd_rsid_harmonized
OUT_DIR=/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026/PCA
OUT_STEM="${OUT_DIR}/hbcd_plink_grm"

mkdir -p "${OUT_DIR}" logs

echo "[$(date)] Starting PLINK GRM computation"
echo "  Input : ${PLINK_STEM}"
echo "  Output: ${OUT_STEM}"

# ---------------------------------------------------------------------------
# PLINK --make-rel computes the genomic relationship matrix
# 'cov' = use covariance form (standard GRM); default is correlation form
# 'bin' = write binary output (much faster to read downstream)
# ---------------------------------------------------------------------------
plink \
    --bfile   "${PLINK_STEM}" \
    --make-rel square \
    --out     "${OUT_STEM}" \
    --memory  60000 \
    --threads "${SLURM_CPUS_PER_TASK:-10}"

echo ""
echo "===== PLINK GRM Summary ====="
echo "  Output matrix : ${OUT_STEM}.rel"
echo "  Sample IDs    : ${OUT_STEM}.rel.id"
echo "  N samples     : $(wc -l < ${OUT_STEM}.rel.id)"
echo "============================="
echo "[$(date)] Done."
