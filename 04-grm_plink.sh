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

# step 7: compute PLINK2 GRM (--make-rel) from harmonized PLINK files
# Output: .rel (lower-triangular), .rel.id (sample order)

set -euo pipefail

module load plink

PLINK_STEM=/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026/AE/rsid_updated/hbcd_rsid_harmonized
DATA_DIR=/projects/standard/basu_hbcd/shared/data
OUT_STEM="${DATA_DIR}/hbcd_plink_grm"

mkdir -p logs

echo "[$(date)] Starting PLINK GRM computation"
echo "  Input : ${PLINK_STEM}"
echo "  Output: ${OUT_STEM}"

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
