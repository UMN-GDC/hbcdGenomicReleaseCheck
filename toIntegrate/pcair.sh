#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=18:00:00
#SBATCH -p agsmall
#SBATCH -o logs/pcair_%j.out
#SBATCH -e logs/pcair_%j.err
#SBATCH --job-name pcair

module load R/4.4.2-openblas-rocky8

WORK=/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026
PCA_DIR="${WORK}/PCA"
NAME=hbcd_rsid_harmonized

mkdir -p "${PCA_DIR}/pca_ir" "${PCA_DIR}/gds" logs

echo "[$(date)] Starting PC-AiR pipeline"

Rscript "${PCA_DIR}/pca_ir_pipeline.R"

echo "[$(date)] Done."
