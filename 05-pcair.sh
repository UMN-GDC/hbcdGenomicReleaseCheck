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

# step 8: SLURM wrapper for PC-AiR pipeline

module load R/4.4.2-openblas-rocky8

DATA_DIR=/projects/standard/basu_hbcd/shared/data

mkdir -p "${DATA_DIR}/gds" logs

echo "[$(date)] Starting PC-AiR pipeline"
Rscript "06-pca_ir_pipeline.R"
echo "[$(date)] Done."
