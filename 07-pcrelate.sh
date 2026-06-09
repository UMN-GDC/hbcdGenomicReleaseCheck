#!/bin/bash -l
#SBATCH --job-name=pcrelate
#SBATCH --output=logs/pcrelate_%j.out
#SBATCH --error=logs/pcrelate_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH -p agsmall

# step 10: SLURM wrapper for PC-Relate pipeline

module load R/4.4.2-openblas-rocky8

DATA_DIR=/projects/standard/basu_hbcd/shared/data

mkdir -p logs

echo "[$(date)] Starting PC-Relate pipeline"
Rscript "08-pca_relate_pipeline.R"
echo "[$(date)] Done."
