#!/bin/bash
#SBATCH --job-name=unzip_c2_%A_%a.out
#SBATCH --array=8-15
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/unzip_c2_%A_%a.out
#SBATCH --error=logs/unzip_c2_%A_%a.err

# step 22: unzip imputed chunk 2 (chrs 8-15)

cd imputed/c2

/home/ood-coffm049/software/7z/7zz x chr_${SLURM_ARRAY_TASK_ID}.zip -p"$(cat password.txt)"
