#!/bin/bash
#SBATCH --job-name=unzip_c3_%A_%a.out
#SBATCH --array=16-22
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/unzip_c3_%A_%a.out
#SBATCH --error=logs/unzip_c3_%A_%a.err

# step 23: unzip imputed chunk 3 (chrs 16-22)

cd imputed/c3

/home/ood-coffm049/software/7z/7zz x chr_${SLURM_ARRAY_TASK_ID}.zip -p"$(cat password.txt)"
