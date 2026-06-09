#!/bin/bash
#SBATCH --job-name=unzip_X
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/unzip_X.out
#SBATCH --error=logs/unzip_X.err

# step 24: unzip imputed chrX

cd imputed/cX

/home/ood-coffm049/software/7z/7zz x chr_X.zip -p"$(cat password.txt)"
