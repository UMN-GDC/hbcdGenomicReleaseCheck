#!/bin/bash
#SBATCH --job-name=unzip_c3_%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/unzip_X.out
#SBATCH --error=logs/unzip_X.err

OUT=~/hbcdData/imputed
cd $OUT/c3

# Grab the password from the TopMed confirmation email
/home/ood-coffm049/software/7z/7zz x chr_X.zip -p"$(cat password.txt)"
