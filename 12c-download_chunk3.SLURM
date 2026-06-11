#!/bin/bash
#SBATCH --job-name=downloadImputed3
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/downloadImputed3.out
#SBATCH --error=logs/downloadImputed3.err

# Get the download link hash from the TopMed Server job results summary
OUT=~/hbcdData/imputed

mkdir -p ${OUT}/c3
cd ${OUT}/c3
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/UApZdP8w5Eq335Pm1namkVXRBWOzAVsWTHKh2vCW | bash
