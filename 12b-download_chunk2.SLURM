#!/bin/bash
#SBATCH --job-name=downloadImputed2
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/downloadImputed2.out
#SBATCH --error=logs/downloadImputed2.err

# Get the download link hash from the TopMed Server job results summary
OUT=~/hbcdData/imputed

mkdir -p ${OUT}/c2
cd ${OUT}/c2
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/sxoj5mxhxxcI0ktpDebF4AvdA0EIH1gk6v22SyFs | bash
