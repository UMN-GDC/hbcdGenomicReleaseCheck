#!/bin/bash
#SBATCH --job-name=downloadImputed1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/downloadImputed1.out
#SBATCH --error=logs/downloadImputed1.err

# step 17: download imputed results chunk 1 from TOPMed
# Get the download link hash from the TopMed Server job results summary
OUT=~/hbcdData/imputed

mkdir -p ${OUT}/c1
cd ${OUT}/c1
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/3U8kWNE4ckhHHlP2elXdHaJ2a1Uv9bLksiPzfiRz | bash
