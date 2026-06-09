#!/bin/bash
#SBATCH --job-name=downloadImputed2
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/downloadImputed2.out
#SBATCH --error=logs/downloadImputed2.err

# step 18: download imputed results chunk 2 from TOPMed

mkdir -p imputed/c2
cd imputed/c2
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/sLOyMzREN3tdbTgBSkSHLO9p4DZOmefqZobT036A | bash
