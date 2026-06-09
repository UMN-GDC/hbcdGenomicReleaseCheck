#!/bin/bash
#SBATCH --job-name=downloadImputed1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/downloadImputed1.out
#SBATCH --error=logs/downloadImputed1.err

# step 17: download imputed results chunk 1 from TOPMed
# Get the download link hash from the TopMed Server job results summary

mkdir -p imputed/c1
cd imputed/c1
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/woDH5ojmb2wA48E1C9TgnNM8nz5SUlUTzHt6Ccbu | bash
