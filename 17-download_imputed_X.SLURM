#!/bin/bash
#SBATCH --job-name=downloadImputed3
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/downloadImputed3.out
#SBATCH --error=logs/downloadImputed3.err

# step 20: download imputed results chrX from TOPMed

mkdir -p imputed/cX
cd imputed/cX
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/5J32bQCgyPOhRNJmcrtesSqpTG83tYFeiHlsWH3m | bash
