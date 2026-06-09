#!/bin/bash
#SBATCH --job-name=downloadImputed3
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=logs/downloadImputed3.out
#SBATCH --error=logs/downloadImputed3.err

# step 19: download imputed results chunk 3 from TOPMed

mkdir -p imputed/c3
cd imputed/c3
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/ZVrtyHaLJ6vcsiLQX8uG7Ry5skC7CjYaxRcpQG4Q | bash
