#!/bin/bash
set -euo pipefail

# step 15: submit VCF chunks to TOPMed Imputation Server (autosomes)

OUTDIR="/home/ood-coffm049/hbcdData/data"
BASE_URL="https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2"

TIS_TOKEN=$(cat ~/topmedKey)
if [ -z "$TIS_TOKEN" ]; then
    echo "Error: ~/topmedKey is empty or missing"
    exit 1
fi

CHUNKS=(
    "1,2,3,4,5,6,7"
    "8,9,10,11,12,13,14"
    "15,16,17,18,19,20,21,22,X"
)

for i in "${!CHUNKS[@]}"; do
    chunk="${CHUNKS[$i]}"
    name="hbcd_chunk$((i+1))"

    IFS=',' read -ra chrs <<< "$chunk"
    file_args=()
    for chr in "${chrs[@]}"; do
        file_args+=(-F "files=@${OUTDIR}/chr${chr}.vcf.gz")
    done

    echo "=== Submitting chunk $((i+1)): chr${chunk} ==="

    curl -H "X-Auth-Token: ${TIS_TOKEN}" \
         -F "refpanel=topmed-r3" \
         -F "build=hg38" \
         -F "phasing=eagle" \
         -F "aesEncryption=true" \
         -F "population=all" \
         -F "job-name=${name}" \
         "${file_args[@]}" \
         "${BASE_URL}/jobs/submit/imputationserver2"
    echo ""
done
