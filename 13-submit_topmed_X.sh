#!/bin/bash
set -euo pipefail

# step 16: submit chrX VCF to TOPMed Imputation Server

OUTDIR="/scratch.global/GDC/hbcdGenomicsPreImputation/data"
BASE_URL="https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2"

TIS_TOKEN=$(cat ~/topmedKey)
if [ -z "$TIS_TOKEN" ]; then
    echo "Error: ~/topmedKey is empty or missing"
    exit 1
fi

echo "=== Submitting X Chr ==="

curl -H "X-Auth-Token: ${TIS_TOKEN}" \
     -F "refpanel=topmed-r3" \
     -F "build=hg38" \
     -F "phasing=eagle" \
     -F "aesEncryption=true" \
     -F "population=all" \
     -F "job-name=hbcdX" \
     -F "files=@${OUTDIR}/chrX.vcf.gz" \
     "${BASE_URL}/jobs/submit/imputationserver2"
