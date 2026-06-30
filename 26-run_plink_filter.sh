#!/usr/bin/env bash
# step 4: PLINK2 --keep → GDA/merged_chroms (release PLINK set)
# Usage: HBCD_RELEASE=br_21p3 ./26-run_plink_filter.sh
#        RELEASE_DIR=/path/to/genotype_microarray ./26-run_plink_filter.sh
# Default derived from HBCD_RELEASE (default br_21p3).

module load plink/2.00-alpha-091019

dataDIR=/projects/standard/basu_hbcd/shared/data
: "${HBCD_RELEASE:=br_21p3}"
: "${RELEASE_DIR:=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_${HBCD_RELEASE}/genotype_microarray}"
releaseRoot="${RELEASE_DIR}"
gdaDir="${releaseRoot}/GDA"
releaseBase="$(dirname "$releaseRoot")"

mkdir -p "$gdaDir"

plink2 \
    --bfile "${dataDIR}/onlyQc" \
    --allow-extra-chr \
    --fam "${releaseBase}/temp.fam" \
    --keep "${releaseBase}/keep_list.txt" \
    --make-bed --out "${gdaDir}/merged_chroms"

# Ensure batch.info contains exactly the same subjects as merged_chroms.fam
awk 'NR==FNR {keep[$2]; next} FNR==1 || $1 in keep' \
    "${gdaDir}/merged_chroms.fam" "${gdaDir}/batch.info" \
    > "${gdaDir}/batch.info.tmp" \
    && mv "${gdaDir}/batch.info.tmp" "${gdaDir}/batch.info"

cd "${gdaDir}" || exit
rm -f *.log
