#!/usr/bin/env bash
# step 4: PLINK2 --keep to filter onlyQc → hbcd (the release set)
# Usage: RELEASE_DIR=/path/to/release/data ./26-run_plink_filter.sh
# Defaults to HBCD_genomics_release_br_21p2/data when HBCD_RELEASE / RELEASE_DIR unset.

module load plink/2.00-alpha-091019

dataDIR=/projects/standard/basu_hbcd/shared/data
: "${RELEASE_DIR:=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_br_21p2/data}"
releaseDir="${RELEASE_DIR}"

plink2 \
    --bfile "${dataDIR}/onlyQc" \
    --allow-extra-chr \
    --fam "${releaseDir}/../temp.fam" \
    --keep "${releaseDir}/../keep_list.txt" \
    --make-bed --out "${releaseDir}/hbcd"

# Ensure batch.info contains exactly the same subjects as hbcd.fam
awk 'NR==FNR {keep[$2]; next} FNR==1 || $1 in keep' \
    "${releaseDir}/hbcd.fam" "${releaseDir}/batch.info" \
    > "${releaseDir}/batch.info.tmp" \
    && mv "${releaseDir}/batch.info.tmp" "${releaseDir}/batch.info"

cd "${releaseDir}" || exit
rm -f *.log
