#!/usr/bin/env bash
# PLINK2 filtering step – called after filter_and_prepare_data.py
# Translates the bash chunk in 01-setupNfilter.qmd

module load plink/2.00-alpha-091019

dataDIR=/projects/standard/basu_hbcd/shared/data
dataPrefix=/projects/standard/basu_hbcd/shared/archive/HBCD_genomics_release_br_20p1/data/hbcd
release=br_21p2
releaseDir=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_${release}/data

awk '{print $2}' ${dataPrefix}.bim > ${releaseDir}/../Extracted_variants.txt

plink2 \
    --bfile $dataPrefix \
    --allow-extra-chr \
    --set-all-var-ids chr@_#_\$r_\$a_b38 \
    --extract ${dataDIR}/Extracted_variants.txt\
    --keep "${releaseDir}/../keep_list.txt" \
    --make-bed --out "${releaseDir}/hbcd"

# Ensure batch.info contains exactly the same subjects as hbcd.fam
awk 'NR==FNR {keep[$2]; next} $1 in keep' \
    "${releaseDir}/hbcd.fam" "${releaseDir}/batch.info" \
    > "${releaseDir}/batch.info.tmp" \
    && mv "${releaseDir}/batch.info.tmp" "${releaseDir}/batch.info"

cd "${releaseDir}" || exit
rm -f *.log
