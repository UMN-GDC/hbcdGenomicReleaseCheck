#!/usr/bin/env bash
# PLINK2 filtering step – called after filter_and_prepare_data.py
# Translates the bash chunk in 01-setupNfilter.qmd

module load plink/2.00-alpha-091019

dataPrefix=/scratch.global/hbcd/full/full.QC8
release=br_20p2
releaseDir=../HBCD_genomics_release_${release}/data/

plink2 \
    --bfile "$dataPrefix" \
    --fam "${releaseDir}../temp.fam" \
    --make-bed --out "${releaseDir}hbcd" \
    --remove "${releaseDir}../Removed_individuals.txt"

cd "${releaseDir}" || exit
rm -f *.log
