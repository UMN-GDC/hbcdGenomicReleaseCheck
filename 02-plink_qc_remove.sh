#!/usr/bin/env bash
# step 2: apply QC removal list to produce onlyQc.{bed,bim,fam}
# Uses the remapped .fam from step 1 so IIDs are release_candid-based.

module load plink/2.00-alpha-091019

HST_DIR=/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_July2025
DATA_DIR=/projects/standard/basu_hbcd/shared/data
REMOVE_FILE="${DATA_DIR}/Remove.txt"
OUT_PREFIX="${DATA_DIR}/onlyQc"

# If Remove.txt is empty or missing, just copy the remapped .fam
if [ ! -s "$REMOVE_FILE" ]; then
    echo "  Remove.txt not found or empty — copying HBCD_remapped.fam as onlyQc.fam"
    cp "${DATA_DIR}/HBCD_remapped.fam" "${OUT_PREFIX}.fam"
    exit 0
fi

plink2 \
    --bfile "${HST_DIR}/HBCD/HBCD" \
    --allow-extra-chr \
    --fam "${DATA_DIR}/HBCD_remapped.fam" \
    --remove "$REMOVE_FILE" \
    --make-bed \
    --out "$OUT_PREFIX"
