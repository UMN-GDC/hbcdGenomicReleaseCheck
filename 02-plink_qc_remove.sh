#!/usr/bin/env bash
# step 2: apply QC removal list + control removal to produce onlyQc.{bed,bim,fam}
# Uses the remapped .fam from step 1 so IIDs are release_candid-based.

set -euo pipefail

module load plink/2.00-alpha-091019

HST_DIR=/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_July2025
DATA_DIR=/projects/standard/basu_hbcd/shared/data
REMOVE_FILE="${DATA_DIR}/Remove.txt"
CONTROLS_FILE="${DATA_DIR}/Removed_controls.txt"
MERGED_REMOVE="${DATA_DIR}/_merged_remove.txt"
OUT_PREFIX="${DATA_DIR}/onlyQc"

# Merge remove lists into one (skip missing/empty)
: > "$MERGED_REMOVE"
for f in "$REMOVE_FILE" "$CONTROLS_FILE"; do
    if [ -s "$f" ]; then
        cat "$f" >> "$MERGED_REMOVE"
    fi
done

if [ ! -s "$MERGED_REMOVE" ]; then
    echo "  No remove files found — copying HBCD_remapped.fam as onlyQc.fam"
    cp "${DATA_DIR}/HBCD_remapped.fam" "${OUT_PREFIX}.fam"
    rm -f "$MERGED_REMOVE"
    exit 0
fi

plink2 \
    --bfile "${HST_DIR}/HBCD/HBCD" \
    --allow-extra-chr \
    --fam "${DATA_DIR}/HBCD_remapped.fam" \
    --remove "$MERGED_REMOVE" \
    --make-bed \
    --out "$OUT_PREFIX"

rm -f "$MERGED_REMOVE"
