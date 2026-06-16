#!/bin/bash
# Auto-detect ancestries from ses-V02 rows and submit 14a-precompute_gp.SLURM.
#
# Usage:
#   bash 14-submit_gp_prep.sh              # submit all ancestries
#   bash 14-submit_gp_prep.sh --setup-only # preview without submitting

set -euo pipefail

ANC_SRC="/shared/release/hbcd/hbcd/rawdata/phenotype/sed_basic_demographics.tsv"
N_CHR=23
SCRIPT_DIR=$(dirname "$0")

# Map numeric race codes to abbreviations (must match 14a-precompute_gp.SLURM)
ANC_GROUPS=$(awk -F'\t' '
  BEGIN{split("WHT BLK AIAN ASN HPI 2PLUS OTH UNK",a); for(i=0;i<8;i++) m[i]=a[i+1]}
  NR>1 && $2=="ses-V02" && $4!="" {print m[$4]}
' "$ANC_SRC" | sort -u | paste -sd, -)
N_ANC=$(echo "$ANC_GROUPS" | awk -F',' '{print NF}')

echo "Detected $N_ANC ancestry groups: $ANC_GROUPS (ses-V02 only)"

if [ "${1:-}" = "--setup-only" ]; then
    echo "sbatch --array=0-$((N_CHR * N_ANC - 1)) --export=ALL,ANCESTRY_GROUPS=\"$ANC_GROUPS\" $SCRIPT_DIR/14a-precompute_gp.SLURM"
    exit 0
fi

sbatch --array="0-$((N_CHR * N_ANC - 1))" \
       --export="ALL,ANCESTRY_GROUPS=$ANC_GROUPS" \
       "$SCRIPT_DIR/14a-precompute_gp.SLURM"
