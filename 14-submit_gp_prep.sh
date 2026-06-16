#!/bin/bash
# Submit 14a-precompute_gp.SLURM with per-ancestry parallelism.
# All 8 race groups are submitted; groups with no data produce empty results.
#
# Usage:
#   bash 14-submit_gp_prep.sh              # submit
#   bash 14-submit_gp_prep.sh --setup-only # preview

set -euo pipefail

N_CHR=23
SCRIPT_DIR=$(dirname "$0")
ANC_GROUPS="WHT,BLK,AIAN,ASN,HPI,2PLUS,OTH,UNK"
N_ANC=8

N_TASKS=$((N_CHR * N_ANC - 1))
echo "Ancestry groups: $ANC_GROUPS"
echo "Submitting array 0-$N_TASKS ($((N_TASKS+1)) tasks)"

if [ "${1:-}" = "--setup-only" ]; then
    echo "sbatch --array=0-${N_TASKS} --export=ALL,ANCESTRY_GROUPS=\"$ANC_GROUPS\" $SCRIPT_DIR/14a-precompute_gp.SLURM"
    exit 0
fi

sbatch --array="0-${N_TASKS}" \
       --export="ALL,ANCESTRY_GROUPS=$ANC_GROUPS" \
       "$SCRIPT_DIR/14a-precompute_gp.SLURM"
