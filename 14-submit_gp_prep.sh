#!/bin/bash
# Auto-detect ancestries and submit 14a-precompute_gp.SLURM with the right array size.
#
# Usage:
#   bash submit_gp_prep.sh              # all ancestries, auto-detect
#   bash submit_gp_prep.sh --setup-only # just build ancestry cache, don't submit

set -euo pipefail

IMPDIR=/home/ood-coffm049/hbcdData/imputed
CACHE_DIR=$IMPDIR/gp_cache
ANC_CACHE="$CACHE_DIR/sample_ancestry.txt"
ANC_SRC="/shared/release/hbcd/hbcd/rawdata/phenotype/sed_basic_demographics.tsv"
N_CHR=23
SCRIPT_DIR=$(dirname "$0")

mkdir -p "$CACHE_DIR"

# ── Detect ancestries ──────────────────────────────────────────────────
if [ -f "$ANC_CACHE" ]; then
    echo "Reading ancestries from cache: $ANC_CACHE"
    ANC_GROUPS=$(awk '!seen[$2]++{printf "%s,",$2}' "$ANC_CACHE" | sed 's/,$//')
    N_ANC=$(awk '!seen[$2]++{c++} END{print c}' "$ANC_CACHE")
else
    echo "Detecting ancestries from: $ANC_SRC"
    if [ -f "$ANC_SRC" ]; then
        ANC_GROUPS=$(awk -F'\t' '
            NR==1 {for(i=1;i<=NF;i++) if($i=="sed_basic_demographics_child_race") col=i}
            NR>1 && $col!="" {print $col}
        ' "$ANC_SRC" | sort -u | paste -sd, -)
    fi
    N_ANC=$(echo "$ANC_GROUPS" | awk -F',' '{print NF}')
    echo "  Found $N_ANC groups: $ANC_GROUPS"
fi

# ── Handle --setup-only ────────────────────────────────────────────────
if [ "${1:-}" = "--setup-only" ]; then
    echo "Setup mode. To submit the array job, run:"
    if [ -n "$ANC_GROUPS" ]; then
        N_TASKS=$((N_CHR * N_ANC - 1))
        echo "  sbatch --array=0-${N_TASKS} \\"
        echo "         --export=ALL,ANCESTRY_GROUPS=\"$ANC_GROUPS\" \\"
        echo "         $SCRIPT_DIR/14a-precompute_gp.SLURM"
    else
        echo "  sbatch $SCRIPT_DIR/14a-precompute_gp.SLURM"
    fi
    exit 0
fi

# ── Submit ─────────────────────────────────────────────────────────────
if [ -n "$ANC_GROUPS" ]; then
    N_TASKS=$((N_CHR * N_ANC - 1))
    echo "Detected $N_ANC ancestry groups: $ANC_GROUPS"
    echo "Submitting array 0-$N_TASKS ($((N_TASKS+1)) tasks)"

    sbatch --array="0-${N_TASKS}" \
           --export="ALL,ANCESTRY_GROUPS=$ANC_GROUPS" \
           "$SCRIPT_DIR/14a-precompute_gp.SLURM"
else
    echo "Submitting per-chromosome array 0-$((N_CHR-1)) ($N_CHR tasks)"
    sbatch "$SCRIPT_DIR/14a-precompute_gp.SLURM"
fi
