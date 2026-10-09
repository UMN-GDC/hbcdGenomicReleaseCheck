#!/bin/bash
# copy_imputation_qc.sh — Copy imputation QC statistics from source to release
# Can be run standalone or sourced by 27-filter_imputed_vcf.SLURM
#
# Usage:
#   ./copy_imputation_qc.sh [--release-tag TAG]
#   or
#   source copy_imputation_qc.sh && copy_imputation_qc

set -euo pipefail

copy_imputation_qc() {
    # Allow override via env vars
    local release="${HBCD_RELEASE:-br_21p3}"
    local data_dir="${HBCD_DATA_DIR:-/projects/standard/basu_hbcd/shared/data}"
    local imp_dir="${HBCD_IMPUTATION_DIR:-/projects/standard/basu_hbcd/shared/hbcdSandboxData}"

    local impdir="${imp_dir}/imputed"
    local release_base="$(dirname "$data_dir")/HBCD_genomics_release_${release}"
    local outdir="${release_base}/genotype_microarray/imputed"

    mkdir -p "$outdir"

    echo "Copying imputation QC statistics to $outdir"
    echo "  Release: $release"
    echo "  Source:  $impdir"
    echo "  Target:  $outdir"

    local copied=0 missing=0

    for i in 1 2 3; do
        local chunk_dir="${impdir}/c${i}"
        local stats_dir="${chunk_dir}/statistics"

        if [ ! -d "$stats_dir" ]; then
            echo "  WARNING: $stats_dir not found, skipping chunk c${i}"
            continue
        fi

        for f in snps-typed-only.txt snps-excluded.txt chunks-excluded.txt; do
            local src="${stats_dir}/${f}"
            if [ -f "$src" ]; then
                local dst="${outdir}/batch${i}-${f}"
                cp "$src" "$dst"
                echo "  Copied: $src → $dst"
                ((copied++))
            else
                echo "  MISSING: $src"
                ((missing++))
            fi
        done

        # quality-control.html — check stats_dir first, then chunk root
        local src="${stats_dir}/quality-control.html"
        if [ ! -f "$src" ]; then
            src="${chunk_dir}/quality-control.html"
        fi
        if [ -f "$src" ]; then
            local dst="${outdir}/batch${i}-quality-control.html"
            cp "$src" "$dst"
            echo "  Copied: $src → $dst"
            ((copied++))
        else
            echo "  MISSING: quality-control.html for chunk c${i}"
            ((missing++))
        fi
    done

    echo "Done. Copied: $copied, Missing: $missing"
    echo "Files in $outdir:"
    ls -la "$outdir"/batch*-*.txt "$outdir"/batch*-quality-control.html 2>/dev/null || true
}

# Allow running standalone
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    copy_imputation_qc
fi