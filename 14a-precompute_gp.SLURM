#!/bin/bash
#SBATCH --job-name=gp_prep
#SBATCH --array=0-23
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/gp_prep_%A_%a.out
#SBATCH --error=logs/gp_prep_%A_%a.err

# Per-chromosome (optionally × ancestry): extracts GP for rare variants.
#
# Options 2-4 implemented:
#   2 — bcftools view -i INFO/R2 filter (no -R regions file)
#   3 — Python/numpy replaces gawk (10-50x faster)
#   4 — Per-ancestry parallelism via ANCESTRY_GROUPS env var
#
# Usage: run via 14-submit_gp_prep.sh (auto-detects ancestries, submits array).
#
# Not intended to be called directly — use the submit script. Example:
#   bash 14-submit_gp_prep.sh

set -euo pipefail

module load bcftools

source ~/miniconda3/etc/profile.d/conda.sh
conda activate python  # default for Python (main pipeline)

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
IMPDIR=/home/ood-coffm049/hbcdData/imputed
ANCESTRY_FILE="/shared/release/hbcd/hbcd/rawdata/phenotype/sed_basic_demographics.tsv"
# sed_basic_demographics_child_race: 0=White, 1=Black, 2=AIAN, 3=Asian,
#                                     4=Hawaiian/PI, 5=Two+, 6=Other, 7=Unknown
IMPUTED_DIRS=($IMPDIR/c1 $IMPDIR/c2 $IMPDIR/c3)
CACHE_DIR=$IMPDIR/gp_cache
SCRIPT_DIR=$(dirname "$0")
mkdir -p "$CACHE_DIR" logs

# ── Determine chromosome and ancestry mode ────────────────────────────
if [ -n "${ANCESTRY_GROUPS:-}" ]; then
    IFS=',' read -ra ANC_ARY <<< "$ANCESTRY_GROUPS"
    N_ANC=${#ANC_ARY[@]}
    chr_idx=$((SLURM_ARRAY_TASK_ID / N_ANC))
    anc_idx=$((SLURM_ARRAY_TASK_ID % N_ANC))
    chr="${CHROMS[$chr_idx]}"
    ANC="${ANC_ARY[$anc_idx]}"
    MODE_TAG="anc=$ANC"
else
    chr="${CHROMS[$SLURM_ARRAY_TASK_ID]}"
    ANC=""
    MODE_TAG=""
fi

# ── Find VCF and INFO files ──────────────────────────────────────────
VCF_FILE=$(find "${IMPUTED_DIRS[@]}" \
    -name "chr${chr}.dose.vcf.gz" ! -name "*.tbi" 2>/dev/null | head -1)
INFO_FILE=$(find "${IMPUTED_DIRS[@]}" \
    -name "chr${chr}.info.gz" 2>/dev/null | head -1)

if [ -z "$VCF_FILE" ]; then
    echo "Skipping chr${chr}: VCF not found"
    exit 0
fi

echo "--- Processing chr${chr} ${MODE_TAG}---"

# ── Rename VCF samples (remove FID_ prefix) ──────────────────────────
first_sample=$(bcftools query -l "$VCF_FILE" | head -1) || true
if echo "$first_sample" | grep -q '_'; then
    # Lock so ancestry-split jobs don't race on the same VCF
    LOCK="$CACHE_DIR/.reheader_${chr}"
    if mkdir "$LOCK" 2>/dev/null; then
        echo "  Renaming samples"
        bcftools reheader --samples \
            <(bcftools query -l "$VCF_FILE" \
                | awk '{o=$0; sub(/^[^_]+_/,"",$0); print o,$0}') \
            "$VCF_FILE" -o "${VCF_FILE}.tmp"
        mv "${VCF_FILE}.tmp" "$VCF_FILE"
        tabix -f -p vcf "$VCF_FILE"
        rmdir "$LOCK"
    else
        # Another job is handling the reheader; wait for it
        while [ -d "$LOCK" ]; do sleep 5; done
    fi
else
    echo "  Sample names already clean"
fi

# ── Ancestry mapping (cached, one-off per release) ───────────────────
if [ ! -f "$CACHE_DIR/sample_ancestry.txt" ]; then
    FIRST_VCF=$(find "${IMPUTED_DIRS[@]}" \
        -name "*.vcf.gz" ! -name "*.tbi" 2>/dev/null | head -1)

    first_samp=$(bcftools query -l "$FIRST_VCF" | head -1) || true
    if echo "$first_samp" | grep -q '_'; then
        bcftools reheader --samples \
            <(bcftools query -l "$FIRST_VCF" \
                | awk '{o=$0; sub(/^[^_]+_/,"",$0); print o,$0}') \
            "$FIRST_VCF" -o "${FIRST_VCF}.tmp"
        mv "${FIRST_VCF}.tmp" "$FIRST_VCF"
        tabix -f -p vcf "$FIRST_VCF"
    fi

    bcftools query -l "$FIRST_VCF" > "$CACHE_DIR/sample_order.txt"

    # Create typed sample list: <clean_sid> <type>  (type = c or m)
    # Used for ancestry mapping: C-type → child race (ses-V02),
    # M-type → mother race (ses-V01)
    awk '{type=tolower(substr($0, length($0))); \
          clean=substr($0, 1, length($0)-1); \
          print tolower(clean), type}' "$CACHE_DIR/sample_order.txt" \
        > "$CACHE_DIR/sample_order_typed.txt"

    conda run -n R Rscript --vanilla -e '
        library(tidyverse, quietly = TRUE)
        si <- read_table("'"$CACHE_DIR/sample_order_typed.txt"'",
                         col_names = c("sid", "type"), col_types = "cc")
        dm <- read_delim("'"$ANCESTRY_FILE"'",
                         show_col_types = FALSE, delim = "\t")

        cat("--- DEBUG: dm columns:", paste(names(dm), collapse = ", "), "\n")
        cat("--- DEBUG: session_id values:", paste(unique(dm$session_id), collapse = ", "), "\n")
        cat("--- DEBUG: si head:", nrow(si), "rows, types:", paste(unique(si$type), collapse = ","), "\n")
        print(head(si, 3))

        race_names <- c("0" = "WHT", "1" = "BLK", "2" = "AIAN",
                        "3" = "ASN", "4" = "HPI", "5" = "2PLUS",
                        "6" = "OTH", "7" = "UNK")

        child_lookup <- dm %>%
            filter(!is.na(.data[["sed_basic_demographics_child_race"]]),
                   session_id == "ses-V02") %>%
            count(participant_id, sed_basic_demographics_child_race) %>%
            group_by(participant_id) %>%
            slice_max(n, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            transmute(sid = tolower(participant_id),
                      ancestry = race_names[as.character(sed_basic_demographics_child_race)])

        mother_lookup <- dm %>%
            filter(!is.na(.data[["sed_basic_demographics_screen_mother_race"]]),
                   session_id == "ses-V01") %>%
            count(participant_id, sed_basic_demographics_screen_mother_race) %>%
            group_by(participant_id) %>%
            slice_max(n, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            transmute(sid = tolower(participant_id),
                      ancestry = race_names[as.character(sed_basic_demographics_screen_mother_race)])

        cat("--- DEBUG: child_lookup:", nrow(child_lookup), "rows ---\n")
        print(head(child_lookup, 3))
        cat("--- DEBUG: mother_lookup:", nrow(mother_lookup), "rows ---\n")
        print(head(mother_lookup, 3))

        mapping <- si %>%
            mutate(ancestry = case_when(
                type == "c" ~ child_lookup$ancestry[match(sid, child_lookup$sid)],
                type == "m" ~ mother_lookup$ancestry[match(sid, mother_lookup$sid)],
                TRUE ~ NA_character_
            ), ancestry = replace_na(ancestry, "UNKNOWN"))

        cat("--- DEBUG: mapping ancestry values:",
            paste(unique(mapping$ancestry), collapse = ","), "\n")

        write_tsv(mapping %>% transmute(pos = row_number(), ancestry),
                  file.path("'"$CACHE_DIR"'", "sample_ancestry.txt"),
                  col_names = FALSE)
    '
fi

n_unknown=$(awk '$2 == "UNKNOWN" {c++} END {print c+0}' \
    "$CACHE_DIR/sample_ancestry.txt")
n_total=$(wc -l < "$CACHE_DIR/sample_ancestry.txt")
echo "  Ancestry: $n_total samples, $n_unknown UNKNOWN"

# ── Per-ancestry sample lists ────────────────────────────────────────
ANCESTRY_FILE_USE="$CACHE_DIR/sample_ancestry.txt"
if [ -n "$ANC" ]; then
    ANC_FILE="$CACHE_DIR/anc_${ANC}_${chr}.txt"
    SAMPLE_FILE="$CACHE_DIR/samples_${ANC}_${chr}.txt"

    paste "$CACHE_DIR/sample_order.txt" "$CACHE_DIR/sample_ancestry.txt" \
        | awk -v anc="$ANC" '$3 == anc {print $1}' \
        > "$SAMPLE_FILE"
    paste "$CACHE_DIR/sample_order.txt" "$CACHE_DIR/sample_ancestry.txt" \
        | awk -v anc="$ANC" '$3 == anc {pos++; print pos, $3}' \
        > "$ANC_FILE"

    ANCESTRY_FILE_USE="$ANC_FILE"
    echo "  Filtering to $ANC: $(wc -l < "$SAMPLE_FILE") samples"
fi

# ── Count rare variants ──────────────────────────────────────────────
if [ -z "$INFO_FILE" ]; then
    n_rare=$(bcftools view -i 'R2<0.1 && R2>0' "$VCF_FILE" 2>/dev/null \
        | bcftools query -f '%CHROM\t%POS\n' | wc -l)
else
    n_rare=$(zcat "$INFO_FILE" \
        | awk -F'\t' 'NR>1 && $4!="" && $4+0<0.1 && $4+0>0 {c++} END {print c+0}')
fi

echo "  Rare variants: $n_rare"

OUT_RDS="$CACHE_DIR/chr${chr}${ANC:+_$ANC}.rds"

if [ "$n_rare" -eq 0 ]; then
    echo "  No rare variants, saving empty RDS"
    conda run -n R Rscript --vanilla -e '
        saveRDS(list(sums = NULL, hist_bins = NULL, ridge = NULL, call_rate = NULL),
                file = "'"$OUT_RDS"'")
    '
    exit 0
fi

# ── Main processing ──────────────────────────────────────────────────
TMPDIR=$(mktemp -d)

# Assemble bcftools filter arguments
BC_ARGS=()
[ -n "$ANC" ] && BC_ARGS+=(-S "$SAMPLE_FILE")

# Try VCF INFO/R2 filter (fastest); fall back to info file variant list
if bcftools view -h "$VCF_FILE" 2>/dev/null | grep -q '##INFO.*<ID=R2'; then
    echo "  Filtering: VCF INFO/R2"
    bcftools view -i 'R2<0.1 && R2>0' "${BC_ARGS[@]}" "$VCF_FILE" -Ou 2>/dev/null \
        | bcftools query -f '%CHROM\t%POS[\t%GP]\n' \
        | python "$SCRIPT_DIR/compute_gp_stats.py" \
            --ancestry "$ANCESTRY_FILE_USE" \
            --chr "$chr" \
            --outdir "$TMPDIR"
else
    echo "  Filtering: info file variant list"
    zcat "$INFO_FILE" \
        | awk -F'\t' 'NR>1 && $4!="" && $4+0<0.1 && $4+0>0 {print $1}' \
        > "$CACHE_DIR/chr${chr}_ids.txt"
    bcftools view -i "ID=@$CACHE_DIR/chr${chr}_ids.txt" \
        "${BC_ARGS[@]}" "$VCF_FILE" -Ou 2>/dev/null \
        | bcftools query -f '%CHROM\t%POS[\t%GP]\n' \
        | python "$SCRIPT_DIR/compute_gp_stats.py" \
            --ancestry "$ANCESTRY_FILE_USE" \
            --chr "$chr" \
            --outdir "$TMPDIR"
    rm -f "$CACHE_DIR/chr${chr}_ids.txt"
fi

# ── Assemble RDS ─────────────────────────────────────────────────────
conda run -n R Rscript --vanilla -e '
    library(tidyverse, quietly = TRUE)
    td <- "'"$TMPDIR"'"
    per_sample <- tryCatch(
        read_tsv(file.path(td, "per_sample.tsv"), show_col_types = FALSE),
        error = function(e) NULL)
    hist_bins <- tryCatch(
        read_tsv(file.path(td, "hist_bins.tsv"), show_col_types = FALSE),
        error = function(e) NULL)
    ridge <- tryCatch(
        read_tsv(file.path(td, "ridge.tsv"), show_col_types = FALSE),
        error = function(e) data.frame())
    call_rate <- tryCatch(
        read_tsv(file.path(td, "call_rate.tsv"), show_col_types = FALSE),
        error = function(e) data.frame())
    saveRDS(list(sums = per_sample, hist_bins = hist_bins,
                 ridge = ridge, call_rate = call_rate),
            file = "'"$OUT_RDS"'")
    cat("  Saved", basename("'"$OUT_RDS"'"), "—",
        if (!is.null(per_sample)) nrow(per_sample) else 0, "samples,",
        nrow(ridge), "ridge,", nrow(call_rate), "call-rate\n")
'

# ── Save INFO stats ──────────────────────────────────────────────────
if [ -n "$INFO_FILE" ]; then
    conda run -n R Rscript --vanilla -e '
        library(tidyverse, quietly = TRUE)
        info <- read_tsv("'"$INFO_FILE"'", show_col_types = FALSE)
        if ("Rsq" %in% names(info) && "MAF" %in% names(info)) {
            info <- info %>%
                mutate(
                    chr = str_extract(SNP, "^[^:]+"),
                    pos = as.numeric(str_extract(SNP, "(?<=:)\\d+"))
                ) %>%
                select(chr, pos, maf = MAF, rsq = Rsq)
            saveRDS(info, file.path("'"$CACHE_DIR"'",
                                    "'"chr${chr}_info.rds"'"))
            cat("  Saved chr", "'"${chr}"'", "_info.rds —",
                nrow(info), "variants\n")
        }
    '
fi

rm -rf "$TMPDIR"
echo "--- Finished chr${chr} ${MODE_TAG}---"
