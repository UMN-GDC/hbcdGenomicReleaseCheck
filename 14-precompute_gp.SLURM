#!/bin/bash
#SBATCH --job-name=gp_prep
#SBATCH --array=0-23
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/gp_prep_%A_%a.out
#SBATCH --error=logs/gp_prep_%A_%a.err

# step 25: precompute genotype probabilities for rare imputed variants
# Per-chromosome: renames VCF samples, maps ancestry, extracts GP for rare variants

set -euo pipefail

module load bcftools
module load rstudio/4.5.3

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
chr="${CHROMS[$SLURM_ARRAY_TASK_ID]}"
IMPDIR=/home/ood-coffm049/hbcdData/imputed
ANCESTRY_FILE="/shared/release/hbcd/hbcd/rawdata/phenotype/sed_bm_demo.tsv"
IMPUTED_DIRS=($IMPDIR/c1 $IMPDIR/c2 $IMPDIR/c3)
CACHE_DIR=$IMPDIR/gp_cache
mkdir -p "$CACHE_DIR" logs

VCF_FILE=$(find "${IMPUTED_DIRS[@]}" -name "chr${chr}.dose.vcf.gz" ! -name "*.tbi" 2>/dev/null | head -1)
INFO_FILE=$(find "${IMPUTED_DIRS[@]}" -name "chr${chr}.info.gz" 2>/dev/null | head -1)

if [ -z "$VCF_FILE" ] || [ -z "$INFO_FILE" ]; then
  echo "Skipping chr${chr}: VCF or INFO not found"
  exit 0
fi

echo "--- Processing chr${chr} ---"

first_sample=$(bcftools query -l "$VCF_FILE" | head -1) || true
if echo "$first_sample" | grep -q '_'; then
  echo "  Renaming samples (removing FID_ prefix)"
  bcftools reheader --samples \
    <(bcftools query -l "$VCF_FILE" | awk '{o=$0; sub(/^[^_]+_/,"",$0); print o,$0}') \
    "$VCF_FILE" | bcftools view -Oz -o "${VCF_FILE}.tmp"
  mv "${VCF_FILE}.tmp" "$VCF_FILE"
  tabix -f -p vcf "$VCF_FILE"
else
  echo "  Sample names already clean"
fi

if [ ! -f "$CACHE_DIR/sample_ancestry.txt" ]; then
  FIRST_VCF=$(find "${IMPUTED_DIRS[@]}" -name "*.vcf.gz" ! -name "*.tbi" 2>/dev/null | head -1)

  first_samp=$(bcftools query -l "$FIRST_VCF" | head -1) || true
  if echo "$first_samp" | grep -q '_'; then
    bcftools reheader --samples \
      <(bcftools query -l "$FIRST_VCF" | awk '{o=$0; sub(/^[^_]+_/,"",$0); print o,$0}') \
      "$FIRST_VCF" | bcftools view -Oz -o "${FIRST_VCF}.tmp"
    mv "${FIRST_VCF}.tmp" "$FIRST_VCF"
    tabix -f -p vcf "$FIRST_VCF"
  fi

  bcftools query -l "$FIRST_VCF" > "$CACHE_DIR/sample_order.txt"
  Rscript --vanilla -e '
    library(tidyverse, quietly = TRUE)
    samples <- read_lines("'"$CACHE_DIR/sample_order.txt"'")

    ancestry <- read_delim("'"$ANCESTRY_FILE"'", show_col_types = FALSE, delim = "\t")
    ancestry <- ancestry %>%
      filter(!is.na(.data[["sed_basic_demographics_child_race"]])) %>%
      count(participant_id, sed_basic_demographics_child_race) %>%
      group_by(participant_id) %>%
      slice_max(n, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(sid = participant_id, ancestry = sed_basic_demographics_child_race)

    mapping <- tibble(sid = tolower(samples)) %>%
      left_join(ancestry %>% mutate(sid = tolower(sid)), by = "sid")
    mapping$ancestry[is.na(mapping$ancestry)] <- "UNKNOWN"
    write_tsv(mutate(mapping, pos = row_number()) %>% select(pos, ancestry),
              file.path("'"$CACHE_DIR"'", "sample_ancestry.txt"), col_names = FALSE)
  '
fi

if [ -f "$CACHE_DIR/sample_ancestry.txt" ]; then
  n_unknown=$(awk '$2 == "UNKNOWN" {c++} END {print c+0}' "$CACHE_DIR/sample_ancestry.txt")
  n_total=$(wc -l < "$CACHE_DIR/sample_ancestry.txt")
  echo "  Ancestry mapping: $n_total samples, $n_unknown UNKNOWN"
  if [ "$n_unknown" -gt "$((n_total / 2))" ]; then
    echo "WARNING: >50% of samples are UNKNOWN — check ANCESTRY_FILE path and sample ID match"
  fi
fi

zcat "$INFO_FILE" | awk -F'\t' '
  NR == 1 { next }
  $4 != "" && $4+0 < 0.1 && $4+0 > 0 {
    split($1, a, ":"); pos = a[2]
    if (pos == "") { split($1, a, "_"); pos = a[2] }
    if (pos != "") print a[1] ":" pos
  }
' > "$CACHE_DIR/chr${chr}_regions.txt" 2>/dev/null || true

n_rare=$(wc -l < "$CACHE_DIR/chr${chr}_regions.txt" 2>/dev/null || echo 0)
echo "  Rare variants: $n_rare"

if [ "$n_rare" -eq 0 ]; then
  echo "  No rare variants, saving empty"
  Rscript --vanilla -e '
    saveRDS(list(sums = NULL, hist_bins = NULL, ridge = NULL),
            file = "'"$CACHE_DIR/chr${chr}.rds"'")
  '
  rm -f "$CACHE_DIR/chr${chr}_regions.txt"
  exit 0
fi

TMPDIR=$(mktemp -d)

bcftools query -R "$CACHE_DIR/chr${chr}_regions.txt" \
  -f '%CHROM\t%POS[\t%GP]\n' "$VCF_FILE" 2>/dev/null | \
gawk -F'\t' -v chr="$chr" -v od="$CACHE_DIR" -v td="$TMPDIR" '
BEGIN {
  while ((getline < (od "/sample_ancestry.txt")) > 0) anc[$1] = $2
  close(od "/sample_ancestry.txt")
  srand(42)
  ridge_rate = 0.002
  call_rate_rate = 0.01
  n_samples = 0
}
{
  if (n_samples == 0) n_samples = NF - 2
  delete anc_tot; delete anc_call

  for (i = 3; i <= NF; i++) {
    split($i, gp, ",")
    maxgp = gp[1] > gp[2] ? (gp[1] > gp[3] ? gp[1] : gp[3]) : (gp[2] > gp[3] ? gp[2] : gp[3])
    a = anc[i]; if (a == "") a = "UNKNOWN"

    sum_gp[i] += maxgp; cnt_gp[i]++
    if (maxgp > 0.9) conf_gp[i]++

    bin = int(maxgp * 50) + 1; if (bin > 50) bin = 50
    h[a, bin]++

    anc_tot[a]++
    if (maxgp > 0.9) anc_call[a]++

    if (rand() < ridge_rate) print chr, a, maxgp >> (td "/ridge.txt")
  }

  if (rand() < call_rate_rate) {
    pos = $2
    for (a in anc_tot)
      print chr, pos, a, anc_call[a], anc_tot[a] >> (td "/call_rate.txt")
  }
}
END {
  for (i = 3; i <= NF; i++) {
    a = anc[i]; if (a == "") a = "UNKNOWN"
    printf "%d\t%s\t%s\t%.10f\t%d\t%d\n", i-2, chr, a, sum_gp[i], cnt_gp[i], conf_gp[i] >> (td "/per_sample.txt")
  }
  for (k in h) { split(k, p, SUBSEP); print p[1], chr, p[2], h[k] >> (td "/hist_bins.txt") }
}'

Rscript --vanilla -e '
  library(tidyverse, quietly = TRUE)

  per_sample <- read_tsv(file.path("'"$TMPDIR"'", "per_sample.txt"),
                         col_names = c("sample_pos", "chr", "ancestry", "sum_maxgp", "n_variants", "n_confident"),
                         col_types = "iccdii", show_col_types = FALSE)

  hist_bins <- read_tsv(file.path("'"$TMPDIR"'", "hist_bins.txt"),
                        col_names = c("ancestry", "chr", "bin", "count"),
                        col_types = "ccii", show_col_types = FALSE)

  ridge <- data.frame()
  ridge_file <- file.path("'"$TMPDIR"'", "ridge.txt")
  if (file.exists(ridge_file) && file.size(ridge_file) > 0) {
    ridge <- read_tsv(ridge_file, col_names = c("chromosome", "ancestry", "maxgp"),
                      col_types = "ccd", show_col_types = FALSE)
  }

  call_rate <- data.frame()
  cr_file <- file.path("'"$TMPDIR"'", "call_rate.txt")
  if (file.exists(cr_file) && file.size(cr_file) > 0) {
    call_rate <- read_tsv(cr_file, col_names = c("chr", "pos", "ancestry", "called", "total"),
                          col_types = "cicii", show_col_types = FALSE)
  }

  info_cols <- tryCatch({
    info <- read_tsv(file.path("'"$INFO_FILE"'"), show_col_types = FALSE)
    names(info)
  }, error = function(e) NULL)
  if ("Rsq" %in% info_cols && "MAF" %in% info_cols) {
    info <- read_tsv(file.path("'"$INFO_FILE"'"), show_col_types = FALSE) %>%
      mutate(
        chr = str_extract(SNP, "^[^:]+"),
        pos = as.numeric(str_extract(SNP, "(?<=:)\\d+"))
      ) %>%
      select(chr, pos, maf = MAF, rsq = Rsq)
    saveRDS(info, file.path("'"$CACHE_DIR"'", "'"chr${chr}_info.rds"'"))
  }

  saveRDS(list(sums = per_sample, hist_bins = hist_bins, ridge = ridge, call_rate = call_rate),
          file = "'"$CACHE_DIR/chr${chr}.rds"'")
  cat("  Saved chr${chr}.rds —", nrow(per_sample), "samples,",
      nrow(ridge), "ridge samples,", nrow(call_rate), "call-rate samples\n")
'

rm -rf "$TMPDIR" "$CACHE_DIR/chr${chr}_regions.txt"
echo "--- Finished chr${chr} ---"
