#!/bin/bash
#SBATCH --job-name=gp_aggregate
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/gp_aggregate_%j.out
#SBATCH --error=logs/gp_aggregate_%j.err

# step 26: combine per-chromosome GP results from step 25 into final CSVs
# Run AFTER all array tasks in 25-precompute_gp.SLURM complete:
#   sbatch --dependency=afterok:<jobid> 26-aggregate_gp.SLURM

set -euo pipefail

IMPUTATION_DIR="${HBCD_IMPUTATION_DIR:-/projects/standard/basu_hbcd/shared/hbcdSandboxData}"
CACHE_DIR="${IMPUTATION_DIR}/imputed/gp_cache"
IMPDIR="${IMPUTATION_DIR}/imputed"
IMPUTED_DIRS=("${IMPDIR}/c1" "${IMPDIR}/c2" "${IMPDIR}/c3")

echo "--- Aggregating per-chromosome RDS files from $CACHE_DIR ---"

Rscript --vanilla -e '
library(tidyverse, quietly = TRUE)
library(data.table)

rds_files <- list.files("'"$CACHE_DIR"'", pattern = "^chr.*\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "chromosome files\n")

all_sums     <- list()
all_hist     <- list()
all_ridge    <- list()
all_callrate <- list()

for (f in rds_files) {
  dat <- readRDS(f)
  if (!is.null(dat$sums) && nrow(dat$sums) > 0) all_sums[[length(all_sums) + 1]] <- dat$sums
  if (!is.null(dat$hist_bins) && nrow(dat$hist_bins) > 0) all_hist[[length(all_hist) + 1]] <- dat$hist_bins
  if (!is.null(dat$ridge) && nrow(dat$ridge) > 0) all_ridge[[length(all_ridge) + 1]] <- dat$ridge
  if (!is.null(dat$call_rate) && nrow(dat$call_rate) > 0) all_callrate[[length(all_callrate) + 1]] <- dat$call_rate
}

per_sample <- bind_rows(all_sums) %>%
  group_by(sample_pos, ancestry) %>%
  summarise(
    avg_maxgp    = sum(sum_maxgp) / sum(n_variants),
    n_variants   = sum(n_variants),
    n_confident  = sum(n_confident),
    pct_confident = n_confident / sum(n_variants) * 100,
    .groups = "drop"
  ) %>%
  arrange(sample_pos)

write_csv(per_sample, file.path("'"$CACHE_DIR"'", "per_sample_maxgp.csv.gz"))
cat("Per-sample results:", nrow(per_sample), "samples\n")

hist_bins <- bind_rows(all_hist) %>%
  group_by(ancestry, chr, bin) %>%
  summarise(count = sum(count), .groups = "drop")

write_csv(hist_bins, file.path("'"$CACHE_DIR"'", "histogram_bins.csv.gz"))
cat("Histogram bins:", nrow(hist_bins), "rows\n")

ridge_raw <- bind_rows(all_ridge)
if (nrow(ridge_raw) > 0) {
  ridge_raw <- ridge_raw %>% mutate(chromosome = fct_inorder(chromosome))
}
write_csv(ridge_raw, file.path("'"$CACHE_DIR"'", "ridge_samples.csv.gz"))
cat("Ridge samples:", nrow(ridge_raw), "rows\n")

cr <- bind_rows(all_callrate)
if (nrow(cr) > 0) {
  info_files <- list.files("'"$CACHE_DIR"'", pattern = "chr.*_info\\.rds$", full.names = TRUE)
  info_data <- map_dfr(info_files, readRDS)

  call_rate_all <- cr %>%
    left_join(info_data, by = c("chr", "pos")) %>%
    mutate(call_rate = called / total) %>%
    filter(!is.na(rsq), !is.na(maf))
} else {
  call_rate_all <- data.frame()
}
write_csv(call_rate_all, file.path("'"$CACHE_DIR"'", "call_rate_by_ancestry.csv.gz"))
cat("Call-rate entries:", nrow(call_rate_all), "\n")

imputed_dirs <- c(file.path("'"$IMPUTATION_DIR"'", "imputed", "c1"),
                  file.path("'"$IMPUTATION_DIR"'", "imputed", "c2"),
                  file.path("'"$IMPUTATION_DIR"'", "imputed", "c3"))
info_gz_files <- list.files(imputed_dirs, pattern = "\\.info\\.gz$", recursive = TRUE, full.names = TRUE)
if (length(info_gz_files) > 0) {
  yield <- map_dfr(info_gz_files, ~{
    chr <- str_extract(basename(.x), "\\d+|X|Y")
    fread(.x, select = c("Rsq")) %>%
      summarise(
        chr = chr,
        total     = n(),
        gt_03     = sum(Rsq > 0.3, na.rm = TRUE),
        gt_05     = sum(Rsq > 0.5, na.rm = TRUE),
        gt_08     = sum(Rsq > 0.8, na.rm = TRUE),
        gt_09     = sum(Rsq > 0.9, na.rm = TRUE)
      )
  })
  write_csv(yield, file.path("'"$CACHE_DIR"'", "imputation_yield.csv.gz"))
  cat("Imputation yield:", nrow(yield), "chromosomes\n")
} else {
  cat("No info files found for yield calculation\n")
}

cat("--- All done! CSVs in", "'"$CACHE_DIR"' ---\n")
'
