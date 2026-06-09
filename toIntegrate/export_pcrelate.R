#!/usr/bin/env Rscript
# export_pcrelate.R
# Loads PC-Relate RDS outputs and exports to CSV and TSV

library(GENESIS)
library(Matrix)

WORK    <- "/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026"
out_dir <- file.path(WORK, "PCA/pca_ir")

# ---------------------------------------------------------------------------
# Helper to write both CSV and TSV
# ---------------------------------------------------------------------------
save_both <- function(df, stem) {
  write.csv(df,
            file      = paste0(stem, ".csv"),
            row.names = FALSE,
            quote     = FALSE)
  write.table(df,
              file      = paste0(stem, ".tsv"),
              sep       = "\t",
              row.names = FALSE,
              quote     = FALSE)
  cat("  Saved:", basename(paste0(stem, ".csv")), "\n")
  cat("  Saved:", basename(paste0(stem, ".tsv")), "\n")
}

# ---------------------------------------------------------------------------
# 1. PC-Relate pairwise kinship table
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Loading pcrelate.RDS\n")
relate <- readRDS(file.path(out_dir, "hbcd_rsid_harmonized_pcrelate.RDS"))

# Extract pairwise results table
kin_pairs <- relate$kinBtwn
cat("  Pairwise kinship pairs:", nrow(kin_pairs), "\n")
cat("  Columns:", paste(colnames(kin_pairs), collapse=", "), "\n")

save_both(kin_pairs,
          file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_pairs"))

# Also export within-individual estimates if present
if (!is.null(relate$kinSelf)) {
  kin_self <- relate$kinSelf
  cat("  Within-individual rows:", nrow(kin_self), "\n")
  save_both(kin_self,
            file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_self"))
}

# ---------------------------------------------------------------------------
# 2. Kinship matrix (sparse -> dense -> long format)
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Loading pcrelate_kinmat.RDS\n")
kinmat <- readRDS(file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_kinmat.RDS"))

# Convert sparse matrix to dense
kinmat_dense <- as.matrix(kinmat)
cat("  Kinship matrix dimensions:", dim(kinmat_dense), "\n")

# Save as wide matrix (rows = samples, columns = samples)
kinmat_df <- as.data.frame(kinmat_dense)
kinmat_df <- cbind(SampleID = rownames(kinmat_dense), kinmat_df)

save_both(kinmat_df,
          file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_kinmat_wide"))

# Also save as long format (ID1, ID2, kinship) — more portable for large matrices
cat("[", format(Sys.time()), "] Converting kinship matrix to long format\n")
long_df <- data.frame(
  ID1     = rownames(kinmat_dense)[row(kinmat_dense)],
  ID2     = colnames(kinmat_dense)[col(kinmat_dense)],
  Kinship = as.vector(kinmat_dense)
)
# Keep only upper triangle + diagonal to avoid duplication
long_df <- long_df[row(kinmat_dense) <= col(kinmat_dense), ]
long_df <- long_df[long_df$Kinship != 0, ]   # drop structural zeros
cat("  Long format rows (non-zero):", nrow(long_df), "\n")

save_both(long_df,
          file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_kinmat_long"))

cat("[", format(Sys.time()), "] All exports complete.\n")
cat("\nFiles written to:", out_dir, "\n")

# ---------------------------------------------------------------------------
# 3. PC-AiR results
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Loading pcaobj.RDS\n")
pcaobj <- readRDS(file.path(out_dir, "hbcd_rsid_harmonized_pcaobj.RDS"))

# PC scores (eigenvectors) — one row per sample, one column per PC
pcs <- as.data.frame(pcaobj$vectors)
colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
pcs <- cbind(SampleID = rownames(pcs), pcs)

save_both(pcs,
          file.path(out_dir, "hbcd_rsid_harmonized_pcair_scores"))
cat("  PC scores dimensions:", nrow(pcs), "samples x", ncol(pcs)-1, "PCs\n")

# Eigenvalues
evals <- data.frame(
  PC         = paste0("PC", seq_along(pcaobj$values)),
  Eigenvalue = pcaobj$values,
  VarExplained = pcaobj$values / sum(pcaobj$values)
)
save_both(evals,
          file.path(out_dir, "hbcd_rsid_harmonized_pcair_eigenvalues"))

# Unrelated sample IDs
unrels_df <- data.frame(SampleID = pcaobj$unrels)
save_both(unrels_df,
          file.path(out_dir, "hbcd_rsid_harmonized_pcair_unrelated_ids"))

# Related sample IDs
rels_df <- data.frame(SampleID = pcaobj$rels)
save_both(rels_df,
          file.path(out_dir, "hbcd_rsid_harmonized_pcair_related_ids"))

cat("[", format(Sys.time()), "] PC-AiR exports complete.\n")
cat("\nAll files written to:", out_dir, "\n")
