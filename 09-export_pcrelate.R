#!/usr/bin/env Rscript
# step 09: export PC-Relate results to CSV and TSV (including IBD probabilities)

library(GENESIS)
library(Matrix)

DATA_DIR <- "/projects/standard/basu_hbcd/shared/data"
out_dir  <- DATA_DIR

save_both <- function(df, stem) {
  write.csv(df, file = paste0(stem, ".csv"), row.names = FALSE, quote = FALSE)
  write.table(df, file = paste0(stem, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Saved:", basename(paste0(stem, ".csv")), "\n")
  cat("  Saved:", basename(paste0(stem, ".tsv")), "\n")
}

cat("[", format(Sys.time()), "] Loading pcrelate.RDS\n")
relate <- readRDS(file.path(out_dir, "hbcd_rsid_harmonized_pcrelate.RDS"))

kin_pairs <- relate$kinBtwn
cat("  Pairwise kinship pairs:", nrow(kin_pairs), "\n")
cat("  Columns:", paste(colnames(kin_pairs), collapse=", "), "\n")
save_both(kin_pairs, file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_pairs"))

# IBD probabilities (available when ibd.probs = TRUE)
ibd_cols <- intersect(c("ID1", "ID2", "IBD0", "IBD1", "IBD2", "kinship", "k0", "k1"), colnames(kin_pairs))
if (length(ibd_cols) > 3) {
  ibd_pairs <- kin_pairs[, ibd_cols]
  save_both(ibd_pairs, file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_ibd"))
  cat("  IBD columns exported:", paste(setdiff(ibd_cols, c("ID1", "ID2")), collapse=", "), "\n")
}

if (!is.null(relate$kinSelf)) {
  kin_self <- relate$kinSelf
  cat("  Within-individual rows:", nrow(kin_self), "\n")
  save_both(kin_self, file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_self"))
}

cat("[", format(Sys.time()), "] Loading pcrelate_kinmat.RDS\n")
kinmat <- readRDS(file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_kinmat.RDS"))
kinmat_dense <- as.matrix(kinmat)
cat("  Kinship matrix dimensions:", dim(kinmat_dense), "\n")

kinmat_df <- as.data.frame(kinmat_dense)
kinmat_df <- cbind(SampleID = rownames(kinmat_dense), kinmat_df)
save_both(kinmat_df, file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_kinmat_wide"))

cat("[", format(Sys.time()), "] Converting kinship matrix to long format\n")
long_df <- data.frame(
  ID1     = rownames(kinmat_dense)[row(kinmat_dense)],
  ID2     = colnames(kinmat_dense)[col(kinmat_dense)],
  Kinship = as.vector(kinmat_dense)
)
long_df <- long_df[row(kinmat_dense) <= col(kinmat_dense), ]
long_df <- long_df[long_df$Kinship != 0, ]
cat("  Long format rows (non-zero):", nrow(long_df), "\n")
save_both(long_df, file.path(out_dir, "hbcd_rsid_harmonized_pcrelate_kinmat_long"))

cat("[", format(Sys.time()), "] Loading pcaobj.RDS\n")
pcaobj <- readRDS(file.path(out_dir, "hbcd_rsid_harmonized_pcaobj.RDS"))

pcs <- as.data.frame(pcaobj$vectors)
colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
pcs <- cbind(SampleID = rownames(pcs), pcs)
save_both(pcs, file.path(out_dir, "hbcd_rsid_harmonized_pcair_scores"))
cat("  PC scores dimensions:", nrow(pcs), "samples x", ncol(pcs)-1, "PCs\n")

evals <- data.frame(
  PC            = paste0("PC", seq_along(pcaobj$values)),
  Eigenvalue    = pcaobj$values,
  VarExplained  = pcaobj$values / sum(pcaobj$values)
)
save_both(evals, file.path(out_dir, "hbcd_rsid_harmonized_pcair_eigenvalues"))

unrels_df <- data.frame(SampleID = pcaobj$unrels)
save_both(unrels_df, file.path(out_dir, "hbcd_rsid_harmonized_pcair_unrelated_ids"))

rels_df <- data.frame(SampleID = pcaobj$rels)
save_both(rels_df, file.path(out_dir, "hbcd_rsid_harmonized_pcair_related_ids"))

cat("[", format(Sys.time()), "] All exports complete.\n")
cat("\nFiles written to:", out_dir, "\n")
