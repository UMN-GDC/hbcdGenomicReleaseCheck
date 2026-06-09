#!/usr/bin/env Rscript
# pca_ir_pipeline.R
# PC-AiR pipeline for hbcd_rsid_harmonized
# Uses existing KING output from relatedness step

library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(GENESIS)

# ---------------------------------------------------------------------------
# Hardcoded paths
# ---------------------------------------------------------------------------
WORK         <- "/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026"
NAME         <- "hbcd_rsid_harmonized"
plink_prefix <- file.path(WORK, "AE/rsid_updated", NAME)
king_kin0    <- file.path(WORK, "relatedness", "kinships.kin0")
king_kin     <- file.path(WORK, "relatedness", "kinships.kin")
gds_file     <- file.path(WORK, "PCA/gds", paste0(NAME, ".gds"))
seq_gds_file <- sub("\\.gds$", "_seq.gds", gds_file)
out_dir      <- file.path(WORK, "PCA/pca_ir")

dir.create(dirname(gds_file), showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir,           showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Step 1 — Convert PLINK to SNPRelate GDS
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 1: Converting PLINK to GDS\n")

if (!file.exists(gds_file)) {
  snpgdsBED2GDS(
    bed.fn   = paste0(plink_prefix, ".bed"),
    bim.fn   = paste0(plink_prefix, ".bim"),
    fam.fn   = paste0(plink_prefix, ".fam"),
    out.gdsfn = gds_file
  )
} else {
  cat("  GDS file already exists — skipping conversion\n")
}

# ---------------------------------------------------------------------------
# Step 2 — Open GDS and get sample IDs
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 2: Opening GDS file\n")
genofile   <- snpgdsOpen(gds_file)
sample.id  <- read.gdsn(index.gdsn(genofile, "sample.id"))
cat("  Samples in GDS:", length(sample.id), "\n")

# ---------------------------------------------------------------------------
# Step 3 — Build KING sparse matrix from existing KING output
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 3: Loading KING relatedness output\n")

read_king_file <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read.table(path, header = TRUE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  df
}

kin0 <- read_king_file(king_kin0)
kinw <- read_king_file(king_kin)

# Normalise column names across kin and kin0
norm_king <- function(df, within = FALSE) {
  if (is.null(df)) return(NULL)
  cn <- colnames(df)
  if ("IID1"    %in% cn) names(df)[names(df) == "IID1"]    <- "ID1"
  if ("IID2"    %in% cn) names(df)[names(df) == "IID2"]    <- "ID2"
  if ("KINSHIP" %in% cn) names(df)[names(df) == "KINSHIP"] <- "Kinship"
  # For within-family (kin), FID is a single column
  if (within && "FID" %in% cn) {
    df$FID1 <- df$FID
    df$FID2 <- df$FID
  }
  df
}

kin0 <- norm_king(kin0, within = FALSE)
kinw <- norm_king(kinw, within = TRUE)

all_pairs <- do.call(rbind, Filter(Negate(is.null), list(
  if (!is.null(kin0)) data.frame(ID1 = kin0$ID1, ID2 = kin0$ID2,
                                  Kinship = as.numeric(kin0$Kinship)),
  if (!is.null(kinw)) data.frame(ID1 = kinw$ID1, ID2 = kinw$ID2,
                                  Kinship = as.numeric(kinw$Kinship))
)))

cat("  Total KING pairs loaded:", nrow(all_pairs), "\n")

# Build symmetric sparse kinship matrix
# Only keep pairs where both IDs are in the GDS sample list
keep <- all_pairs$ID1 %in% sample.id & all_pairs$ID2 %in% sample.id
all_pairs <- all_pairs[keep, ]
cat("  Pairs with both IDs in GDS:", nrow(all_pairs), "\n")

# Build sparse symmetric kinship matrix manually from KING pairs
# Matrix package provides sparseMatrix
library(Matrix)

n <- length(sample.id)
id_idx <- setNames(seq_len(n), sample.id)

i_idx <- id_idx[all_pairs$ID1]
j_idx <- id_idx[all_pairs$ID2]

# Remove any pairs where ID lookup failed
valid <- !is.na(i_idx) & !is.na(j_idx)
i_idx <- i_idx[valid]
j_idx <- j_idx[valid]
kvals <- all_pairs$Kinship[valid]

# Build symmetric sparse matrix (upper + lower + diagonal)
king_sparse <- sparseMatrix(
  i    = c(i_idx, j_idx, seq_len(n)),
  j    = c(j_idx, i_idx, seq_len(n)),
  x    = c(kvals, kvals, rep(0.5, n)),
  dims = c(n, n),
  dimnames = list(sample.id, sample.id)
)

cat("  KING sparse matrix dimensions:", dim(king_sparse), "\n")

# ---------------------------------------------------------------------------
# Step 4 — PC-AiR
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 4: Running PC-AiR\n")

pcaobj <- pcair(
  genofile,
  kinobj  = king_sparse,
  divobj  = king_sparse,
  verbose = TRUE
)

cat("  PC-AiR complete\n")
cat("  Unrelated set size:", length(pcaobj$unrels), "\n")
cat("  Related set size  :", length(pcaobj$rels), "\n")

# ---------------------------------------------------------------------------
# Step 5 — Convert to SeqArray GDS for PC-Relate
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 5: Converting to SeqArray GDS\n")

snpgdsClose(genofile)

if (!file.exists(seq_gds_file)) {
  seqSNP2GDS(gds.fn = gds_file, out.fn = seq_gds_file, verbose = TRUE)
} else {
  cat("  SeqArray GDS already exists — skipping\n")
}

seqfile  <- seqOpen(seq_gds_file)
seqData  <- SeqVarData(seqfile)

# ---------------------------------------------------------------------------
# Step 6 — Save results
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 6: Saving results\n")

saveRDS(pcaobj,
        file = file.path(out_dir, paste0(NAME, "_pcaobj.RDS")))

# Save PC scores
pc_scores <- as.data.frame(pcaobj$vectors)
pc_scores <- cbind(sample.id = rownames(pc_scores), pc_scores)
write.table(pc_scores,
            file  = file.path(out_dir, paste0(NAME, "_pc_scores.txt")),
            quote = FALSE, row.names = FALSE, sep = "\t")

# Save unrelated IDs
write.table(pcaobj$unrels,
            file  = file.path(out_dir, paste0(NAME, "_unrelated_ids.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Save related IDs
write.table(pcaobj$rels,
            file  = file.path(out_dir, paste0(NAME, "_related_ids.txt")),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("  Saved: ", NAME, "_pcaobj.RDS\n", sep = "")
cat("  Saved: ", NAME, "_pc_scores.txt\n", sep = "")
cat("  Saved: ", NAME, "_unrelated_ids.txt\n", sep = "")
cat("  Saved: ", NAME, "_related_ids.txt\n", sep = "")

# ---------------------------------------------------------------------------
# Step 7 — Cleanup
# ---------------------------------------------------------------------------
seqClose(seqfile)

cat("[", format(Sys.time()), "] Done.\n")
