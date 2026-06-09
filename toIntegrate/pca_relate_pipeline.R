#!/usr/bin/env Rscript
# pca_relate_pipeline.R
# PC-Relate pipeline for hbcd_rsid_harmonized
# Uses PC-AiR output from pca_ir_pipeline.R
# LD pruned, 20 PCs, 32 cores

suppressPackageStartupMessages({
  library(gdsfmt)
  library(SNPRelate)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
  library(BiocParallel)
})

# ---------------------------------------------------------------------------
# Hardcoded paths
# ---------------------------------------------------------------------------
WORK         <- "/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026"
NAME         <- "hbcd_rsid_harmonized"
N_PCS        <- 20L
N_CORES      <- 32L
VARIANT_BLOCK <- 50000L

gds_file     <- file.path(WORK, "PCA/gds", paste0(NAME, ".gds"))
seq_gds_file <- sub("\\.gds$", "_seq.gds", gds_file)
pcair_rds    <- file.path(WORK, "PCA/pca_ir", paste0(NAME, "_pcaobj.RDS"))
out_dir      <- file.path(WORK, "PCA/pca_ir")
out_rds      <- file.path(out_dir, paste0(NAME, "_pcrelate.RDS"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# MHC region to exclude (GRCh38)
drop_MHC  <- TRUE
mhc_chr   <- "6"
mhc_start <- 25e6
mhc_end   <- 34e6

# ---------------------------------------------------------------------------
# Step 0 — Close any open GDS handles
# ---------------------------------------------------------------------------
gdsfmt::showfile.gds(closeall = TRUE)

# ---------------------------------------------------------------------------
# Step 1 — LD pruning on SNPRelate GDS
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 1: LD pruning\n")

if (!file.exists(gds_file)) stop("SNPRelate GDS not found: ", gds_file)

genofile <- snpgdsOpen(gds_file)

# Get all variant info for MHC exclusion
all_chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
all_pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
all_snp <- read.gdsn(index.gdsn(genofile, "snp.id"))

# Exclude MHC region
if (drop_MHC) {
  chr_char <- as.character(all_chr)
  mhc_mask <- chr_char == mhc_chr & all_pos >= mhc_start & all_pos <= mhc_end
  snps_no_mhc <- all_snp[!mhc_mask]
  cat("  Variants after MHC drop:", length(snps_no_mhc), "\n")
} else {
  snps_no_mhc <- all_snp
}

# LD pruning: window 500kb, r^2 threshold 0.1
set.seed(42)
pruned_snps <- snpgdsLDpruning(
  genofile,
  snp.id        = snps_no_mhc,
  autosome.only = TRUE,
  ld.threshold  = sqrt(0.1),   # r threshold = sqrt(r^2 threshold)
  slide.max.bp  = 500000L,
  verbose       = TRUE
)

pruned_ids <- unlist(pruned_snps, use.names = FALSE)
cat("  Variants after LD pruning:", length(pruned_ids), "\n")
snpgdsClose(genofile)

# ---------------------------------------------------------------------------
# Step 2 — Open SeqArray GDS and apply filters
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 2: Opening SeqArray GDS\n")

if (!file.exists(seq_gds_file)) stop("SeqArray GDS not found: ", seq_gds_file)

seqfile <- seqOpen(seq_gds_file, allow.duplicate = TRUE)
on.exit({ try(seqClose(seqfile), silent = TRUE) }, add = TRUE)

seqResetFilter(seqfile)
all_samples  <- seqGetData(seqfile, "sample.id")
all_seq_vars <- seqGetData(seqfile, "variant.id")

cat("  SeqArray samples:", length(all_samples), "\n")
cat("  SeqArray variants:", length(all_seq_vars), "\n")

# Map pruned SNP IDs to SeqArray variant IDs
# Both GDS files are derived from the same PLINK file so positions match
seq_chr <- seqGetData(seqfile, "chromosome")
seq_pos <- seqGetData(seqfile, "position")

# Build position key for matching
snp_pos_key <- paste(all_chr[match(pruned_ids, all_snp)],
                     all_pos[match(pruned_ids, all_snp)], sep = ":")
seq_pos_key <- paste(seq_chr, seq_pos, sep = ":")
keep_seq_vars <- all_seq_vars[seq_pos_key %in% snp_pos_key]

cat("  SeqArray variants after LD prune mapping:", length(keep_seq_vars), "\n")

if (length(keep_seq_vars) < 1000) {
  stop("Fewer than 1000 variants after pruning — check GDS position matching")
}

seqSetFilter(seqfile, variant.id = keep_seq_vars, verbose = FALSE)

# ---------------------------------------------------------------------------
# Step 3 — Load PC-AiR object and extract first 20 PCs
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 3: Loading PC-AiR PCs\n")

if (!file.exists(pcair_rds)) stop("PC-AiR RDS not found: ", pcair_rds)
pcair_obj <- readRDS(pcair_rds)

pcs_all <- pcair_obj$vectors
if (is.null(pcs_all)) pcs_all <- pcair_obj$eigenvectors
if (is.null(pcs_all)) stop("Cannot find eigenvectors in PC-AiR object")

# Use first N_PCS only
pcs_mat <- pcs_all[, seq_len(min(N_PCS, ncol(pcs_all))), drop = FALSE]
cat("  Using", ncol(pcs_mat), "PCs from PC-AiR object\n")

# Align to active sample order in SeqArray
active_samples <- seqGetData(seqfile, "sample.id")
match_idx <- match(active_samples, rownames(pcs_mat))
if (any(is.na(match_idx))) {
  stop("Some GDS sample IDs not found in PC-AiR eigenvectors")
}
pcs_mat <- pcs_mat[match_idx, , drop = FALSE]

# Training set = unrelated samples
unrels       <- pcair_obj$unrels
training_set <- intersect(unrels, active_samples)
cat("  Training set (unrelateds):", length(training_set), "\n")

# ---------------------------------------------------------------------------
# Step 4 — Build iterator and run PC-Relate
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 4: Running PC-Relate\n")
cat("  Samples:", length(active_samples), "\n")
cat("  Variants:", length(keep_seq_vars), "\n")
cat("  Cores:", N_CORES, "\n")

seqData <- SeqVarData(seqfile)
seqIter <- SeqVarBlockIterator(seqData, variantBlock = VARIANT_BLOCK, verbose = FALSE)

BPPARAM <- MulticoreParam(workers = N_CORES, progressbar = TRUE)

relate <- pcrelate(
  seqIter,
  pcs                = pcs_mat,
  training.set       = training_set,
  ibd.probs          = FALSE,
  scale              = "variant",
  small.samp.correct = TRUE,
  BPPARAM            = BPPARAM,
  verbose            = TRUE
)

# ---------------------------------------------------------------------------
# Step 5 — Save results
# ---------------------------------------------------------------------------
cat("[", format(Sys.time()), "] Step 5: Saving results\n")

saveRDS(relate, file = out_rds)
cat("  Saved:", out_rds, "\n")

# Also save kinship pairs as a readable table
kin_pairs <- pcrelateToMatrix(relate, thresh = 2^(-9/2), scaleKin = 2)
saveRDS(kin_pairs, file = file.path(out_dir, paste0(NAME, "_pcrelate_kinmat.RDS")))

cat("[", format(Sys.time()), "] PC-Relate pipeline complete.\n")
