#!/usr/bin/env Rscript
# step 11: PC-Relate pipeline for hbcd_rsid_harmonized
# Uses PC-AiR output from step 9 (pca_ir_pipeline.R)

suppressPackageStartupMessages({
  library(gdsfmt)
  library(SNPRelate)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
  library(BiocParallel)
})

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

drop_MHC  <- TRUE
mhc_chr   <- "6"
mhc_start <- 25e6
mhc_end   <- 34e6

gdsfmt::showfile.gds(closeall = TRUE)

cat("[", format(Sys.time()), "] Step 1: LD pruning\n")
if (!file.exists(gds_file)) stop("SNPRelate GDS not found: ", gds_file)

genofile <- snpgdsOpen(gds_file)
all_chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
all_pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
all_snp <- read.gdsn(index.gdsn(genofile, "snp.id"))

if (drop_MHC) {
  chr_char <- as.character(all_chr)
  mhc_mask <- chr_char == mhc_chr & all_pos >= mhc_start & all_pos <= mhc_end
  snps_no_mhc <- all_snp[!mhc_mask]
  cat("  Variants after MHC drop:", length(snps_no_mhc), "\n")
} else {
  snps_no_mhc <- all_snp
}

set.seed(42)
pruned_snps <- snpgdsLDpruning(
  genofile,
  snp.id        = snps_no_mhc,
  autosome.only = TRUE,
  ld.threshold  = sqrt(0.1),
  slide.max.bp  = 500000L,
  verbose       = TRUE
)
pruned_ids <- unlist(pruned_snps, use.names = FALSE)
cat("  Variants after LD pruning:", length(pruned_ids), "\n")
snpgdsClose(genofile)

cat("[", format(Sys.time()), "] Step 2: Opening SeqArray GDS\n")
if (!file.exists(seq_gds_file)) stop("SeqArray GDS not found: ", seq_gds_file)

seqfile <- seqOpen(seq_gds_file, allow.duplicate = TRUE)
on.exit({ try(seqClose(seqfile), silent = TRUE) }, add = TRUE)

seqResetFilter(seqfile)
all_samples  <- seqGetData(seqfile, "sample.id")
all_seq_vars <- seqGetData(seqfile, "variant.id")

cat("  SeqArray samples:", length(all_samples), "\n")
cat("  SeqArray variants:", length(all_seq_vars), "\n")

seq_chr <- seqGetData(seqfile, "chromosome")
seq_pos <- seqGetData(seqfile, "position")

snp_pos_key <- paste(all_chr[match(pruned_ids, all_snp)],
                     all_pos[match(pruned_ids, all_snp)], sep = ":")
seq_pos_key <- paste(seq_chr, seq_pos, sep = ":")
keep_seq_vars <- all_seq_vars[seq_pos_key %in% snp_pos_key]

cat("  SeqArray variants after LD prune mapping:", length(keep_seq_vars), "\n")
if (length(keep_seq_vars) < 1000) {
  stop("Fewer than 1000 variants after pruning — check GDS position matching")
}
seqSetFilter(seqfile, variant.id = keep_seq_vars, verbose = FALSE)

cat("[", format(Sys.time()), "] Step 3: Loading PC-AiR PCs\n")
if (!file.exists(pcair_rds)) stop("PC-AiR RDS not found: ", pcair_rds)
pcair_obj <- readRDS(pcair_rds)

pcs_all <- pcair_obj$vectors
if (is.null(pcs_all)) pcs_all <- pcair_obj$eigenvectors
if (is.null(pcs_all)) stop("Cannot find eigenvectors in PC-AiR object")

pcs_mat <- pcs_all[, seq_len(min(N_PCS, ncol(pcs_all))), drop = FALSE]
cat("  Using", ncol(pcs_mat), "PCs from PC-AiR object\n")

active_samples <- seqGetData(seqfile, "sample.id")
match_idx <- match(active_samples, rownames(pcs_mat))
if (any(is.na(match_idx))) {
  stop("Some GDS sample IDs not found in PC-AiR eigenvectors")
}
pcs_mat <- pcs_mat[match_idx, , drop = FALSE]

unrels       <- pcair_obj$unrels
training_set <- intersect(unrels, active_samples)
cat("  Training set (unrelateds):", length(training_set), "\n")

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
  ibd.probs          = TRUE,
  scale              = "variant",
  small.samp.correct = TRUE,
  BPPARAM            = BPPARAM,
  verbose            = TRUE
)

cat("[", format(Sys.time()), "] Step 5: Saving results\n")
saveRDS(relate, file = out_rds)
cat("  Saved:", out_rds, "\n")

kin_pairs <- pcrelateToMatrix(relate, thresh = 2^(-9/2), scaleKin = 2)
saveRDS(kin_pairs, file = file.path(out_dir, paste0(NAME, "_pcrelate_kinmat.RDS")))

cat("[", format(Sys.time()), "] PC-Relate pipeline complete.\n")
