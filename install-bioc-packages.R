#!/usr/bin/env Rscript
# Install R packages for the HBCD release pipeline.
#
# Two usage modes:
#   1) Conda-based:  conda activate hbcd-release && Rscript install-bioc-packages.R
#   2) HPC-module:   module load R/4.4.2-openblas-rocky8 && Rscript install-bioc-packages.R
#
# Packages are installed to the personal/site R library.
# BiocManager auto-selects the correct Bioconductor version for the R version.

# module load R/4.4.2-openblas-rocky8

pkgs <- c(
    "Matrix",
    "SNPRelate",
    "SeqArray",
    "SeqVarTools",
    "GENESIS",
    "BiocParallel"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.r-project.org")
}

BiocManager::install(pkgs, ask = FALSE, update = TRUE, Ncpus = 4)

for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' failed to install", pkg))
    }
}
cat(sprintf("All %d packages installed successfully.\n", length(pkgs)))
