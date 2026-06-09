#!/usr/bin/env Rscript
# Install Bioconductor packages for the HBCD release pipeline.
# Run after: conda env create -f environment.yaml && conda activate hbcd-release
# Usage:   Rscript install-bioc-packages.R

pkgs <- c(
    "SNPRelate",
    "SeqArray",
    "SeqVarTools",
    "GENESIS",
    "BiocParallel"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cran.r-project.org")
}

BiocManager::install(pkgs, ask = FALSE, update = TRUE)

for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' failed to install", pkg))
    }
}
cat(sprintf("All %d packages installed successfully.\n", length(pkgs)))
