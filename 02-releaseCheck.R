# title: Fam file release check
# author: Christian Coffman
# Author: Christian Coffman
# Date Created: 25-10-09
# Last Updated: 25-11-19
# [Further information](https://docs.google.com/document/d/1dYBvSBbSGM_-6453IT4EZkIUkvk6vaVQcLL-FPADz5Y/edit?tab=t.0)

library(tidyverse)
library(testthat)
library(kableExtra)
dataPrefix="/scratch.global/hbcd/full/full.QC8"
release="br_20p1"
releaseDir=paste0("../HBCD_genomics_release_", release, "/data/")

fam <- read_table(paste0(releaseDir, "hbcd.fam"), col_names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
batch <- read_table(paste0(releaseDir, "batch.info"))
excluded <- read_delim(paste0(releaseDir, "../excluded.txt"), delim = "\t")

# fam and batch order
mean(fam$IID == batch$IID)

# deIDed checks
combined <- full_join(fam, batch, by = "IID") |>
  full_join(excluded, by = "IID")

combined |>
  mutate(FID = nchar(FID), IID = nchar(IID)) |>
  count(FID, IID)
# only 10s and 11s means they were deID

mean(is.numeric(combined$FID))


iidTest <- full_join(fam, batch, by = "IID") |>
  full_join(excluded, by = "IID") |>
  mutate(
    IID2 = as.numeric(str_sub(IID, 1, 10)),
    Relation = str_sub(IID, 11, 11),
  ) |>
  filter(! IID %in% excluded$IID)

# missing FID should be in the excluded file

# FID matches IID
mean(iidTest$IID2 == iidTest$FID)

# test relationships appended to FID to create IID
mean((iidTest$Relation == "C") | (iidTest$Relation == "M"))
