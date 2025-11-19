# title: Fam file release check
# author: Christian Coffman
# Author: Christian Coffman
# Date Created: 25-10-09
# Last Updated: 25-11-19
# [Further information](https://docs.google.com/document/d/1dYBvSBbSGM_-6453IT4EZkIUkvk6vaVQcLL-FPADz5Y/edit?tab=t.0)

library(tidyverse)
library(testthat)
dataPrefix="/scratch.global/hbcd/full/full.QC8"
release="br_20p1"
releaseDir=paste0("../../../HBCD_genomics_release_", release, "/data/")

fam <- read_table(paste0(releaseDir, "hbcd.fam"), col_names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
batch <- read_table(paste0(releaseDir, "batch.info"))
excluded <- read_delim(paste0(releaseDir, "../excluded.txt"), delim = "\t")

combined <- full_join(fam, batch, by = "IID") |>
  full_join(excluded, by = "IID")
idLengths <- combined |>
  mutate(FID = nchar(FID), IID = nchar(IID)) |>
  count(FID, IID)
iidTest <- full_join(fam, batch, by = "IID") |>
  full_join(excluded, by = "IID") |>
  mutate(
    IID2 = as.numeric(str_sub(IID, 1, 10)),
    Relation = str_sub(IID, 11, 11),
  ) |>
  filter(! IID %in% excluded$IID)
# missing FID are in the excluded file


test_that("Fam and batch order match", {
  expect_equal(mean(fam$IID == batch$IID), 1)
})
test_that("DeID check: ID lengths are expected", {
  expect_equal(mean(idLengths$FID %in% c(NA, 10)), 1)
  expect_equal(mean(idLengths$IID %in% c(11)), 1)
  expect_equal(mean(is.numeric(combined$FID)), 1)
})
test_that("DeID check: IID is composed of FID and either C or M", {
  expect_equal(mean(iidTest$IID2 == iidTest$FID), 1)
  expect_equal(mean((iidTest$Relation == "C") | (iidTest$Relation == "M")), 1)
})
