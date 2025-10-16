# Generate a delimited batch info filea
# Author: Christian Coffman
# Date Created: 25-10-09
# Last Updated: 25-10-09

library(tidyverse)

bf <- read_table("/projects/standard/basu_hbcd/shared/batch.info")
identifiers <- read_csv("../release_identifiers_20251006.csv") %>% select(release_candid, pscid) %>%
  mutate(release_candid = as.numeric(release_candid)) %>%
  drop_na()
fam <- read_table("HBCD_BR20.0.fam", col_names = c("FID", "IID", "PAT","MAT","SEX","PHENO"))

release <- fam %>%
  left_join(identifiers, by = c("FID" = "release_candid")) %>%
  mutate(IID.new = paste0(pscid, str_sub(IID, -1, -1))) %>%
  left_join(bf, by = c("IID.new" = "IID")) %>%
  select(FID, IID, plate_number, visit)

write_delim(release, "batch.info", delim = "\t")
