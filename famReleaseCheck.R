library(tidyverse)
library(testthat)
library(kableExtra)

info <- read_table("/projects/standard/basu_hbcd/shared/HBCD_genomics_release_excludecontrols/batch.info")
info <- read_table("/projects/standard/basu_hbcd/shared/HBCD_genomics_release_br_20p0/concat/genetics/genotype_microarray/GDA/batch.info")
info |>
  mutate(
    IID = nchar(IID), 
  ) |>
  count(IID)

info |> count(visit)
info |> filter(is.na(visit))

release <- read_table("/projects/standard/basu_hbcd/shared/HBCD_genomics_release_br_20p0/concat/genetics/genotype_microarray/GDA/HBCD.fam", col_names=c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))

release |>
  mutate(
    across(c(FID, IID), nchar), 
  ) |>
  count(FID, IID)


convTable <- release %>%
  left_join(identifiers, by = c("FID" = "release_candid")) %>%
  mutate(
    pscidR = paste0(pscid, str_sub(IID, -1, -1)),
  ) %>%
  full_join(original, by = "pscidR", suffix= c(".release", ".original")) %>%
  mutate(
    relation.release = str_sub(IID.release, -1, -1),
    relation.original = str_sub(IID.original, -1, -1)
  )


# make sure everyone is deidentified
# temp <- release %>%
#   filter(!grepl("_Control_", FID)) %>%
#   filter(!grepl("^NTC", FID)) %>%
#   mutate(
#     fidLength = nchar(FID),
#     iidLength = nchar(IID),
#   ) %>%
#   arrange(desc(fidLength)) %>%
#   head(1)

colSums(is.na(convTable)) %>%
  kable(caption = "Missingness count after joining all of the data sources.")
```
- The missing values are from the individual who wasn't deidentified
  - additionally, this individual has an IID starting with "QI" (note that everyone else has IID starting with "CH")


# Comparisons 
```{r}
# Same 
print("Proportion same Sex")
mean(convTable$SEX.original == convTable$SEX.release, na.rm = T) 
print("Proportion same relation")
mean(convTable$relation.original == convTable$relation.release, na.rm = T)
print("Proportion same sex/relation")
mean(interaction(convTable$relation.original, convTable$SEX.original) == interaction(convTable$relation.release, convTable$SEX.release), na.rm = T)


print("Proportion same ID")
mean(convTable$pscid.original == convTable$pscid.release, na.rm = T)
# convTable %>%
#   count(pscid.original, pscid.release) %>%
#   group_by(pscid.original, pscid.release) %>%
#   filter(n() > 1)
# empty table shows that FID is unique mapping between two data tables
```
  - [x] Counts
    - [x] Sex
    - [x] Mother/Child
    - [x] Mother/Child * Sex
    - [x] FID
    - All these counts will contain missingness counts if there are any


  - [x] Order
```{r}
print("Proportion same comparing new vs old order")
mean(convTable$pscid.release == original$pscid)

#temp <- convTable[which(convTable$pscid.release != original$pscid), ]

```

- Everything looks good once we solve the issues associated with the not-yet-deidentified individual

