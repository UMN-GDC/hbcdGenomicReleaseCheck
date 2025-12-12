HBCD genomics QC validation

# Running DEID and validate
- Update release version in 01-setupNfilter.qmd for the shell portion and the R portion
- Updated release version in tests/testthat/test-ids.R
- Run 01-setupNfilter.qmd
- Run 02-validate.R 
- Checdk to see if test-log.txt is clean



# Data for release
- .bed
- .bim
- .fam
  - [ ] Six columns
    - Family ID (FID)
    - [ ] Within-family ID (IID) i.e. FID plus relationship in family C = child, M = mother
    - Within-family ID of father
    - Within-family ID of mother
    - Sex
    - Phenotype
  - [ ] Space delimited
- batch.info
  - [ ] Tab delimited
  - [ ] 3 columns
    - [ ] IID (character, 10 digits + 1 letter either M or C)
    - [ ] plate_number (float)
    - [ ] visit (character)
  - [ ] .fam and batch.info have the same ordering of IID column
