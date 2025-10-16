HBCD genomics QC validation

PLINK1 files in "data"
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

- [ ] Genomic Derivatives (in "data/derivatives" Under Construction)
    - .eigenvec
    - .eigenval
    - Inferred Sex
    - Inferred Ancestry
