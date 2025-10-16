HBCD genomics QC validation

PLINK files
- .bed
- .bim
- .fam
  - Six columns
    - Family ID (FID)
    - Within-family ID (IID) i.e. FID plus relationship in family C = child, M = mother
    - Within-family ID of father
    - Within-family ID of mother
    - Sex
    - Phenotype
  - Space delimited
- batch.info
  - Tab delimited
  - 3 columns
    - IID (character, 10 digits + 1 letter either M or C)
    - plate_number (float)
    - visit (character)
  - .fam and batch.info have the same ordering of IID column
- Genomic Derivatives (Under Construction)
- .eigenvec
- .eigenval
- Inferred Sex
- Inferred Ancestry


# Prerelease staging file structure
- Data to be released will be in the “data” directory
- QC’d genomic data will be stored as Plink1 binary (.bed/.bim/.fam) format in the first level of the “data” directory
Derivatives associated with genomics release will be stored in “data/derivatives”





:w
gg:w
