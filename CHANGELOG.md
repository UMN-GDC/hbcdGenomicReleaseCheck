# Changelog

All notable changes to the HBCD Genomic Release Pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial pipeline structure for HBCD genomic data release
- Phase A: QC & derivative computation (scripts 01-17)
- Phase B: De-identification (batch processing)
- Phase C: Release filtering & validation (scripts 25-30)
- Comprehensive test suite for release validation
- Automated exclusion validation
- Cross-form subject count consistency checks

### Changed
- N/A (initial development)

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A

## [0.1.0] - 2024-XX-XX

### Added
- **Step 01** (`01-onlyQCremoved.py`): QC removal - maps HBCD.fam IIDs to release_candid, produces QC_removed.txt
- **Step 02** (`02-plink_qc_remove.sh`): PLINK2 --remove to produce onlyQc.{bed,bim,fam}
- **Step 03** (`03-grm_gcta.SLURM`): GCTA GRM computation (per-chromosome parts, merged)
- **Step 04** (`04-grm_plink.SLURM`): PLINK2 GRM (--make-rel)
- **Step 05** (`05-pcair.SLURM`): PC-AiR via KING kinship estimation
- **Step 06** (`06-pca_ir_pipeline.R`): R script for PC-AiR ancestry-adjusted PCs
- **Step 07** (`07-pcrelate.sh`): PC-Relate pipeline wrapper
- **Step 08** (`08-pca_relate_pipeline.R`): R script for PC-Relate relatedness estimation
- **Step 09** (`09-export_pcrelate.R`): Export kinship/IBD to CSV/TSV
- **Step 10a** (`10a-prepare_imputation.SLURM`): Prepare autosomes 1-22 for TOPMed imputation
- **Step 10b** (`10b-prepare_imputation_X.SLURM`): Prepare chrX for TOPMed imputation
- **Step 11** (`11-submit_topmed.sh`): Submit VCFs to TOPMed imputation server
- **Step 12a-c** (`12a-download_chunk1.SLURM`, `12b-download_chunk2.SLURM`, `12c-download_chunk3.SLURM`): Download imputation results
- **Step 13a-d** (`13a-unzip.SLURM` through `13d-unzip.SLURM`): Unzip imputation chunks
- **Step 14** (`futureScripts/14a-precompute_gp.SLURM`): GP precomputation per chromosome × ancestry
- **Step 15** (`futureScripts/15-aggregate_gp.SLURM`): Aggregate GP results
- **Step 16** (`futureScripts/16-imputation-quality-report.qmd`): Imputation QC report
- **Step 17** (`futureScripts/17-cnv-qc-report.qmd`): CNV genomic profile report
- **Step 25** (`25-filterGenotypeFiles.py`): De-identification + filtering of genotype data
- **Step 26** (`26-run_plink_filter.sh`): PLINK2 --keep to produce release PLINK files
- **Step 27** (`27-filter_imputed_vcf.SLURM`): Filter imputed VCFs to release subjects
- **Step 28** (`28-cnv-deid.py`): CNV de-identification (pscid → release_candid)
- **Step 29** (`29-filter_release_outputs.py`): Filter all derivatives to release subjects
- **Step 30** (`30-validateExclusions.py`): Validate no excluded subjects in any output
- **Validation tests** (`tests/test_release_data.py`, `tests/test_exclusions.py`): 35+ tests covering de-ID integrity, filter correctness, derivative integrity, imputed VCF validation, exclusion compliance
- **Wrapper script** (`run_release_pipeline.sh`): Complete Phase C pipeline orchestration
- **Documentation**: Comprehensive README.md with pipeline overview and usage instructions
- **Filter concept document** (`00-filterConcept.qmd`): Venn diagrams and pipeline flowchart

### Infrastructure
- Apache 2.0 license
- pyproject.toml with dependencies and metadata
- Semantic versioning (hbcd_genomic_release/_version.py)
- pytest configuration with markers
- Ruff linting configuration
- Type checking with mypy
- GitHub Actions ready structure

### Documentation
- README.md: Full pipeline documentation with phases, steps, and commands
- imputation-README.md: Imputation-specific pipeline documentation
- Inline script documentation

## Release Targets

| Pipeline Version | HBCD Release | Status |
|------------------|--------------|--------|
| 0.1.0 | br31p2 | Target |

## Migration Guide

### From Pre-Release to 0.1.0
This is the initial release. No migration needed.

## Support Policy

- **Current release**: Full support, bug fixes, security patches
- **Previous release (n-1)**: Security patches only
- **Older releases**: No official support

See [README.md](README.md) for installation and usage instructions.