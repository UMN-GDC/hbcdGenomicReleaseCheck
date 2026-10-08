# Installation

## Requirements

### System Requirements

- **OS**: Linux (tested on Rocky Linux 8/9, Ubuntu 20.04+)
- **Python**: ≥ 3.10
- **R**: ≥ 4.3 (with tidyverse, ggh4x, zoo, scales)
- **PLINK**: 2.00-alpha (2019-09-10 or later)
- **bcftools**: ≥ 1.15
- **GCTA**: ≥ 1.94.1 (for GRM computation)
- **KING**: ≥ 2.2.5 (for PC-AiR kinship)
- **TOPMed API access**: Valid TOPMed imputation server account

### HPC Environment (MSI/UMN)

The pipeline is designed to run on the Minnesota Supercomputing Institute (MSI) clusters with SLURM job scheduler.

```bash
# Required modules on MSI
module load plink/2.00-alpha-091019
module load bcftools
module load R/4.4.2-openblas-rocky8
module load gcta/1.94.1
module load king/2.2.5
```

## Conda Environment

### Option 1: Shared HBCD Environment (Recommended)

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline
```

The `gdcPipeline` environment includes:
- Python 3.10+ with pandas, numpy, openpyxl, pytest, matplotlib
- R with tidyverse, ggh4x, zoo, scales
- All required system tools via module system

### Option 2: Local Conda Environment

```bash
# Create from pyproject.toml
conda create -n hbcd-genomic-release python=3.10
conda activate hbcd-genomic-release
pip install -e ".[dev]"

# Or install core dependencies only
pip install pandas numpy openpyxl pytest matplotlib scipy
```

### Option 3: Docker/Apptainer (Planned)

```bash
# Future: Containerized execution
apptainer pull docker://ghcr.io/hbcd-genomics/hbcd-genomic-release:0.1.0
```

## Directory Structure

```
hbcd-genomic-release/
├── *.py, *.sh, *.SLURM, *.R, *.qmd    # Pipeline scripts
├── tests/                              # Validation test suite
├── futureScripts/                      # Future/planned scripts
├── setup/                              # Installation scripts
├── docs/                               # Documentation (this site)
├── hbcd_genomic_release/               # Package with version info
├── pyproject.toml                      # Package metadata & dependencies
├── README.md                           # Main documentation
├── CHANGELOG.md                        # Version history
└── LICENSE                             # Apache 2.0 license
```

## Data Prerequisites

The pipeline expects the following data structure in `$HBCD_DATA_DIR`:

```
/projects/standard/basu_hbcd/shared/data/
├── HBCD.{bed,bim,fam}                  # Raw genotype data
├── HBCD.fam                            # Original PLINK .fam with pscid IIDs
├── batch.info                          # Batch info with IID={pscid}{C|M}
├── onlyQc.{bed,bim,fam}                # Post-QC genotype data (Phase A output)
├── release_identifiers_20260628.csv    # PSCID → release_candid crosswalk
├── HBCD_genetics_QC1_missing_race_LORIS.xlsx  # Excel exclusion list
├── HBCDexclusions.csv                  # Additional PSCID exclusions
├── par_visit_data_br21_1.tsv           # par_visit participation data
└── data_handoff/                       # External derivative files
    ├── hbcd_pcrelate_grm.grm.*         # PC-Relate GRM files
    ├── hbcd_pcair_32PCs_clean.tsv      # PC-AiR weights
    ├── hbcd_pcrelate_grm_pairwise.tsv  # PC-Relate pairwise
    ├── CNV_slim_clean.txt              # CNV calls (pscid IIDs)
    └── HBCD_CNV_bookmark_metrics_clean.csv  # CNV bookmarks
```

And imputation inputs in `$HBCD_IMPUTATION_DIR`:

```
/projects/standard/basu_hbcd/shared/hbcdSandboxData/
├── imputed/
│   ├── c1/chr*.dose.vcf.gz             # Imputation chunk 1
│   ├── c2/chr*.dose.vcf.gz             # Imputation chunk 2
│   ├── c3/chr*.dose.vcf.gz             # Imputation chunk 3
│   └── chrX.dose.vcf.gz                # Chr X imputation
```

## Verification

Verify installation:

```bash
# Check Python environment
python -c "import pandas, numpy, openpyxl, pytest; print('OK')"

# Check system tools
plink2 --version
bcftools --version
R --version

# Run test suite (requires data)
pytest tests/ -v
```

## Troubleshooting Installation

| Issue | Solution |
|-------|----------|
| `ModuleNotFoundError: pandas` | Activate conda env: `conda activate gdcPipeline` |
| `plink2: command not found` | Load module: `module load plink/2.00-alpha-091019` |
| `bcftools: command not found` | Load module: `module load bcftools` |
| `R: command not found` | Load module: `module load R/4.4.2-openblas-rocky8` |
| Permission denied on SLURM | Check SLURM account: `sbatch --account=your_account` |
| TOPMed upload fails | Verify `~/topmedKey` exists and is valid |