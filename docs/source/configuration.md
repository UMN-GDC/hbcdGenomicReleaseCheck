# Configuration

## Environment Variables

The pipeline is configured primarily through environment variables. Set these before running any pipeline steps.

### Required Variables

| Variable | Description | Default | Example |
|----------|-------------|---------|---------|
| `HBCD_RELEASE` | Release tag (used in output directory name) | `br_21p3` | `br31p2` |
| `HBCD_DATA_DIR` | Root directory for HBCD shared data | `/projects/standard/basu_hbcd/shared/data` | `/projects/standard/basu_hbcd/shared/data` |
| `HBCD_IMPUTATION_DIR` | Root for imputation inputs (VCFs, GP cache) | `/projects/standard/basu_hbcd/shared/hbcdSandboxData` | `/projects/standard/basu_hbcd/shared/hbcdSandboxData` |

### Optional Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `HBCD_SLURM_ACCOUNT` | SLURM account for job submission | None (uses default) |
| `HBCD_SLURM_PARTITION` | SLURM partition for jobs | None (uses default) |
| `HBCD_TOPMED_KEY` | Path to TOPMed API key | `~/topmedKey` |
| `CONDA_ENV` | Conda environment name | `gdcPipeline` |

### Input File Overrides (New)

| Variable | Description | Default |
|----------|-------------|---------|
| `HBCD_IDENTIFIERS_FILE` | Path to PSCID ↔ release_candid crosswalk CSV | `$HBCD_DATA_DIR/release_identifiers_20260628.csv` |
| `HBCD_PAR_VISIT_FILE` | Path to par_visit participation TSV | `$HBCD_DATA_DIR/par_visit_data_br21_1.tsv` |

These allow using different crosswalk/par_visit files per release without code changes.
Set via `run_release_pipeline.sh` or manually before running steps.

### Using Updated Crosswalk and par_visit Files

For each new release (e.g., `br31p2`), place updated files in `$HBCD_DATA_DIR`:

```bash
# Standard location
ls /projects/standard/basu_hbcd/shared/data/
# release_identifiers_20261201.csv   # Updated crosswalk
# par_visit_data_br31_2.tsv          # Updated par_visit
# HBCDexclusions.csv                 # Updated exclusions (if changed)
# HBCD_genetics_QC1_missing_race_LORIS.xlsx  # Updated QC Excel (if changed)
```

Then run with overrides (or let `run_release_pipeline.sh` export them):

```bash
export HBCD_RELEASE=br31p2
export HBCD_IDENTIFIERS_FILE=/projects/standard/basu_hbcd/shared/data/release_identifiers_20261201.csv
export HBCD_PAR_VISIT_FILE=/projects/standard/basu_hbcd/shared/data/par_visit_data_br31_2.tsv
./run_release_pipeline.sh
```

**Expected file formats:**

| File | Key Columns |
|------|-------------|
| `release_identifiers_YYYYMMDD.csv` | `pscid`, `candid`, `release_candid` |
| `par_visit_data_brXX_X.tsv` | `participant_id`, `par_visit_data_visit_missed` |

The pipeline reads these via `_lib.py` helpers which respect the env vars.

### Example Configuration

```bash
# ~/.bashrc or job script
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData
export HBCD_SLURM_ACCOUNT=hbcd_genomics
export HBCD_SLURM_PARTITION=msismall
```

## Release Directory Structure

The pipeline creates a release-specific directory:

```
HBCD_genomics_release_${HBCD_RELEASE}/
├── genotype_microarray/           # Main release outputs
│   ├── GDA/                       # Genotype data (PLINK)
│   ├── genesis/                   # GRM, PC-AiR, PC-Relate
│   ├── cnv/                       # CNV calls
│   ├── imputed/                   # Filtered imputed VCFs
│   ├── keep_list.txt              # Release subject whitelist
│   ├── temp.fam                   # De-ID fam (all QC subjects)
│   └── staging/                   # Intermediate files
│       └── cnv/
└── logs/                          # Step execution logs
```

## Input File Configuration

### Required Input Files (in `$HBCD_DATA_DIR`)

| File | Purpose | Format |
|------|---------|--------|
| `HBCD.{bed,bim,fam}` | Raw genotype data | PLINK binary |
| `onlyQc.{bed,bim,fam}` | Post-QC genotypes (Phase A output) | PLINK binary |
| `batch.info` | Batch/plate metadata | TSV (IID, visit, plate_number) |
| `release_identifiers_20260628.csv` | PSCID ↔ release_candid mapping | CSV |
| `HBCD_genetics_QC1_missing_race_LORIS.xlsx` | Excel exclusion list | XLSX (Exclude_Summary sheet) |
| `HBCDexclusions.csv` | Additional PSCID exclusions | CSV (reason columns) |
| `par_visit_data_br21_1.tsv` | par_visit participation | TSV |
| `data_handoff/hbcd_pcrelate_grm.grm.*` | PC-Relate GRM | Binary GRM |
| `data_handoff/hbcd_pcair_32PCs_clean.tsv` | PC-AiR weights | TSV |
| `data_handoff/hbcd_pcrelate_grm_pairwise.tsv` | PC-Relate pairwise | TSV |
| `data_handoff/CNV_slim_clean.txt` | CNV calls (pscid IIDs) | TSV |
| `data_handoff/HBCD_CNV_bookmark_metrics_clean.csv` | CNV bookmarks | CSV |

### Imputation Input Files (in `$HBCD_IMPUTATION_DIR/imputed`)

| Path | Description |
|------|-------------|
| `c1/chr{1-22}.dose.vcf.gz` | Imputation chunk 1 (autosomes) |
| `c2/chr{1-22}.dose.vcf.gz` | Imputation chunk 2 (autosomes) |
| `c3/chr{1-22}.dose.vcf.gz` | Imputation chunk 3 (autosomes) |
| `chrX.dose.vcf.gz` | Chr X imputation |
| `c1/chr{1-22}.info.gz` | Imputation INFO scores (chunk 1) |
| `c2/chr{1-22}.info.gz` | Imputation INFO scores (chunk 2) |
| `c3/chr{1-22}.info.gz` | Imputation INFO scores (chunk 3) |

## Script-Level Configuration

### Step 25: `25-filterGenotypeFiles.py`

Configuration via `_lib.py` constants:

```python
# In _lib.py
PAR_VISIT_FILE = "par_visit_data_br21_1.tsv"
IDENTIFIERS_FILE = "release_identifiers_20260628.csv"
EXCLUSION_FILE = "HBCD_genetics_QC1_missing_race_LORIS.xlsx"
ADDITIONAL_EXCLUSIONS_FILE = "HBCDexclusions.csv"
```

To use different files, set environment variables before running:

```bash
export PAR_VISIT_FILE=par_visit_data_br31_2.tsv
export IDENTIFIERS_FILE=release_identifiers_20261201.csv
python 25-filterGenotypeFiles.py
```

### Step 27: `27-filter_imputed_vcf.SLURM`

SLURM array configuration (edit the SLURM file):

```bash
# In 27-filter_imputed_vcf.SLURM
#SBATCH --array=0-23          # 24 tasks (chr 1-22 + X)
#SBATCH --mem=8G              # Memory per task
#SBATCH --time=02:00:00       # Time limit
#SBATCH --account=hbcd_genomics
#SBATCH --partition=msismall
```

### Step 26: `26-run_plink_filter.sh`

Reads `HBCD_RELEASE` directly:

```bash
export HBCD_RELEASE=br31p2
./26-run_plink_filter.sh
```

## SLURM Configuration

### Default SLURM Settings (in .SLURM files)

| Script | Array | Memory | Time | Partition |
|--------|-------|--------|------|-----------|
| `03-grm_gcta.SLURM` | 1-22 | 16G | 4:00:00 | msismall |
| `04-grm_plink.SLURM` | 1-22 | 8G | 2:00:00 | msismall |
| `05-pcair.SLURM` | - | 32G | 8:00:00 | msilarge |
| `10a-prepare_imputation.SLURM` | 0-22 | 8G | 2:00:00 | msismall |
| `12a-download_chunk1.SLURM` | 1-7 | 4G | 1:00:00 | msismall |
| `13a-unzip.SLURM` | 1-7 | 8G | 2:00:00 | msismall |
| `27-filter_imputed_vcf.SLURM` | 0-23 | 8G | 2:00:00 | msismall |

### Customizing SLURM Resources

Create a local override file:

```bash
# slurm_config.sh
export SLURM_ACCOUNT=${HBCD_SLURM_ACCOUNT:-hbcd_genomics}
export SLURM_PARTITION=${HBCD_SLURM_PARTITION:-msismall}
export SLURM_MEM_DEFAULT=8G
export SLURM_TIME_DEFAULT=2:00:00
```

Source before submitting:

```bash
source slurm_config.sh
sbatch 27-filter_imputed_vcf.SLURM
```

## Conda Environment

### `gdcPipeline` Environment (Shared)

```yaml
# Approximate environment specification
name: gdcPipeline
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - pandas>=2.0
  - numpy>=1.24
  - openpyxl>=3.1
  - pytest>=7.0
  - matplotlib>=3.7
  - scipy>=1.10
  - r-base=4.4
  - r-tidyverse
  - r-ggh4x
  - r-zoo
  - r-scales
  - plink2=2.00a091019
  - bcftools=1.18
  - gcta=1.94.1
  - king=2.2.5
```

### Local Development Environment

```bash
# From pyproject.toml
pip install -e ".[dev]"

# Or minimal
pip install pandas numpy openpyxl pytest matplotlib scipy
```

## TOPMed Configuration

### API Key Setup

```bash
# Get API key from TOPMed Imputation Server
# Save to ~/topmedKey (chmod 600)
echo "your_api_key_here" > ~/topmedKey
chmod 600 ~/topmedKey

# Or set custom location
export HBCD_TOPMED_KEY=/path/to/topmedKey
```

### Submission Script (`11-submit_topmed.sh`)

Configurable variables at top of script:

```bash
# In 11-submit_topmed.sh
TOPMED_KEY="${HBCD_TOPMED_KEY:-~/topmedKey}"
VCF_DIR="${HBCD_IMPUTATION_DIR}/imputed/prepared"
EMAIL="your.email@umn.edu"  # Notification email
```

## Version-Specific Configuration

### For Release `br31p2`

```bash
export HBCD_RELEASE=br31p2
export PAR_VISIT_FILE=par_visit_data_br31_2.tsv
export IDENTIFIERS_FILE=release_identifiers_20261201.csv
# ... other br31p2-specific files
```

### For Future Releases

Update the following files for each release:
- `PAR_VISIT_FILE` in `_lib.py`
- `IDENTIFIERS_FILE` in `_lib.py`
- `EXCLUSION_FILE` in `_lib.py`
- `par_visit_data_brXX_X.tsv` in data directory
- `release_identifiers_YYYYMMDD.csv` in data directory

## Validation

Verify configuration before running:

```bash
# Check environment
echo "Release: $HBCD_RELEASE"
echo "Data dir: $HBCD_DATA_DIR"
echo "Imputation dir: $HBCD_IMPUTATION_DIR"

# Check required files
for f in HBCD.fam onlyQc.fam batch.info release_identifiers_20260628.csv \
         HBCD_genetics_QC1_missing_race_LORIS.xlsx HBCDexclusions.csv \
         par_visit_data_br21_1.tsv; do
    test -f "$HBCD_DATA_DIR/$f" && echo "✓ $f" || echo "✗ $f MISSING"
done

# Check imputation files
for c in c1 c2 c3; do
    test -d "$HBCD_IMPUTATION_DIR/imputed/$c" && echo "✓ imputed/$c" || echo "✗ imputed/$c MISSING"
done
```