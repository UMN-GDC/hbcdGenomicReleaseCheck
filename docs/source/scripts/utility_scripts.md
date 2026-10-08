# Utility Scripts

## run_release_pipeline.sh

**Purpose**: Complete Phase C wrapper — runs steps 25→30 in order with SLURM monitoring

**Usage**:
```bash
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData
./run_release_pipeline.sh
```

**Environment Variables** (with defaults):
```bash
HBCD_RELEASE=br_21p3
HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData
```

**Flow**:
```bash
#!/bin/bash
set -euo pipefail

# 1. Activate environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline

# 2. Step 25: De-ID + filter genotypes
python 25-filterGenotypeFiles.py

# 3. Step 26: PLINK2 --keep
./26-run_plink_filter.sh

# 4. Step 27: Filter imputed VCFs (SLURM array + polling)
sbatch 27-filter_imputed_vcf.SLURM
wait_for_slurm_array "filter_imputed"

# 5. Step 28: CNV de-identification
python 28-cnv-deid.py

# 6. Step 29: Filter all derivatives
python 29-filter_release_outputs.py

# 7. Step 30: Validate exclusions
python 30-validateExclusions.py
```

**SLURM Wait Function**:
```bash
wait_for_slurm_array() {
    local job_name=$1
    while squeue -u $USER --name=$job_name --noheader | grep -q .; do
        sleep 30
    done
}
```

**Output**: Logs each step, final validation summary

---

## _lib.py

**Purpose**: Shared library for paths, data loaders, and common utilities

**Location**: Repository root (importable as `from _lib import ...`)

### Key Functions

#### Path Helpers
```python
DATA_DIR = Path(os.environ.get("HBCD_DATA_DIR", "/projects/standard/basu_hbcd/shared/data"))
DEFAULT_RELEASE = "br_21p3"
_RELEASE_ENV = os.environ.get("HBCD_RELEASE", DEFAULT_RELEASE)

def get_release_dir(release=None):
    """Return genotype_microarray directory for release."""
    if release is None:
        release = _RELEASE_ENV
    return DATA_DIR.parent / f"HBCD_genomics_release_{release}" / "genotype_microarray"

def get_release_base(release=None):
    """Return release base directory (parent of genotype_microarray)."""
    return get_release_dir(release).parent

def get_release_staging_dir(release=None):
    """Return staging directory for intermediate files."""
    return get_release_base(release) / "staging"
```

#### Data Loaders
```python
PAR_VISIT_FILE = "par_visit_data_br21_1.tsv"
IDENTIFIERS_FILE = "release_identifiers_20260628.csv"
EXCLUSION_FILE = "HBCD_genetics_QC1_missing_race_LORIS.xlsx"

def load_par_visit_candids(path=None):
    """Return set of candid from par_visit (completed visits only)."""
    if path is None:
        path = DATA_DIR / PAR_VISIT_FILE
    pv = pd.read_csv(path, sep="\t").query("par_visit_data_visit_missed == 'No'")
    def _to_int(v):
        s = str(v)
        return int(s[4:] if s.startswith("sub-") else s)
    return set(int(x) for x in pv["participant_id"].dropna().apply(_to_int).unique())

def load_identifiers(path=None):
    """Load release-identifiers CSV, parse numeric columns."""
    if path is None:
        path = DATA_DIR / IDENTIFIERS_FILE
    return (pd.read_csv(path)
            .query("release_candid != 'release_candid'")
            .assign(release_candid=lambda x: pd.to_numeric(x["release_candid"]))
            .assign(candid=lambda x: pd.to_numeric(x["candid"]))
            .dropna(subset=["release_candid"])
            .drop_duplicates(subset=["pscid", "candid", "release_candid"]))

def load_excluded_release_candids(path=None):
    """Return set of release_candid from Excel Exclude_Summary sheet."""
    if path is None:
        path = DATA_DIR / EXCLUSION_FILE
    exc = pd.read_excel(path, sheet_name="Exclude_Summary")
    exc[["_", "release_candid", "_pscid"]] = exc["Sampled ID"].str.split("_", n=2, expand=True)
    exc["release_candid"] = pd.to_numeric(exc["release_candid"])
    return set(int(x) for x in exc.dropna(subset=["release_candid"])["release_candid"].unique())

def load_excluded_with_relationship(path=None):
    """Return DataFrame with release_candid + relationship from Excel."""
    if path is None:
        path = DATA_DIR / EXCLUSION_FILE
    exc = pd.read_excel(path, sheet_name="Exclude_Summary")
    exc[["_", "release_candid", "_pscid"]] = exc["Sampled ID"].str.split("_", n=2, expand=True)
    exc["release_candid"] = pd.to_numeric(exc["release_candid"])
    exc["relationship"] = exc["Study_ID"].str[-1]
    return exc.dropna(subset=["release_candid"])[["release_candid", "relationship"]]

ADDITIONAL_EXCLUSIONS_FILE = "HBCDexclusions.csv"

def load_additional_excluded_pscids():
    """Load HBCDexclusions.csv → set of all pscids to exclude."""
    raw = pd.read_csv(DATA_DIR / ADDITIONAL_EXCLUSIONS_FILE)
    excluded = set()
    for col in raw.columns:
        vals = raw[col].dropna().astype(str).str.strip()
        excluded.update(vals[vals != "nan"])
    return excluded
```

### Usage in Scripts
```python
from _lib import (
    DATA_DIR,
    get_release_dir,
    get_release_base,
    load_par_visit_candids,
    load_identifiers,
    load_excluded_release_candids,
    load_excluded_with_relationship,
    load_additional_excluded_pscids,
)
```

---

## compute_gp_stats.py

**Purpose**: Compute genotype probability statistics from imputed VCFs

**Usage**:
```bash
python compute_gp_stats.py [options]
```

**Note**: Currently in development / future use

---

## test_ancestry.R

**Purpose**: Ancestry analysis and visualization

**Usage**:
```bash
Rscript test_ancestry.R
```

**Dependencies**: `tidyverse`, `ggplot2`, `GENESIS`

**Outputs**: Ancestry plots, population structure analysis

---

## Pipeline Flow Summary

```
_lib.py
    ├── Path management (get_release_dir, get_release_base, get_release_staging_dir)
    ├── Data loading (par_visit, identifiers, exclusions)
    └── Used by: 25, 28, 29, 30, tests, run_release_pipeline.sh

run_release_pipeline.sh
    ├── Orchestrates: 25 → 26 → 27(SLURM) → 28 → 29 → 30
    ├── Handles SLURM polling
    └── Uses env vars for configuration

25-filterGenotypeFiles.py
    ├── Uses _lib for all data loading
    ├── Produces: temp.fam, keep_list.txt, batch.info, removed_individuals.txt
    └── Core filtering logic

26-run_plink_filter.sh
    ├── Reads HBCD_RELEASE env var
    ├── PLINK2 --keep with temp.fam
    └── Produces: merged_chroms.{bed,bim,fam}

27-filter_imputed_vcf.SLURM
    ├── SLURM array 0-23 (chr 1-22 + X)
    ├── Reads keep_list.txt
    ├── bcftools filter VCFs
    └── Produces: imputed/chr*_dose.vcf.gz + QC stats

28-cnv-deid.py
    ├── De-identifies CNV_slim_clean.txt (pscid → release_candid)
    ├── Copies CNV bookmarks (already de-IDed)
    └── Produces: staging/cnv/CNV_slim_deid.txt, CNV_bookmarks_deid.csv

29-filter_release_outputs.py
    ├── Reads keep_list.txt → release_iids
    ├── Filters all handoff derivatives + CNV
    └── Produces: genesis/ + cnv/ in release dir

30-validateExclusions.py
    ├── Loads all exclusion sources
    ├── Checks every output file for excluded IIDs
    ├── Cross-form consistency check
    └── Logs detailed validation report
```