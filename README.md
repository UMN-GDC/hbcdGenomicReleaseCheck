# HBCD Genomic Release Pipeline

**De-identification, filtering, and validation of genomic data for HBCD Study releases.**

[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://hbcd-genomic-release.readthedocs.io/)
[![License](https://img.shields.io/badge/license-Apache%202.0-green)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)

---

## Quick Reference — Phase C (Release Filtering & Validation)

> **Prerequisites**: Phases A & B complete (see [docs](https://hbcd-genomic-release.readthedocs.io/en/latest/running.html#phase-a-qc-derivative-computation))

### Environment Setup

The pipeline runs as a **SLURM batch job** (12hr, 8GB). The wrapper script handles environment activation internally.

**Edit `run_release_pipeline.sh` first** — set your release config at the top:

```bash
# User-configurable variables (EDIT THESE before submitting)
HBCD_RELEASE="br31p2"
HBCD_DATA_DIR="/projects/standard/basu_hbcd/shared/data"
HBCD_IMPUTATION_DIR="/projects/standard/basu_hbcd/shared/hbcdSandboxData"

# Input file paths (update filenames per release)
HBCD_IDENTIFIERS_FILE="${HBCD_DATA_DIR}/release_identifiers_20261201.csv"
HBCD_PAR_VISIT_FILE="${HBCD_DATA_DIR}/par_visit_data_br31_2.tsv"
```

### Automated (Recommended)

```bash
# Submit as SLURM batch job (runs steps 25→30, including step 27 array)
sbatch run_release_pipeline.sh
```

The script:
1. Initializes environment (modules, .venv, conda) for SLURM batch
2. Runs step 25 (de-ID + filter genotypes)
3. Runs step 26 (PLINK2 --keep)
4. Submits step 27 (SLURM array, 24 tasks) and polls until complete
5. Runs step 28 (CNV de-identification)
6. Runs step 29 (filter all derivatives)
7. Runs step 30 (validate exclusions)

Output/error logs: `/projects/standard/basu_hbcd/shared/HBCD_genomics_release_hbcd_release_pipeline_<jobid>.out/.err`

### Manual Step-by-Step
```bash
# 1. Build release whitelist (de-ID + filter genotypes)
python 25-filterGenotypeFiles.py
# → keep_list.txt, temp.fam, GDA/batch.info, GDA/removed_individuals.txt

# 2. PLINK2 --keep → release PLINK files
./26-run_plink_filter.sh
# → GDA/merged_chroms.{bed,bim,fam}

# 3. Filter imputed VCFs (SLURM array, 24 tasks)
sbatch 27-filter_imputed_vcf.SLURM
# Wait: squeue -u $USER --name=filter_imputed

# 4. CNV de-identification (pscid → release_candid)
python 28-cnv-deid.py
# → staging/cnv/CNV_slim_deid.txt, CNV_bookmarks_deid.csv

# 5. Filter all derivatives to release subjects
python 29-filter_release_outputs.py
# → genesis/ + cnv/

# 6. Validate no excluded subjects in outputs
python 30-validateExclusions.py

# 7. Run test suite
python -m pytest tests/ -v
```

### Re-run Partial Pipeline (Individual Steps)

To run only a specific step or subset, set the required environment variables and run the script directly. All scripts read from `HBCD_RELEASE`, `HBCD_DATA_DIR`, `HBCD_IMPUTATION_DIR`, `HBCD_IDENTIFIERS_FILE`, `HBCD_PAR_VISIT_FILE`.

**Example: Run just the QC stats copy (after step 27 VCFs done):**
```bash
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData
./copy_imputation_qc.sh
```

**Example: Run just step 28 (CNV de-ID):**
```bash
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
source /projects/standard/basu_hbcd/shared/.venv/bin/activate
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/standard/gdc/public/envs/gdcPipeline
python 28-cnv-deid.py
```

**Example: Run just step 29 (filter derivatives):**
```bash
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData
source /projects/standard/basu_hbcd/shared/.venv/bin/activate
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/standard/gdc/public/envs/gdcPipeline
python 29-filter_release_outputs.py
```

**Example: Run just step 30 (validate exclusions):**
```bash
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
source /projects/standard/basu_hbcd/shared/.venv/bin/activate
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/standard/gdc/public/envs/gdcPipeline
python 30-validateExclusions.py
```

**Key environment variables:**
| Variable | Purpose | Default |
|----------|---------|---------|
| `HBCD_RELEASE` | Release tag (e.g., br31p2) | br_21p3 |
| `HBCD_DATA_DIR` | Base data directory | /projects/standard/basu_hbcd/shared/data |
| `HBCD_IMPUTATION_DIR` | Imputation data root | /projects/standard/basu_hbcd/shared/hbcdSandboxData |
| `HBCD_IDENTIFIERS_FILE` | Crosswalk CSV path | $HBCD_DATA_DIR/release_identifiers_20260628.csv |
| `HBCD_PAR_VISIT_FILE` | par_visit TSV path | $HBCD_DATA_DIR/par_visit_data_br21_1.tsv |

**Required environment activation for Python steps:**
```bash
# .venv for Python packages
source /projects/standard/basu_hbcd/shared/.venv/bin/activate
# Conda for PLINK/bcftools/R
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/standard/gdc/public/envs/gdcPipeline
```

---

### Re-run After Exclusion Changes
```bash
python 25-filterGenotypeFiles.py
./26-run_plink_filter.sh
sbatch 27-filter_imputed_vcf.SLURM
python 28-cnv-deid.py
python 29-filter_release_outputs.py
python 30-validateExclusions.py
python -m pytest tests/ -v
```

---

## Documentation

**Full documentation on ReadTheDocs:** https://hbcd-genomic-release.readthedocs.io/

| Topic | Link |
|-------|------|
| Installation & requirements | [Installation](https://hbcd-genomic-release.readthedocs.io/en/latest/installation.html) |
| Complete pipeline guide | [Running the Pipeline](https://hbcd-genomic-release.readthedocs.io/en/latest/running.html) |
| Configuration reference | [Configuration](https://hbcd-genomic-release.readthedocs.io/en/latest/configuration.html) |
| Output file formats | [Output Files](https://hbcd-genomic-release.readthedocs.io/en/latest/outputs.html) |
| Exclusion validation logic | [Exclusion Validation](https://hbcd-genomic-release.readthedocs.io/en/latest/exclusion_validation.html) |
| Testing guide | [Testing](https://hbcd-genomic-release.readthedocs.io/en/latest/testing.html) |
| Troubleshooting | [Troubleshooting](https://hbcd-genomic-release.readthedocs.io/en/latest/troubleshooting.html) |
| Script reference | [Script Reference](https://hbcd-genomic-release.readthedocs.io/en/latest/scripts/index.html) |
| Contributing | [Contributing](https://hbcd-genomic-release.readthedocs.io/en/latest/contributing.html) |

---

## Release Output Structure

```
HBCD_genomics_release_${HBCD_RELEASE}/
├── genotype_microarray/
│   ├── GDA/
│   │   ├── merged_chroms.{bed,bim,fam}   # Release genotypes
│   │   ├── batch.info                     # Visit + plate metadata
│   │   └── removed_individuals.txt        # Excluded IIDs
│   ├── genesis/
│   │   ├── pcair_weights.tsv              # PC-AiR (32 PCs)
│   │   ├── pcrelate_relatedness.grm.*     # GRM (binary + text)
│   │   └── pcrelate_relatedness.tsv       # Pairwise relatedness
│   ├── cnv/
│   │   ├── CNV_slim.txt                   # CNV calls
│   │   └── CNV_bookmarks.csv              # CNV bookmarks
│   ├── imputed/
│   │   ├── chr{1-22}_dose.vcf.gz + .tbi   # Filtered imputed VCFs
│   │   ├── chrX_dose.vcf.gz + .tbi
│   │   └── batch*-*.txt/.html             # QC statistics
│   ├── keep_list.txt                      # Release whitelist
│   ├── temp.fam                           # De-ID fam (all QC subjects)
│   └── staging/cnv/                       # Intermediates
```

---

## Key Principles

| Principle | Description |
|-----------|-------------|
| **De-ID first** | Map `pscid` → `release_candid` before any filtering |
| **Filter second** | Apply `keep_list.txt` to all derivatives |
| **Validate always** | Step 30 + test suite are mandatory quality gates |
| **Re-run Phase C only** | Exclusion changes don't require Phase A/B re-run |

---

## Support

- **Issues**: [GitHub Issues](https://github.com/hbcd-genomics/hbcd-genomic-release/issues)
- **Documentation**: [ReadTheDocs](https://hbcd-genomic-release.readthedocs.io/)
- **HBCD Genomics Team**: hbcd-genomics@umn.edu

---

## License

Apache License 2.0 — see [LICENSE](LICENSE) for details.
