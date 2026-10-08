# HBCD Genomic Release Pipeline

**De-identification, filtering, and validation of genomic data for HBCD Study releases.**

[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue)](https://hbcd-genomic-release.readthedocs.io/)
[![License](https://img.shields.io/badge/license-Apache%202.0-green)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)

---

## Quick Reference — Phase C (Release Filtering & Validation)

> **Prerequisites**: Phases A & B complete (see [docs](https://hbcd-genomic-release.readthedocs.io/en/latest/running.html#phase-a-qc-derivative-computation))

### Environment Setup
best to run this in an srun
```bash
srun --mem=8gb --time=8:00:00 --pty bash
```


```bash
# MSI modules
module load plink/2.00-alpha-091019
module load bcftools
module load R/4.4.2-openblas-rocky8

# Python virtual env (for pipeline Python scripts)
source /projects/standard/basu_hbcd/shared/.venv/bin/activate

# Conda env — activate with FULL PATH before running pipeline
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/standard/gdc/public/envs/gdcPipeline

# Release config (adjust for your release)
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData

# Input file overrides (place updated files in HBCD_DATA_DIR)
export HBCD_IDENTIFIERS_FILE=$HBCD_DATA_DIR/release_identifiers_20261201.csv
export HBCD_PAR_VISIT_FILE=$HBCD_DATA_DIR/par_visit_data_br31_2.tsv
```

### Automated (Recommended)
```bash
./run_release_pipeline.sh
```
Runs steps 25→30 with SLURM monitoring.

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
