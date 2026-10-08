# Testing Guide

## Overview

The pipeline includes a comprehensive test suite validating de-identification integrity, filter correctness, derivative integrity, imputed VCF validation, and exclusion compliance.

## Test Structure

```
tests/
├── conftest.py              # Pytest configuration & fixtures
├── test_release_data.py     # 31 tests: de-ID, data integrity, derivatives, VCFs
├── test_exclusions.py       # 4 tests: exclusion compliance
└── fixtures/                # Test data (future)
```

## Running Tests

### Prerequisites

```bash
# Activate environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline

# Or install locally
pip install -e ".[dev]"

# Required system tools for full test suite
module load plink/2.00-alpha-091019
module load bcftools
```

### Basic Run

```bash
# All tests
python -m pytest tests/ -v

# Specific file
python -m pytest tests/test_release_data.py -v
python -m pytest tests/test_exclusions.py -v
```

### Selective Execution

```bash
# Skip tests requiring external tools
python -m pytest tests/ -v -m "not (requires_bcftools or requires_plink2)"

# Only fast tests
python -m pytest tests/ -v -m "not slow"

# Specific test
python -m pytest tests/test_release_data.py::test_filter_correctness -v
```

### With Coverage

```bash
python -m pytest tests/ --cov=. --cov-report=term-missing --cov-report=html
```

## Test Categories

### 1. De-identification Integrity
Validates that all output IDs follow the de-identified format `{release_candid}{C|M}` with 10-digit release_candid.

**Key Tests**:
- `test_all_iids_match_deid_pattern` — Strictest check
- `test_fid_is_numeric` — Catches raw ID leaks
- `test_no_fid_zero_in_output` — No unmatched subjects

### 2. Filter Correctness
Re-derives the expected subject set from source data and compares to actual output.

**Key Test**:
- `test_filter_correctness` — End-to-end filter logic validation

### 3. Derivative Integrity
Ensures GRM, PC-AiR, PC-Relate, CNV outputs contain only release subjects with correct dimensions.

**Key Tests**:
- `test_pcrelate_grm_dimensions` — GRM matrix size matches sample count
- `test_pcair_32pcs_clean_ids_are_release` — PC-AiR subjects filtered correctly

### 4. Imputed VCF Validation
Validates filtered VCFs have correct samples, indices, and genotype concordance.

**Key Tests**:
- `test_imputed_vcf_subject_set_matches_release` — Exact IID set match
- `test_imputed_vcf_genotype_concordance` — PLINK/VCF genotype concordance

### 5. Exclusion Compliance
Verifies no excluded subject appears in any output file.

**Key Tests**:
- `test_hbcdcsv_excluded_pscids_absent` — HBCDexclusions.csv compliance
- `test_excel_excluded_release_candids_absent` — Excel exclusion compliance
- `test_removed_individuals_absent` — Pipeline exclusion list compliance

## Test Data Requirements

Tests require a completed Phase C release directory:

```bash
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData

# Tests auto-discover release dir via _lib.get_release_dir()
python -m pytest tests/ -v
```

### Required Files

| Test | Required Files |
|------|----------------|
| All | `keep_list.txt`, `temp.fam`, `GDA/merged_chroms.{bed,bim,fam}`, `GDA/batch.info` |
| Derivatives | `genesis/pcair_weights.tsv`, `genesis/pcrelate_relatedness.grm.*` |
| CNV | `cnv/CNV_slim.txt`, `cnv/CNV_bookmarks.csv` |
| VCF | `imputed/chr*_dose.vcf.gz`, `.tbi`, `.info.gz` |
| Concordance | `plink2`, `bcftools`, `bgzip`, `tabix` |

## Skipping Behavior

Tests skip gracefully with informative messages:

```python
# Missing file
if not p.exists():
    pytest.skip("genesis/pcair_weights.tsv not found")

# Missing external tool
if not _bcftools_available():
    pytest.skip("bcftools not on PATH")
```

## Writing New Tests

### Template

```python
def test_new_validation():
    """Description of what this test validates."""
    # Load required data
    fam, batch, excluded = _load_data()
    
    # Perform validation
    result = some_validation(fam)
    
    # Assert with helpful error message
    assert result.expected == result.actual, (
        f"Validation failed: {result.details}"
    )
```

### Best Practices

1. **Use helpers**: `_load_data()`, `_release_iids()`, `_release_n_subjects()`
2. **Skip gracefully**: Check file existence, tool availability
3. **Clear assertions**: Include actual values in error messages
4. **Mark appropriately**: `@pytest.mark.requires_bcftools`, `@pytest.mark.slow`

### Markers

```python
@pytest.mark.requires_bcftools
def test_vcf_something():
    ...

@pytest.mark.requires_plink2
def test_plink_something():
    ...

@pytest.mark.slow
def test_large_computation():
    ...
```

## Continuous Integration

### Local Pre-commit

```bash
# Install pre-commit
pip install pre-commit
pre-commit install

# Run manually
pre-commit run --all-files
```

### GitHub Actions

See [CI/CD Integration](scripts/test_scripts.md#ci-cd-integration) for example workflow.

## Troubleshooting Tests

### Common Failures

| Error | Cause | Solution |
|-------|-------|----------|
| `pytest.skip: merged_chroms.fam not found` | Phase C not run | Run `./run_release_pipeline.sh` first |
| `bcftools not on PATH` | Module not loaded | `module load bcftools` |
| `plink2 not on PATH` | Module not loaded | `module load plink/2.00-alpha-091019` |
| `AssertionError: X IIDs not in release set` | Filter failed | Check step 29 logs, re-run |
| `Genotype concordance failed` | Sample mismatch | Verify keep_list.txt, re-run step 27 |

### Debugging Tips

```bash
# Run single test with maximum output
python -m pytest tests/test_release_data.py::test_filter_correctness -v -s --tb=long

# Inspect test data
python -c "
from _lib import get_release_dir
import pandas as pd
r = get_release_dir()
print('Release dir:', r)
print('Files:', list(r.glob('**/*')))
"
```

## Performance

| Test Suite | Time (MSI) | Notes |
|------------|------------|-------|
| `test_release_data.py` | ~2 min | Most time in VCF concordance |
| `test_exclusions.py` | ~30 sec | Fast |
| Full suite | ~3 min | With all tools available |

Run without external tools for quick validation:
```bash
python -m pytest tests/ -v -m "not (requires_bcftools or requires_plink2)"
```