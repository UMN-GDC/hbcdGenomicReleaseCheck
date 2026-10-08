# Test Scripts Reference

## tests/test_release_data.py

**Purpose**: Validate de-identification integrity, data integrity, filter correctness, and per-file IID checks for release derivative outputs.

**Run**:
```bash
python -m pytest tests/test_release_data.py -v
```

**Dependencies**: `pytest`, `pandas`, `numpy`, `bcftools` (for VCF tests), `plink2` (for concordance test)

### Test Categories

#### De-identification Integrity (8 tests)

| Test | Description |
|------|-------------|
| `test_no_fid_zero_in_output` | No unmatched subjects (FID=0) in release output |
| `test_all_iids_match_deid_pattern` | Every IID matches `^\d{10}[CM]$` |
| `test_all_fids_are_positive_10_digit_integers` | FIDs are positive 10-digit integers |
| `test_iid_first_10_digits_equal_fid` | IID prefix = FID |
| `test_relation_is_C_or_M` | Suffix is C or M |
| `test_fid_length_is_10` | Strict 10-digit FID check |
| `test_iid_length_is_11` | IID length = 11 |
| `test_fid_is_numeric` | FID is numeric (no raw ID leak) |

#### Row-Count Integrity (1 test)

| Test | Description |
|------|-------------|
| `test_temp_fam_row_count_matches_onlyqc` | temp.fam rows = onlyQc.fam rows (PLINK .bed size match) |

#### fam/batch Consistency (2 tests)

| Test | Description |
|------|-------------|
| `test_fam_and_batch_order_matches` | fam/batch IID sets + order match 1:1 |
| `test_batch_info_iids_are_deidentified` | batch.info IIDs match de-ID pattern |

#### par_visit + Exclusion Constraint (1 test)

| Test | Description |
|------|-------------|
| `test_all_output_iids_in_par_visit` | All output IIDs map to par_visit candid AND not excluded |

#### Filter Correctness (1 test)

| Test | Description |
|------|-------------|
| `test_filter_correctness` | Re-derives expected subject set from scratch, confirms match with merged_chroms.fam |

#### Variant Integrity (1 test)

| Test | Description |
|------|-------------|
| `test_variant_count_preserved` | .bim variant count unchanged (PLINK --keep preserves variants) |

#### Derivative Integrity (6 tests)

| Test | Description |
|------|-------------|
| `test_pcrelate_grm_ids_are_release` | All IIDs in GRM .id are release IIDs |
| `test_pcrelate_grm_dimensions` | GRM .id row count = merged_chroms.fam rows |
| `test_pcair_32pcs_clean_ids_are_release` | All subject_ids in pcair_weights.tsv are release |
| `test_pcrelate_grm_pairwise_ids_are_release` | All ID1/ID2 in pairwise TSV are release |
| `test_pcrelate_grm_gz_ids_are_release` | All IIDs in GRM .gz (via index mapping) are release |
| `test_cnv_slim_clean_ids_are_release` | All CNV sample_ids are release |
| `test_cnv_bookmark_ids_are_release` | All bookmark sample_ids are release |

#### Imputed VCF Integrity (7 tests)

| Test | Description |
|------|-------------|
| `test_imputed_vcf_sample_ids_deidentified` | VCF samples match `^\d{10}[CM]$` |
| `test_imputed_vcf_subject_count_matches_release` | Per-chr sample count = release total |
| `test_imputed_vcf_subject_set_matches_release` | VCF IID set = release IIDs (no missing/extra) |
| `test_imputed_output_completeness` | All expected VCF outputs exist |
| `test_release_tree_completeness` | Full release directory tree exists |
| `test_imputation_qc_files_content` | QC stats files non-empty, correct columns |
| `test_imputed_vcf_genotype_concordance` | PLINK/VCF genotype concordance (chr22, 3 subjects × 5 variants) |

### Helper Functions

```python
def _load_data():
    """Load merged_chroms.fam, batch.info, removed_individuals.txt"""
    ...

def _load_temp_fam():
    """Load temp.fam from release base"""
    ...

def _prepare_iid_test(fam, batch, excluded):
    """Prepare IID test dataframe with parsed FID/IID"""
    ...

def _id_lengths(iid_test):
    """Compute FID/IID length distributions"""
    ...

def _release_n_subjects():
    """Number of subjects in merged_chroms.fam"""
    ...

def _release_iids():
    """Memoized release IID set from keep_list.txt"""
    ...
```

### Skipping Behavior

Tests skip gracefully if required files missing:
```python
if not p.exists():
    pytest.skip("genesis/pcair_weights.tsv not found")
```

Tests requiring external tools skip if unavailable:
```python
def _bcftools_available():
    try:
        subprocess.run(["bcftools", "--version"], capture_output=True)
        return True
    except FileNotFoundError:
        return False
```

---

## tests/test_exclusions.py

**Purpose**: Ensure all exclusion sources properly applied to final output.

**Run**:
```bash
python -m pytest tests/test_exclusions.py -v
```

### Tests (4)

| Test | Description |
|------|-------------|
| `test_hbcdcsv_excluded_pscids_absent` | No HBCDexclusions.csv pscid leaks into merged_chroms.fam |
| `test_each_hbcdcsv_exclusion_reason_individually` | Each HBCDexclusions column checked separately |
| `test_excel_excluded_release_candids_absent` | No Excel-excluded release_candid leaks |
| `test_removed_individuals_absent` | No GDA/removed_individuals.txt IID leaks |

### Logic

```python
# 1. Load merged_chroms.fam
fam_rc = set(fam[fam.FID != 0].FID.unique())
fam_iids = set(fam.IID.astype(str))

# 2. HBCDexclusions.csv (pscid-level)
excluded_pscids = load_additional_excluded_pscids()
exc_rc = identifiers[identifiers.pscid.isin(excluded_pscids)].release_candid
overlap = fam_rc ∩ exc_rc
assert len(overlap) == 0

# 3. Per-column check
for col in HBCDexclusions.columns:
    col_pscids = set(raw[col].dropna())
    col_rc = identifiers[identifiers.pscid.isin(col_pscids)].release_candid
    overlap = fam_rc ∩ col_rc
    assert len(overlap) == 0

# 4. Excel exclusions (release_candid-level)
exc_rc_excel = load_excluded_release_candids()
overlap = fam_rc ∩ exc_rc_excel
assert len(overlap) == 0

# 5. Pipeline's removed_individuals.txt
removed_iids = load_removed_individuals()
overlap = fam_iids ∩ removed_iids
assert len(overlap) == 0
```

---

## tests/conftest.py

**Purpose**: Pytest configuration and shared fixtures.

**Contents**:
```python
import pytest
from _lib import get_release_dir, DATA_DIR

@pytest.fixture(scope="session")
def release_dir():
    return get_release_dir()

@pytest.fixture(scope="session")
def data_dir():
    return DATA_DIR

@pytest.fixture(scope="session")
def release_iids():
    """Release IID set from keep_list.txt"""
    keep = pd.read_csv(
        get_release_dir().parent / "keep_list.txt",
        sep=r"\s+", header=None, names=["FID", "IID"]
    )
    return set(keep["IID"].astype(str))
```

---

## Running Tests

### Full Suite
```bash
python -m pytest tests/ -v
```

### Specific Test File
```bash
python -m pytest tests/test_release_data.py -v
python -m pytest tests/test_exclusions.py -v
```

### Filter by Marker
```bash
# Skip slow tests
python -m pytest tests/ -v -m "not slow"

# Only VCF tests (requires bcftools)
python -m pytest tests/ -v -m "requires_bcftools"

# Only PLINK tests (requires plink2)
python -m pytest tests/ -v -m "requires_plink2"
```

### With Coverage
```bash
python -m pytest tests/ --cov=. --cov-report=html
```

### Parallel Execution
```bash
python -m pytest tests/ -v -n auto  # Requires pytest-xdist
```

---

## Test Data Requirements

Tests require a completed Phase C run with all outputs present:

```
$RELEASE_DIR/
├── GDA/
│   ├── merged_chroms.{bed,bim,fam}
│   ├── batch.info
│   └── removed_individuals.txt
├── genesis/
│   ├── pcair_weights.tsv
│   ├── pcrelate_relatedness.grm.{id,bin,N.bin,gz,tsv}
├── cnv/
│   ├── CNV_slim.txt
│   └── CNV_bookmarks.csv
├── imputed/
│   ├── chr{1-22}_dose.vcf.gz + .tbi
│   ├── chrX_dose.vcf.gz + .tbi
│   ├── chr{1-22}.info.gz
│   ├── batch*-snps-typed-only.txt
│   ├── batch*-snps-excluded.txt
│   ├── batch*-chunks-excluded.txt
│   └── batch*-quality-control.html
├── keep_list.txt
└── temp.fam
```

### Generating Test Data (Development)

For local development without full pipeline run:

```bash
# Create minimal test fixtures
mkdir -p tests/fixtures/release/GDA
mkdir -p tests/fixtures/release/genesis
mkdir -p tests/fixtures/release/cnv
mkdir -p tests/fixtures/release/imputed

# Generate synthetic data matching schemas
python tests/generate_fixtures.py
```

(Fixture generation script not yet implemented)

---

## CI/CD Integration

### GitHub Actions Example

```yaml
# .github/workflows/test.yml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.10'
      - name: Install dependencies
        run: pip install -e ".[dev]"
      - name: Install system tools
        run: |
          sudo apt-get update
          sudo apt-get install -y plink2 bcftools
      - name: Run tests
        run: python -m pytest tests/ -v
        env:
          HBCD_DATA_DIR: ${{ secrets.HBCD_DATA_DIR }}
          HBCD_IMPUTATION_DIR: ${{ secrets.HBCD_IMPUTATION_DIR }}
          HBCD_RELEASE: br_test
```

---

## Test Maintenance

### Adding New Tests

1. Add test function to appropriate file (`test_release_data.py` or `test_exclusions.py`)
2. Follow naming convention: `test_<description>`
3. Use `_load_data()`, `_release_iids()` helpers
4. Skip gracefully if files missing
5. Add marker if requires external tools (`@pytest.mark.requires_bcftools`)

### Updating Tests for New Release

When release schema changes:
1. Update expected column names in filter tests
2. Adjust file path expectations in completeness tests
3. Update concordance test parameters if needed
4. Regenerate test fixtures if applicable