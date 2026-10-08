# Contributing Guide

## Development Setup

### Fork & Clone
```bash
git clone https://github.com/your-username/hbcd-genomic-release.git
cd hbcd-genomic-release
git remote add upstream https://github.com/hbcd-genomics/hbcd-genomic-release.git
```

### Environment
```bash
# Create development environment
conda create -n hbcd-genomic-dev python=3.10
conda activate hbcd-genomic-dev
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

### Code Style
```bash
# Format with black
black .

# Lint with ruff
ruff check .

# Type check with mypy
mypy .

# Run tests
pytest tests/ -v
```

## Pull Request Process

### 1. Create Feature Branch
```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/issue-description
```

### 2. Make Changes
- Write code following existing patterns
- Add tests for new functionality
- Update documentation
- Run full test suite

### 3. Commit Guidelines
```bash
# Conventional commit format
git commit -m "feat: add new CNV filtering option"
git commit -m "fix: correct GRM binary filter indexing"
git commit -m "docs: update installation guide for MSI"
git commit -m "test: add VCF concordance edge case test"
```

### 4. Push & Create PR
```bash
git push origin feature/your-feature-name
# Create PR via GitHub UI
```

### 5. PR Requirements
- [ ] All tests pass (`pytest tests/ -v`)
- [ ] Code formatted (`black . && ruff check .`)
- [ ] Type checks pass (`mypy .`)
- [ ] Documentation updated
- [ ] CHANGELOG.md updated (Unreleased section)
- [ ] Version bump if applicable

## Code Standards

### Python
- **Formatter**: Black (line length 100)
- **Linter**: Ruff (with pyflakes, pycodestyle, isort, bugbear)
- **Type Hints**: Required for public functions (mypy strict mode)
- **Docstrings**: NumPy style for all public functions/classes

```python
def filter_subjects(df: pd.DataFrame, keep_ids: set[str], id_col: str = "IID") -> pd.DataFrame:
    """Filter DataFrame to keep only specified subject IDs.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with subject ID column.
    keep_ids : set[str]
        Set of subject IDs to retain.
    id_col : str, default "IID"
        Column name containing subject IDs.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only keep_ids subjects.
    """
    return df[df[id_col].isin(keep_ids)].copy()
```

### R Scripts
- Use `tidyverse` style
- Document with `roxygen2` comments
- Use explicit library calls (`library(tidyverse)`)

### Shell Scripts
- Use `set -euo pipefail`
- Document expected environment variables
- Use descriptive variable names

### SLURM Scripts
- Include resource requirements (`--mem`, `--time`, `--partition`)
- Use arrays for parallelization
- Document task-to-chromosome mapping

## Testing Requirements

### New Features
- Add unit tests in `tests/test_release_data.py` or `test_exclusions.py`
- Add integration tests if external tools required
- Mark appropriately: `@pytest.mark.requires_bcftools`, `@pytest.mark.slow`

### Bug Fixes
- Add regression test demonstrating the bug
- Ensure test fails before fix, passes after

### Test Data
- Use fixtures for synthetic test data
- Don't commit real HBCD data
- Document fixture generation in `tests/generate_fixtures.py` (future)

## Documentation Standards

### Inline Documentation
- All public functions: NumPy docstring
- Script headers: Purpose, usage, inputs, outputs
- Complex logic: Inline comments explaining "why"

### User-Facing Docs
- Update relevant `.md` files in `docs/source/`
- Keep examples runnable
- Update version-specific paths

### Changelog
Add entries to `CHANGELOG.md` under `## [Unreleased]`:

```markdown
## [Unreleased]

### Added
- New CNV filtering option for mosaic events

### Fixed
- GRM binary filter off-by-one indexing error

### Changed
- Updated PLINK2 minimum version to 2.00-alpha-2023
```

## Versioning

### Semantic Versioning
- **MAJOR**: Breaking changes to outputs or API
- **MINOR**: New features, backward compatible
- **PATCH**: Bug fixes, backward compatible

### Release Process
```bash
# 1. Update version in pyproject.toml and hbcd_genomic_release/_version.py
# 2. Update CHANGELOG.md with release date
# 3. Create git tag
git tag -a v0.2.0 -m "Release v0.2.0"
git push origin v0.2.0
# 4. GitHub Action builds and publishes to PyPI (future)
# 5. Create Zenodo release (future)
```

## HBCD-Specific Guidelines

### Data Privacy
- **Never commit real participant data**
- Use synthetic fixtures for testing
- All outputs must be de-identified (`^\d{10}[CM]$`)

### Release Compatibility
- Test against target HBCD release (e.g., `br31p2`)
- Document required input file versions
- Maintain backward compatibility where possible

### Review Process
- Code review by at least one Genomics Team member
- HDCC liaison review for release-critical changes
- NMIND checklist review for major versions (future)

## Issue Reporting

### Bug Reports
Include:
1. Pipeline step and script name
2. Environment (MSI, local, container)
3. HBCD release tag
4. Full error message + stack trace
5. Relevant log files
6. Steps to reproduce

### Feature Requests
Include:
1. Use case description
2. Expected behavior
3. Affected pipeline steps
4. Backward compatibility considerations

## Code of Conduct

This project follows the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/). By participating, you agree to uphold this code.

## License

By contributing, you agree that your contributions will be licensed under the Apache License 2.0 (see [LICENSE](../license.md)).