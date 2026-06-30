"""Ensure all exclusion sources are properly applied to the final output.

Three independent exclusion mechanisms are checked:
  1. HBCDexclusions.csv        — pscid-level exclusion (multi-column reasons)
   2. Excel exclusion list          — release_candid-level exclusion
   3. GDA/removed_individuals.txt   — IID-level exclusion written by the pipeline

Set ``HBCD_RELEASE`` env var to test against a different release.
"""
import pandas as pd
import pytest

from _lib import (
    DATA_DIR,
    get_release_dir,
    get_release_base,
    load_additional_excluded_pscids,
    load_excluded_release_candids,
    load_identifiers,
)

RELEASE_DIR = get_release_dir()
GDA_DIR = RELEASE_DIR / "GDA"


def _load_hbcd_fam():
    p = GDA_DIR / "merged_chroms.fam"
    if not p.exists():
        pytest.skip("merged_chroms.fam not found — run pipeline first")
    return pd.read_csv(
        p, sep=r"\s+", header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )


def _fam_release_candids(fam):
    """Return the set of non-zero FID values (= release_candids)."""
    return set(fam.loc[fam["FID"] != 0, "FID"].dropna().astype(int).unique())


def _fam_iids(fam):
    return set(fam["IID"].astype(str))


# ── 1. HBCDexclusions.csv (pscid-level) ─────────────────────────────────


def test_hbcdcsv_excluded_pscids_absent():
    """Every pscid listed in HBCDexclusions.csv maps to a release_candid
    that is NOT present in GDA/merged_chroms.fam."""
    hbcd = _load_hbcd_fam()
    fam_rc = _fam_release_candids(hbcd)

    excluded_pscids = load_additional_excluded_pscids()
    assert len(excluded_pscids) > 0, "HBCDexclusions.csv is empty"

    ids = load_identifiers()
    exc_rc = set(
        ids.loc[ids["pscid"].isin(excluded_pscids), "release_candid"]
        .dropna().astype(int).unique()
    )

    overlap = fam_rc & exc_rc
    assert len(overlap) == 0, (
        f"{len(overlap)} HBCDexclusions pscid(s) map to release_candid "
        f"found in GDA/merged_chroms.fam: {sorted(overlap)[:20]}"
    )


def test_each_hbcdcsv_exclusion_reason_individually():
    """Each column of HBCDexclusions.csv is checked separately — no single
    exclusion reason should leak a subject into the output."""
    hbcd = _load_hbcd_fam()
    fam_rc = _fam_release_candids(hbcd)
    ids = load_identifiers()
    raw = pd.read_csv(DATA_DIR / "HBCDexclusions.csv")

    for col in raw.columns:
        pscids = set(raw[col].dropna().astype(str).str.strip())
        pscids = {p for p in pscids if p and p != "nan"}
        if not pscids:
            continue
        col_rc = set(
            ids.loc[ids["pscid"].isin(pscids), "release_candid"]
            .dropna().astype(int).unique()
        )
        overlap = fam_rc & col_rc
        assert len(overlap) == 0, (
            f"{len(overlap)} subject(s) from HBCDexclusions column "
            f"'{col}' present in GDA/merged_chroms.fam: {sorted(overlap)[:15]}"
        )


# ── 2. Excel exclusion list (release_candid-level) ──────────────────────


def test_excel_excluded_release_candids_absent():
    """Every release_candid marked for exclusion in the Excel spreadsheet
    is absent from GDA/merged_chroms.fam."""
    hbcd = _load_hbcd_fam()
    fam_rc = _fam_release_candids(hbcd)

    exc_rc = load_excluded_release_candids()
    assert len(exc_rc) > 0, "Excel exclusion list is empty"

    overlap = fam_rc & exc_rc
    assert len(overlap) == 0, (
        f"{len(overlap)} Excel-excluded release_candid(s) "
        f"present in GDA/merged_chroms.fam: {sorted(overlap)[:20]}"
    )


# ── 3. GDA/removed_individuals.txt (IID-level, written by pipeline) ──────


def test_removed_individuals_absent():
    """Every IID listed in GDA/removed_individuals.txt is absent from merged_chroms.fam."""
    hbcd = _load_hbcd_fam()
    fam_iids = _fam_iids(hbcd)

    p = GDA_DIR / "removed_individuals.txt"
    if not p.exists():
        pytest.skip("GDA/removed_individuals.txt not found")
    removed = pd.read_csv(p, delim_whitespace=True, header=None, names=["IID"])
    removed_iids = set(removed["IID"].astype(str))

    overlap = fam_iids & removed_iids
    assert len(overlap) == 0, (
        f"{len(overlap)} removed IID(s) present in merged_chroms.fam: "
        f"{sorted(overlap)[:20]}"
    )
