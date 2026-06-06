"""Ensure HBCDexclusions.csv pscids do not appear in the pipeline output."""

import pandas as pd
from pathlib import Path

from _lib import DATA_DIR, load_additional_excluded_pscids

RELEASE_DIR = Path(
    "/projects/standard/basu_hbcd/shared/HBCD_genomics_release_br_21p2/data/"
)


def test_excluded_pscids_absent_from_hbcd_fam():
    hbcd_fam_path = RELEASE_DIR / "hbcd.fam"
    if not hbcd_fam_path.exists():
        import pytest
        pytest.skip("hbcd.fam not found — run pipeline first")

    hbcd_fam = pd.read_csv(
        hbcd_fam_path,
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    fam_rc = set(
        hbcd_fam.loc[hbcd_fam["FID"] != 0, "FID"].dropna().unique()
    )

    excluded = load_additional_excluded_pscids()
    assert len(excluded) > 0, "No pscids loaded from HBCDexclusions.csv"

    identifiers = pd.read_csv(DATA_DIR / "release_identifiers_20260526.csv")
    identifiers = identifiers[identifiers["release_candid"] != "release_candid"]
    identifiers["release_candid"] = pd.to_numeric(identifiers["release_candid"])
    identifiers = identifiers.dropna(subset=["release_candid"])
    identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()

    exc_rc = set(
        identifiers.loc[identifiers["pscid"].isin(excluded), "release_candid"]
        .dropna()
        .astype(int)
        .unique()
    )

    bad = exc_rc & fam_rc
    assert len(bad) == 0, (
        f"{len(bad)} excluded pscid(s) map to release_candid in hbcd.fam: "
        f"{sorted(bad)[:20]}"
    )
