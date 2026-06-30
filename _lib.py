"""Shared helpers for the HBCD genomics release pipeline."""
import os
from pathlib import Path

import pandas as pd

DATA_DIR = Path("/projects/standard/basu_hbcd/shared/data")
DEFAULT_RELEASE = "br_21p3"

# If HBCD_RELEASE env var is set, all scripts use that release
_RELEASE_ENV = os.environ.get("HBCD_RELEASE", DEFAULT_RELEASE)


def get_release_dir(release=None):
    """Return the release genotype_microarray directory, from env var or explicit argument.

    The directory is ``<DATA_DIR.parent>/HBCD_genomics_release_<release>/genotype_microarray/``.
    If *release* is ``None`` (default) the ``HBCD_RELEASE`` environment
    variable is consulted, falling back to ``DEFAULT_RELEASE``.
    """
    if release is None:
        release = _RELEASE_ENV
    return DATA_DIR.parent / f"HBCD_genomics_release_{release}" / "genotype_microarray"


def get_release_base(release=None):
    """Return the release base directory (parent of ``genotype_microarray/``)."""
    return get_release_dir(release).parent


PAR_VISIT_FILE = "par_visit_data_br21_1.tsv"
IDENTIFIERS_FILE = "release_identifiers_20260628.csv"
EXCLUSION_FILE = "HBCD_genetics_QC1_missing_race_LORIS.xlsx"


def load_par_visit_candids(path=None):
    """Return a set of integral participant IDs (*candid*) from par_visit.

    Filters to completed visits (``par_visit_data_visit_missed == 'No'``)
    and strips the ``sub-`` prefix from ``participant_id`` per
    ``int(val[4:])``.
    """
    if path is None:
        path = DATA_DIR / PAR_VISIT_FILE
    pv = pd.read_csv(path, sep="\t").query("par_visit_data_visit_missed == 'No'")
    def _to_int(v):
        s = str(v)
        return int(s[4:] if s.startswith("sub-") else s)
    ids = pv["participant_id"].dropna().apply(_to_int).unique()
    return set(int(x) for x in ids)


def load_identifiers(path=None):
    """Load the release-identifiers CSV, parse numeric columns, return DataFrame."""
    if path is None:
        path = DATA_DIR / IDENTIFIERS_FILE
    df = (
        pd.read_csv(path)
        .query("release_candid != 'release_candid'")  # drop placeholder rows
        .assign(release_candid=lambda x: pd.to_numeric(x["release_candid"]))
        .assign(candid=lambda x: pd.to_numeric(x["candid"]))
        .dropna(subset=["release_candid"])
        .drop_duplicates(subset=["pscid", "candid", "release_candid"])
    )
    return df


def load_excluded_release_candids(path=None):
    """Return a set of ``release_candid`` (int) from the Excel exclusion list."""
    if path is None:
        path = DATA_DIR / EXCLUSION_FILE
    exc = pd.read_excel(path, sheet_name="Exclude_Summary")
    exc[["_", "release_candid", "_pscid"]] = exc["Sampled ID"].str.split(
        "_", n=2, expand=True
    )
    exc["release_candid"] = pd.to_numeric(exc["release_candid"])
    exc = exc.dropna(subset=["release_candid"])
    return set(int(x) for x in exc["release_candid"].unique())


def load_excluded_with_relationship(path=None):
    """Return a DataFrame with ``release_candid`` and ``relationship`` from
    the Excel exclusion list, for constructing IIDs of removed individuals."""
    if path is None:
        path = DATA_DIR / EXCLUSION_FILE
    exc = pd.read_excel(path, sheet_name="Exclude_Summary")
    exc[["_", "release_candid", "_pscid"]] = exc["Sampled ID"].str.split(
        "_", n=2, expand=True
    )
    exc["release_candid"] = pd.to_numeric(exc["release_candid"])
    exc["relationship"] = exc["Study_ID"].str[-1]
    return exc.dropna(subset=["release_candid"])[["release_candid", "relationship"]]


ADDITIONAL_EXCLUSIONS_FILE = "HBCDexclusions.csv"


def load_additional_excluded_pscids():
    """Load HBCDexclusions.csv and return a set of all pscids to exclude.

    CSV has exclusion reasons as column names and pscids as values.
    """
    raw = pd.read_csv(DATA_DIR / ADDITIONAL_EXCLUSIONS_FILE)
    excluded = set()
    for col in raw.columns:
        vals = raw[col].dropna().astype(str).str.strip()
        excluded.update(vals[vals != "nan"])
    return excluded
