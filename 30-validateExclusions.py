#!/usr/bin/env python3

# step 30: validate no excluded subjects appear in any release output file
# Checks hbcd.fam + all derivative files + CNV for excluded-IID contamination.
#
# Usage: python 30-validateExclusions.py [--release-dir /path]

import sys
from pathlib import Path
import pandas as pd
from _lib import DATA_DIR, get_release_dir, load_additional_excluded_pscids

TOTAL_WIDTH = 67

# Accept --release-dir / -r override, default from env/HBCD_RELEASE
if len(sys.argv) > 1 and sys.argv[1] in ("--release-dir", "-r"):
    RELEASE_DIR = Path(sys.argv[2])
else:
    RELEASE_DIR = get_release_dir()

release_base = RELEASE_DIR.parent


def _sep(title=""):
    w = TOTAL_WIDTH
    if title:
        left = (w - len(title) - 2) // 2
        right = w - left - len(title) - 2
        print(f"\n{'=' * left}  {title}  {'=' * right}")
    else:
        print("=" * w)


def _ok(n_overlap):
    return "OK" if n_overlap == 0 else "OVERLAP"


# ══════════════════════════════════════════════════════════════════════════
# 1. Load exclusion sets
# ══════════════════════════════════════════════════════════════════════════

_sep("Loading exclusions")

excluded_pscids = load_additional_excluded_pscids()
print(f"  HBCDexclusions.csv pscids           : {len(excluded_pscids):>6}")

identifiers = pd.read_csv(DATA_DIR / "release_identifiers_20260628.csv")
identifiers = identifiers[identifiers["release_candid"] != "release_candid"]
identifiers["release_candid"] = pd.to_numeric(identifiers["release_candid"])
identifiers = identifiers.dropna(subset=["release_candid"])
identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()

exc_rc = set(
    identifiers.loc[identifiers["pscid"].isin(excluded_pscids), "release_candid"]
    .dropna()
    .astype(int)
    .unique()
)
print(f"  Mapped to release_candid            : {len(exc_rc):>6}")

# Build all possible excluded IIDs (both C and M for each RC)
exc_iids = set()
for rc in exc_rc:
    rc_s = str(rc)
    exc_iids.add(rc_s + "C")
    exc_iids.add(rc_s + "M")
print(f"  Derived excluded IIDs (C + M)       : {len(exc_iids):>6}")

# Load release IID set for comparison
try:
    keep = pd.read_csv(
        release_base / "keep_list.txt",
        sep=r"\s+", header=None, names=["FID", "IID"],
    )
    release_iids = set(keep["IID"].astype(str))
    print(f"  Release IIDs (from keep_list.txt)   : {len(release_iids):>6}")
except FileNotFoundError:
    release_iids = set()
    print("  Release IIDs: keep_list.txt not found — can't check reverse")


# ══════════════════════════════════════════════════════════════════════════
# 2. hbcd.fam — FID check (release_candid level)
# ══════════════════════════════════════════════════════════════════════════

_sep("hbcd.fam (release_candid-level check)")

fam_path = RELEASE_DIR / "hbcd.fam"
if fam_path.exists():
    hbcd = pd.read_csv(
        fam_path,
        sep=r"\s+", header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    fam_rc = set(hbcd.loc[hbcd["FID"] != 0, "FID"].dropna().unique())
    overlap = exc_rc & fam_rc
    s = _ok(len(overlap))
    print(f"  [{s:>7}]  FID (release_candid)")
    print(f"           unique RC in fam  : {len(fam_rc):>6}")
    print(f"           excluded RC in fam: {len(overlap):>6}")
    if overlap:
        print(f"           overlapping RC    : {sorted(overlap)[:20]}")

    # Also check IID level
    fam_iids = set(hbcd["IID"].astype(str))
    iid_overlap = fam_iids & exc_iids
    s2 = _ok(len(iid_overlap))
    n_exc = (hbcd["IID"].astype(str).isin(exc_iids)).sum()
    print(f"  [{s2:>7}]  IID")
    print(f"           total IIDs        : {len(fam_iids):>6}")
    print(f"           excluded IIDs     : {n_exc:>6}")
    if iid_overlap:
        print(f"           overlapping IIDs  : {sorted(iid_overlap)[:20]}")
else:
    print("  SKIP — hbcd.fam not found")


# ══════════════════════════════════════════════════════════════════════════
# 3. HBCDexclusions.csv — per-column check
# ══════════════════════════════════════════════════════════════════════════

_sep("HBCDexclusions.csv — per-column check")

if fam_path.exists():
    raw = pd.read_csv(DATA_DIR / "HBCDexclusions.csv")
    for col in raw.columns:
        pscids = set(raw[col].dropna().astype(str).str.strip())
        pscids = {p for p in pscids if p and p != "nan"}
        col_rc = set(
            identifiers.loc[identifiers["pscid"].isin(pscids), "release_candid"]
            .dropna()
            .astype(int)
            .unique()
        )
        in_fam = col_rc & fam_rc
        s = _ok(len(in_fam))
        print(f"  [{s:>7}]  {col}")
        print(f"           list size    : {len(pscids):>6}")
        print(f"           in hbcd.fam  : {len(in_fam):>6}")
        if in_fam:
            print(f"           overlapping RC: {sorted(in_fam)[:15]}")
else:
    print("  SKIP — hbcd.fam not found, skipping per-column check")


# ══════════════════════════════════════════════════════════════════════════
# 4. Derivative file checks (IID-level)
# ══════════════════════════════════════════════════════════════════════════

_sep("Derivative file checks (IID-level)")


def _check_file(path, label, id_cols, reader="csv"):
    """Read *path* and check that no excluded IID appears in *id_cols*.

    Parameters
    ----------
    path : Path
        Full path to the file.
    label : str
        Display name (e.g. ``hbcd_pcair_32PCs_clean.tsv``).
    id_cols : list of str
        Column name(s) containing IIDs to check.
    reader : str
        ``"csv"`` for ``pd.read_csv``, ``"fwf"`` for ``pd.read_fwf``-style
        (sep/headerless), or ``"tsv"`` for tab-separated.
    """
    if not path.exists():
        print(f"  [  SKIP]  {label}  — not found")
        return

    if reader == "fwf":
        df = pd.read_csv(path, sep=r"\s+", header=None, dtype=str)
    elif reader == "tsv":
        df = pd.read_csv(path, sep="\t", dtype=str)
    else:
        df = pd.read_csv(path, dtype=str)

    total = len(df)
    all_ids = set()
    for col in id_cols:
        if col in df.columns or (reader == "fwf" and col.isdigit()):
            if reader == "fwf":
                col_idx = int(col)
                vals = df[col_idx].dropna().astype(str)
            else:
                vals = df[col].dropna().astype(str)
            all_ids.update(vals)

    overlap = all_ids & exc_iids
    n_overlap = len(overlap)
    s = _ok(n_overlap)
    print(f"  [{s:>7}]  {label}")
    print(f"           total rows   : {total:>6}")
    print(f"           excluded IIDs: {n_overlap:>6}")
    if overlap:
        print(f"           overlapping  : {sorted(overlap)[:20]}")

    if release_iids:
        extra = all_ids - release_iids - exc_iids
        n_extra = len(extra)
        if n_extra:
            print(f"           WARN: {n_extra} IID(s) not in release set "
                  f"nor excluded: {sorted(extra)[:10]}")


# ── PC-Relate GRM binary ──────────────────────────────────────────────────

_check_file(
    RELEASE_DIR / "hbcd_pcrelate_grm.id",
    "hbcd_pcrelate_grm.id",
    ["0", "1"],  # FID, IID — headerless
    reader="fwf",
)

# ── PC-AiR 32 PCs clean ──────────────────────────────────────────────────

_check_file(
    RELEASE_DIR / "hbcd_pcair_32PCs_clean.tsv",
    "hbcd_pcair_32PCs_clean.tsv",
    ["subject_id"],
    reader="tsv",
)

# ── PC-Relate GRM pairwise ───────────────────────────────────────────────

_check_file(
    RELEASE_DIR / "hbcd_pcrelate_grm_pairwise.tsv",
    "hbcd_pcrelate_grm_pairwise.tsv",
    ["ID1", "ID2"],
    reader="tsv",
)

# ── CNV slim clean ───────────────────────────────────────────────────────

_check_file(
    RELEASE_DIR / "CNV_slim_clean.txt",
    "CNV_slim_clean.txt",
    ["sample_id"],
    reader="tsv",
)

# ══════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════

_sep()
print("Validation complete.")
print(f"  Excluded pscids : {len(excluded_pscids)}")
print(f"  Excluded RC     : {len(exc_rc)}")
print(f"  Excluded IIDs   : {len(exc_iids)}")
