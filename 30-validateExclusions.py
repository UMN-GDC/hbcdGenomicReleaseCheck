#!/usr/bin/env python3

# step 30: validate no excluded subjects appear in any release output file
# Checks GDA/merged_chroms.fam + all derivative files + CNV for excluded-IID contamination.
#
# Usage: python 30-validateExclusions.py [--release-dir /path]

import sys
from pathlib import Path
import pandas as pd
from datetime import datetime
from _lib import DATA_DIR, get_release_dir, load_additional_excluded_pscids

TOTAL_WIDTH = 67

# Accept --release-dir / -r override, default from env/HBCD_RELEASE
if len(sys.argv) > 1 and sys.argv[1] in ("--release-dir", "-r"):
    RELEASE_DIR = Path(sys.argv[2])
else:
    RELEASE_DIR = get_release_dir()

release_base = RELEASE_DIR.parent
gda_dir = RELEASE_DIR / "GDA"
genesis_dir = RELEASE_DIR / "genesis"
cnv_dir = RELEASE_DIR / "cnv"

# ── tee stdout to log file ──────────────────────────────────────────────
class _Tee:
    def __init__(self, path):
        self.file = open(path, "w", buffering=1)
    def write(self, data):
        sys.__stdout__.write(data)
        self.file.write(data)
    def flush(self):
        sys.__stdout__.flush()
        self.file.flush()

release_base.mkdir(parents=True, exist_ok=True)
_log_path = release_base / f"log_30_validateExclusions_{datetime.now():%Y%m%d_%H%M%S}.txt"
_tee = _Tee(_log_path)
sys.stdout = _tee
print(f"Log: {_log_path}")


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
# 2. GDA/merged_chroms.fam — FID check (release_candid level)
# ══════════════════════════════════════════════════════════════════════════

_sep("GDA/merged_chroms.fam (release_candid-level check)")

fam_path = gda_dir / "merged_chroms.fam"
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
        print(f"           in merged_chroms.fam: {len(in_fam):>6}")
        if in_fam:
            print(f"           overlapping RC: {sorted(in_fam)[:15]}")
else:
    print("  SKIP — merged_chroms.fam not found, skipping per-column check")


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


# ── PC-Relate GRM binary (genesis/) ────────────────────────────────────────

_check_file(
    genesis_dir / "pcrelate_relatedness.grm.id",
    "genesis/pcrelate_relatedness.grm.id",
    ["1"],  # IID only (col 1) — FID is release_candid, not an IID
    reader="fwf",
)

# ── PC-AiR 32 PCs clean (genesis/) ────────────────────────────────────────

_check_file(
    genesis_dir / "pcair_weights.tsv",
    "genesis/pcair_weights.tsv",
    ["subject_id"],
    reader="tsv",
)

# ── PC-Relate GRM pairwise (genesis/) ─────────────────────────────────────

_check_file(
    genesis_dir / "pcrelate_relatedness.tsv",
    "genesis/pcrelate_relatedness.tsv",
    ["ID1", "ID2"],
    reader="tsv",
)

# ── PC-Relate GRM gzipped text (genesis/) ────────────────────────────────

# Uses 1-based indices into .grm.id; map through the ID file
_grm_gz_path = genesis_dir / "pcrelate_relatedness.grm.gz"
_grm_id_path = genesis_dir / "pcrelate_relatedness.grm.id"
if _grm_gz_path.exists():
    if not _grm_id_path.exists():
        print("  [  SKIP]  genesis/pcrelate_relatedness.grm.gz — .grm.id not found")
    else:
        grm_ids = pd.read_csv(
            _grm_id_path, sep=r"\s+", header=None,
            names=["FID", "IID"], dtype=str,
        )
        iid_list = grm_ids["IID"].tolist()
        df_gz = pd.read_csv(
            _grm_gz_path, sep=r"\s+", header=None,
            names=["IID1", "IID2", "N", "GRM"],
            dtype={"IID1": int, "IID2": int, "N": int, "GRM": float},
        )
        gz_iids = set()
        max_idx1 = df_gz["IID1"].max()
        max_idx2 = df_gz["IID2"].max()
        for _, row in df_gz.iterrows():
            idx1, idx2 = int(row["IID1"]), int(row["IID2"])
            if 1 <= idx1 <= len(iid_list):
                gz_iids.add(iid_list[idx1 - 1])
            if 1 <= idx2 <= len(iid_list):
                gz_iids.add(iid_list[idx2 - 1])
        overlap = gz_iids & exc_iids
        s = _ok(len(overlap))
        print(f"  [{s:>7}]  genesis/pcrelate_relatedness.grm.gz")
        print(f"           total rows       : {len(df_gz):>6}")
        print(f"           max idx1/idx2    : {max_idx1} / {max_idx2}")
        print(f"           grm.id subjects  : {len(iid_list):>6}")
        print(f"           excluded IIDs    : {len(overlap):>6}")
        if overlap:
            print(f"           overlapping      : {sorted(overlap)[:20]}")
        if release_iids:
            extra = gz_iids - release_iids - exc_iids
            if extra:
                print(f"           WARN: {len(extra)} IID(s) not in release "
                      f"set nor excluded: {sorted(extra)[:10]}")
else:
    print("  [  SKIP]  genesis/pcrelate_relatedness.grm.gz — not found")

# ── CNV slim clean (cnv/) ────────────────────────────────────────────────

_check_file(
    cnv_dir / "CNV_slim.txt",
    "cnv/CNV_slim.txt",
    ["sample_id"],
    reader="tsv",
)

# ── CNV bookmark metrics (cnv/) ──────────────────────────────────────────

_check_file(
    cnv_dir / "HBCD_CNV_bookmark_metrics_clean.csv",
    "cnv/HBCD_CNV_bookmark_metrics_clean.csv",
    ["sample_id"],
    reader="csv",
)

# ══════════════════════════════════════════════════════════════════════════
# 5. Cross-form subject count consistency
# ══════════════════════════════════════════════════════════════════════════

_sep("Cross-form subject count consistency")

def _unique_iids_from(path, id_col, sep="\t"):
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep=sep, dtype=str)
    return set(df[id_col].dropna().astype(str))

derivatives = []
if genesis_dir.exists():
    for f in genesis_dir.glob("*"):
        derivatives.append(("genesis/" + f.name, f))

if cnv_dir.exists():
    for f in cnv_dir.glob("*"):
        derivatives.append(("cnv/" + f.name, f))

for label, path in derivatives:
    if path.suffix in (".bin", ".N.bin", ".gz", ".zip"):
        continue
    if ".grm.id" in path.name:
        df = pd.read_csv(path, sep=r"\s+", header=None, dtype=str)
        ids = set(df[1].dropna().astype(str))
        iid_set = set(ids)
    elif path.name == "CNV_slim.txt":
        iid_set = _unique_iids_from(path, "sample_id", sep="\t")
    elif path.name == "HBCD_CNV_bookmark_metrics_clean.csv":
        iid_set = _unique_iids_from(path, "sample_id", sep=",")
    elif path.name.endswith(".tsv"):
        df = pd.read_csv(path, sep="\t", dtype=str)
        ids = set()
        for col in df.columns:
            if col in ("ID1", "ID2", "subject_id"):
                ids.update(df[col].dropna().astype(str))
        iid_set = set(ids)
    elif path.name.endswith(".id"):
        df = pd.read_csv(path, sep=r"\s+", header=None, dtype=str)
        ids = set(df[1].dropna().astype(str))
        iid_set = set(ids)
    else:
        continue

    n = len(iid_set)
    match_release = "✓" if iid_set == release_iids else ("DIFF" if release_iids else "?")
    print(f"  [{match_release}]  {label}")
    print(f"           unique IIDs     : {n:>6}")

# ══════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════

_sep()
print("Validation complete.")
n_fam = len(fam_iids) if fam_path.exists() else 0
print(f"  Release IIDs (merged_chroms.fam): {n_fam}")
print(f"  Release IIDs (keep_list.txt)    : {len(release_iids)}")
print(f"  Excluded pscids                 : {len(excluded_pscids)}")
print(f"  Excluded RC                     : {len(exc_rc)}")
print(f"  Excluded IIDs                   : {len(exc_iids)}")
