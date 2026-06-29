#!/usr/bin/env python3

# step 29: filter all de-identified handoff derivatives to release subjects
# Reads keep_list.txt (release IIDs) + de-IDed CNV from data_handoff/,
# writes filtered copies to the release data directory.
#
# Filtered outputs:
#   - PC-Relate GRM (binary)  →  release_dir/hbcd_pcrelate_grm.{grm.id,grm.bin,grm.N.bin}
#   - PC-Relate GRM (text gz) →  release_dir/hbcd_pcrelate_grm.gz
#   - PC-Relate GRM pairwise  →  release_dir/hbcd_pcrelate_grm_pairwise.tsv
#   - PC-AiR 32 PCs clean     →  release_dir/hbcd_pcair_32PCs_clean.tsv
#   - CNV slim clean          →  release_dir/CNV_slim_clean.txt
#
# Usage: python 29-filter_release_outputs.py [--release-dir /path]
# Default release from HBCD_RELEASE env var or br_21p2.

import sys
import pandas as pd
import numpy as np
from pathlib import Path
from _lib import DATA_DIR, get_release_dir

HANDOFF_DIR = DATA_DIR.parent / "data_handoff"

if len(sys.argv) > 1 and sys.argv[1] in ("--release-dir", "-r"):
    RELEASE_DIR = Path(sys.argv[2])
else:
    RELEASE_DIR = get_release_dir()

release_base = RELEASE_DIR.parent

keep = pd.read_csv(
    release_base / "keep_list.txt",
    sep=r"\s+", header=None, names=["FID", "IID"],
)
release_iids = set(keep["IID"].astype(str))
print(f"Release IIDs: {len(release_iids)}")


# ── filter helpers ─────────────────────────────────────────────────────


def filter_text(in_path, id_col, out_path):
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    before = len(df)
    df = df[df[id_col].isin(release_iids)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  {in_path.name}: {before} → {len(df)} rows")


def filter_binary_grm(in_prefix, out_prefix, id_ext=".grm.id", bin_ext=".grm.bin", n_ext=".grm.N.bin"):
    id_in = Path(f"{in_prefix}{id_ext}")
    bin_in = Path(f"{in_prefix}{bin_ext}")
    n_in = Path(f"{in_prefix}{n_ext}")
    id_out = Path(f"{out_prefix}{id_ext}")
    bin_out = Path(f"{out_prefix}{bin_ext}")
    n_out = Path(f"{out_prefix}{n_ext}")

    if not id_in.exists() or not bin_in.exists():
        print(f"  SKIP (not found): {in_prefix}.grm.*")
        return

    ids = pd.read_csv(id_in, sep="\t", header=None, names=["FID", "IID"], dtype=str)
    n = len(ids)
    keep_mask = ids["IID"].isin(release_iids)
    keep_idx = [i for i, k in enumerate(keep_mask) if k]
    n_keep = len(keep_idx)
    print(f"  GRM: {n} → {n_keep} subjects")

    id_out.parent.mkdir(parents=True, exist_ok=True)
    ids_filtered = ids[keep_mask].reset_index(drop=True)
    ids_filtered.to_csv(id_out, sep="\t", header=False, index=False)

    dt = np.dtype("f4")
    grm_flat = np.fromfile(str(bin_in), dtype=dt)
    grm_full = np.zeros((n, n), dtype=dt)
    grm_full[np.tril_indices(n)] = grm_flat
    grm_full += grm_full.T
    np.fill_diagonal(grm_full, np.diag(grm_full) / 2)
    grm_sub = grm_full[np.ix_(keep_idx, keep_idx)]
    grm_sub[np.tril_indices(n_keep)].tofile(str(bin_out))
    print(f"  {bin_in.name}: {n}x{n} → {n_keep}x{n_keep}")

    if n_in.exists():
        n_flat = np.fromfile(str(n_in), dtype=dt)
        n_full = np.zeros((n, n), dtype=dt)
        n_full[np.tril_indices(n)] = n_flat
        n_full += n_full.T
        np.fill_diagonal(n_full, np.diag(n_full) / 2)
        n_sub = n_full[np.ix_(keep_idx, keep_idx)]
        n_sub[np.tril_indices(n_keep)].tofile(str(n_out))
        print(f"  {n_in.name}: {n}x{n} → {n_keep}x{n_keep}")


def filter_csv(in_path, out_path, id_cols, sep=","):
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(in_path, sep=sep, dtype=str)
    before = len(df)
    mask = df[id_cols].apply(lambda c: c.isin(release_iids)).any(axis=1)
    df = df[mask]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep=sep, index=False)
    print(f"  {in_path.name}: {before} → {len(df)} rows")


def filter_grm_text_gz(in_path, out_path):
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(
        in_path, sep=r"\s+", header=None,
        names=["FID1", "IID1", "FID2", "IID2", "N", "GRM"],
        dtype={"FID1": str, "IID1": str, "FID2": str, "IID2": str},
    )
    before = len(df)
    mask = df["IID1"].isin(release_iids) & df["IID2"].isin(release_iids)
    df = df[mask]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep=" ", header=False, index=False, compression="gzip")
    print(f"  {in_path.name}: {before} → {len(df)} rows")


# ===========================================================================
# PC-Relate GRM (binary)
# ===========================================================================
print("\n--- PC-Relate GRM (binary) ---")
filter_binary_grm(
    str(HANDOFF_DIR / "hbcd_pcrelate_grm"),
    str(RELEASE_DIR / "hbcd_pcrelate_grm"),
    id_ext=".grm.id", bin_ext=".grm.bin", n_ext=".grm.N.bin",
)

# ===========================================================================
# PC-Relate GRM (gzipped text)
# ===========================================================================
print("\n--- PC-Relate GRM (text gz) ---")
filter_grm_text_gz(
    HANDOFF_DIR / "hbcd_pcrelate_grm.gz",
    RELEASE_DIR / "hbcd_pcrelate_grm.gz",
)

# ===========================================================================
# PC-Relate GRM pairwise
# ===========================================================================
print("\n--- PC-Relate GRM pairwise ---")
filter_csv(
    HANDOFF_DIR / "hbcd_pcrelate_grm_pairwise.tsv",
    RELEASE_DIR / "hbcd_pcrelate_grm_pairwise.tsv",
    ["ID1", "ID2"], sep="\t",
)

# ===========================================================================
# PC-AiR 32 PCs clean
# ===========================================================================
print("\n--- PC-AiR 32 PCs clean ---")
filter_text(
    HANDOFF_DIR / "hbcd_pcair_32PCs_clean.tsv",
    "subject_id",
    RELEASE_DIR / "hbcd_pcair_32PCs_clean.tsv",
)

# ===========================================================================
# CNV slim clean (de-identified by 28-cnv-deid.py)
# ===========================================================================
print("\n--- CNV slim clean ---")
cnv_path = release_base / "CNV_slim_clean_deid.txt"
if not cnv_path.exists():
    print(f"  SKIP (not found — run 28-cnv-deid.py first): {cnv_path}")
else:
    cnv = pd.read_csv(cnv_path, sep="\t", dtype=str)
    before = len(cnv)
    cnv = cnv[cnv["sample_id"].isin(release_iids)]
    cnv_out = RELEASE_DIR / "CNV_slim_clean.txt"
    cnv_out.parent.mkdir(parents=True, exist_ok=True)
    cnv.to_csv(cnv_out, sep="\t", index=False)
    print(f"  CNV_slim_clean.txt: {before} → {len(cnv)} rows ({((before - len(cnv)) / before * 100):.1f}% removed)")

print("\nDone.")
