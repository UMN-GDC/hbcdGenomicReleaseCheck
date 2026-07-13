#!/usr/bin/env python3

# step 29: filter all de-identified handoff derivatives to release subjects
# Reads keep_list.txt (release IIDs) + de-IDed CNV, writes filtered copies
# to genotype_microarray subdirectories (genesis/, cnv/).
#
# Filtered outputs:
#   - PC-Relate GRM (binary)  →  genesis/pcrelate_relatedness.grm.{id,bin,N.bin}
#   - PC-Relate GRM (text gz) →  genesis/pcrelate_relatedness.grm.gz
#   - PC-Relate GRM pairwise  →  genesis/pcrelate_relatedness.tsv
#   - PC-AiR 32 PCs clean     →  genesis/pcair_weights.tsv
#   - CNV slim clean          →  cnv/CNV_slim.txt
#
# Usage: python 29-filter_release_outputs.py [--release-dir /path]
# Default release from HBCD_RELEASE env var or br_21p3.

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
genesis_dir = RELEASE_DIR / "genesis"
cnv_dir = RELEASE_DIR / "cnv"

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
    print(f"  Columns: {list(df.columns)}")
    mask = df[id_cols].apply(lambda c: c.isin(release_iids)).all(axis=1)
    df = df[mask]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep=sep, index=False)
    print(f"  {in_path.name}: {before} → {len(df)} rows")


def filter_grm_text_gz(in_path, orig_id_path, out_path):
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    if not orig_id_path.exists():
        print(f"  SKIP (not found): {orig_id_path}")
        return

    # Read original .grm.id to map 1-based indices → IIDs
    orig_ids = pd.read_csv(
        orig_id_path, sep="\t", header=None, names=["FID", "IID"], dtype=str
    )
    keep_indices = {i + 1 for i, row in orig_ids.iterrows()
                    if row["IID"] in release_iids}
    print(f"  GRM text: {len(orig_ids)} subjects → {len(keep_indices)} kept")

    # Text GRM format: IID1 IID2 N GRM  (1-based indices, not IID strings)
    df = pd.read_csv(
        in_path, sep=r"\s+", header=None,
        names=["IID1", "IID2", "N", "GRM"],
        dtype={"IID1": int, "IID2": int, "N": int, "GRM": float},
    )
    before = len(df)
    mask = df["IID1"].isin(keep_indices) & df["IID2"].isin(keep_indices)
    df = df[mask]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep=" ", header=False, index=False, compression="gzip")
    print(f"  {in_path.name}: {before} → {len(df)} rows")


# ===========================================================================
# PC-Relate GRM (binary)
# ===========================================================================
print("\n--- PC-Relate GRM (binary) → genesis/ ---")
filter_binary_grm(
    str(HANDOFF_DIR / "hbcd_pcrelate_grm"),
    str(genesis_dir / "pcrelate_relatedness"),
    id_ext=".grm.id", bin_ext=".grm.bin", n_ext=".grm.N.bin",
)

# ===========================================================================
# PC-Relate GRM (gzipped text)
# ===========================================================================
print("\n--- PC-Relate GRM (text gz) → genesis/ ---")
filter_grm_text_gz(
    HANDOFF_DIR / "hbcd_pcrelate_grm.grm.gz",
    HANDOFF_DIR / "hbcd_pcrelate_grm.grm.id",
    genesis_dir / "pcrelate_relatedness.grm.gz",
)

# ===========================================================================
# PC-Relate GRM pairwise
# ===========================================================================
print("\n--- PC-Relate GRM pairwise → genesis/ ---")
filter_csv(
    HANDOFF_DIR / "hbcd_pcrelate_grm_pairwise.tsv",
    genesis_dir / "pcrelate_relatedness.tsv",
    ["ID1", "ID2"], sep="\t",
)

# ===========================================================================
# PC-AiR 32 PCs clean
# ===========================================================================
print("\n--- PC-AiR 32 PCs clean → genesis/ ---")
filter_text(
    HANDOFF_DIR / "hbcd_pcair_32PCs_clean.tsv",
    "subject_id",
    genesis_dir / "pcair_weights.tsv",
)

# ===========================================================================
# CNV slim clean (de-identified by 28-cnv-deid.py)
# ===========================================================================
print("\n--- CNV slim clean → cnv/ ---")
cnv_path = cnv_dir / "CNV_slim_deid.txt"
cnv = pd.read_csv(cnv_path, sep="\t", dtype=str)
before = len(np.unique(cnv.sample_id))
cnv = cnv[cnv["sample_id"].isin(release_iids)]
cnv_out = cnv_dir / "CNV_slim.txt"
cnv_out.parent.mkdir(parents=True, exist_ok=True)
cnv.to_csv(cnv_out, sep="\t", index=False)
print(f"  CNV_slim.txt: {before} → {len(np.unique(cnv.sample_id))} unique Subjects")

# ===========================================================================
# CNV bookmark metrics (de-identified by 28-cnv-deid.py)
# ===========================================================================
print("\n--- CNV bookmark metrics → cnv/ ---")
bm_deid_path = cnv_dir / "HBCD_CNV_bookmark_metrics_clean_deid.csv"
bm = pd.read_csv(bm_deid_path, sep=",", dtype=str)
bm_before = len(np.unique(bm.sample_id))
bm = bm[bm["sample_id"].isin(release_iids)]
bm_out = cnv_dir / "HBCD_CNV_bookmark_metrics_clean.csv"
bm_out.parent.mkdir(parents=True, exist_ok=True)
bm.to_csv(bm_out, sep=",", index=False)
print(f"  HBCD_CNV_bookmark_metrics_clean.csv: {bm_before} → {len(np.unique(bm.sample_id))} unique Subjects")

# ===========================================================================
# Cross-form subject count consistency
# ===========================================================================
print("\n--- Cross-form subject count consistency ---")

def _unique_iids_from(path, id_col, sep="\t"):
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep=sep, dtype=str)
    return set(df[id_col].dropna().astype(str))

cnv_iids = _unique_iids_from(cnv_dir / "CNV_slim.txt", "sample_id", sep="\t")
bm_iids = _unique_iids_from(cnv_dir / "HBCD_CNV_bookmark_metrics_clean.csv", "sample_id", sep=",")

rc_cnv = set(i[:-1] for i in cnv_iids)
print(f"  CNV_slim.txt unique subjects        : {len(rc_cnv):>6}")
rc_bm = set(i[:-1] for i in bm_iids)
print(f"  Bookmarks unique subjects           : {len(rc_bm):>6}")

common = rc_cnv & rc_bm
only_cnv = rc_cnv - rc_bm
only_bm = rc_bm - rc_cnv
print(f"  Subjects in BOTH CNV+Bookmarks     : {len(common):>6}")
print(f"  Subjects ONLY in CNV               : {len(only_cnv):>6}")
print(f"  Subjects ONLY in Bookmarks          : {len(only_bm):>6}")

print(f"\nReference release_iids (keep_list.txt): {len(release_iids):>6}")
missing_from_release = cnv_iids - release_iids
if missing_from_release:
    print(f"  WARN: {len(missing_from_release)} CNV IID(s) not in release set")
missing_from_release = bm_iids - release_iids
if missing_from_release:
    print(f"  WARN: {len(missing_from_release)} bookmark IID(s) not in release set")

# Compare to fam (release GDA set)
print("\nDone.")
