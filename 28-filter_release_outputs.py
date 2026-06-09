#!/usr/bin/env python3

# step 28: filter all derivative outputs to only include release subjects
# Reads keep_list.txt (release IIDs) and writes filtered copies to the
# release data directory, preserving the full-data originals in DATA_DIR.
#
# Filtered outputs:
#   - GCTA binary GRM  →  release_dir/hbcd_gcta_grm.{grm.id,grm.bin,grm.N.bin}
#   - PLINK GRM        →  release_dir/hbcd_plink_grm.{rel,rel.id}
#   - PC-AiR scores    →  release_dir/hbcd_rsid_harmonized_pc_scores.txt
#   - PC-Relate        →  release_dir/hbcd_rsid_harmonized_pcrelate_pairs.{csv,tsv}
#   - etc.
#
# Usage: python 28-filter_release_outputs.py [--release-dir /path/to/release/data]
# Default release from HBCD_RELEASE env var or br_21p2.

import sys
import shutil
import pandas as pd
import numpy as np
from pathlib import Path
from _lib import DATA_DIR, get_release_dir

NAME = "hbcd_rsid_harmonized"

# Accept --release-dir / -r override, default from env/HBCD_RELEASE
if len(sys.argv) > 1 and sys.argv[1] in ("--release-dir", "-r"):
    RELEASE_DIR = Path(sys.argv[2])
else:
    RELEASE_DIR = get_release_dir()

release_base = RELEASE_DIR.parent

# -- load release IIDs --
keep = pd.read_csv(
    release_base / "keep_list.txt",
    sep=r"\s+",
    header=None,
    names=["FID", "IID"],
)
release_iids = set(keep["IID"].astype(str))
print(f"Release IIDs: {len(release_iids)}")


# ── filter helpers (always write to out_path; leave source untouched) ────


def filter_text(in_path, id_col, out_path):
    """Read a delimited text file, keep rows where *id_col* is in
    ``release_iids``, write filtered result to *out_path*."""
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    before = len(df)
    df = df[df[id_col].isin(release_iids)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  {in_path.name}: {before} → {len(df)} rows")


def filter_matrix(in_rel, in_id, out_rel, out_id):
    """Filter a square matrix (.rel + .rel.id) to only release subjects."""
    if not in_rel.exists() or not in_id.exists():
        print(f"  SKIP (not found): {in_rel.name} / {in_id.name}")
        return
    ids = pd.read_csv(in_id, sep=r"\s+", header=None, names=["FID", "IID"])
    orig_n = len(ids)
    keep_mask = ids["IID"].astype(str).isin(release_iids)
    keep_idx = [i for i, k in enumerate(keep_mask) if k]
    ids_filtered = ids[keep_mask].reset_index(drop=True)
    out_id.parent.mkdir(parents=True, exist_ok=True)
    ids_filtered.to_csv(out_id, sep=" ", header=False, index=False)

    mat = np.loadtxt(str(in_rel))
    mat_filtered = mat[np.ix_(keep_idx, keep_idx)]
    np.savetxt(str(out_rel), mat_filtered, fmt="%.10g")
    print(f"  {in_rel.name}: {orig_n}x{orig_n} → {len(ids_filtered)}x{len(ids_filtered)}")


def filter_binary_grm(in_prefix, out_prefix):
    """Filter GCTA binary GRM (.grm.id, .grm.bin, .grm.N.bin) to release subjects."""
    id_in = Path(f"{in_prefix}.grm.id")
    bin_in = Path(f"{in_prefix}.grm.bin")
    n_in = Path(f"{in_prefix}.grm.N.bin")
    id_out = Path(f"{out_prefix}.grm.id")
    bin_out = Path(f"{out_prefix}.grm.bin")
    n_out = Path(f"{out_prefix}.grm.N.bin")

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


def filter_csv(in_path, out_path, id_cols):
    """Filter a CSV keeping rows where any *id_cols* matches a release IID."""
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(in_path, dtype=str)
    before = len(df)
    mask = df[id_cols].apply(lambda c: c.isin(release_iids)).any(axis=1)
    df = df[mask]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"  {in_path.name}: {before} → {len(df)} rows")


def filter_kinmat_wide(in_path, out_path):
    """Filter a wide-format kinship matrix CSV keeping only release columns/rows."""
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(in_path, dtype=str)
    cols = ["SampleID"] + [c for c in df.columns[1:] if c in release_iids]
    df = df[df["SampleID"].isin(release_iids)][cols]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"  {in_path.name}: {len(df)} samples x {len(cols)-1} cols")


# ===========================================================================
# GCTA binary GRM
# ===========================================================================
print("\n--- GCTA binary GRM ---")
filter_binary_grm(
    str(DATA_DIR / "hbcd_gcta_grm"),
    str(RELEASE_DIR / "hbcd_gcta_grm"),
)

# ===========================================================================
# PLINK GRM
# ===========================================================================
print("\n--- PLINK GRM ---")
filter_matrix(
    DATA_DIR / "hbcd_plink_grm.rel",
    DATA_DIR / "hbcd_plink_grm.rel.id",
    RELEASE_DIR / "hbcd_plink_grm.rel",
    RELEASE_DIR / "hbcd_plink_grm.rel.id",
)

# ===========================================================================
# PC-AiR text exports
# ===========================================================================
print("\n--- PC-AiR ---")
filter_text(
    DATA_DIR / f"{NAME}_pc_scores.txt", "sample.id",
    RELEASE_DIR / f"{NAME}_pc_scores.txt",
)
filter_text(
    DATA_DIR / f"{NAME}_unrelated_ids.txt", "SampleID",
    RELEASE_DIR / f"{NAME}_unrelated_ids.txt",
)
filter_text(
    DATA_DIR / f"{NAME}_related_ids.txt", "SampleID",
    RELEASE_DIR / f"{NAME}_related_ids.txt",
)

# ===========================================================================
# PC-AiR eigenvalues (no sample IDs — just copy)
# ===========================================================================
print("\n--- PC-AiR eigenvalues ---")
for fname in [f"{NAME}_pcair_eigenvalues.csv", f"{NAME}_pcair_eigenvalues.tsv"]:
    src = DATA_DIR / fname
    if src.exists():
        shutil.copy2(src, RELEASE_DIR / fname)
        print(f"  Copied: {fname}")

# ===========================================================================
# PC-AiR scores (CSV/TSV from export_pcrelate.R)
# ===========================================================================
print("\n--- PC-AiR scores (export) ---")
for ext in ["csv", "tsv"]:
    filter_csv(
        DATA_DIR / f"{NAME}_pcair_scores.{ext}",
        RELEASE_DIR / f"{NAME}_pcair_scores.{ext}",
        ["SampleID"],
    )

# ===========================================================================
# PC-AiR unrelated/related IDs (export format)
# ===========================================================================
print("\n--- PC-AiR unrelated/related (export) ---")
for fstem in ["pcair_unrelated_ids", "pcair_related_ids"]:
    for ext in ["csv", "tsv"]:
        filter_csv(
            DATA_DIR / f"{NAME}_{fstem}.{ext}",
            RELEASE_DIR / f"{NAME}_{fstem}.{ext}",
            ["SampleID"],
        )

# ===========================================================================
# PC-Relate pairs + IBD + self
# ===========================================================================
print("\n--- PC-Relate pairs ---")
for ext in ["csv", "tsv"]:
    filter_csv(
        DATA_DIR / f"{NAME}_pcrelate_pairs.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_pairs.{ext}",
        ["ID1", "ID2"],
    )
    filter_csv(
        DATA_DIR / f"{NAME}_pcrelate_ibd.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_ibd.{ext}",
        ["ID1", "ID2"],
    )
    filter_csv(
        DATA_DIR / f"{NAME}_pcrelate_self.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_self.{ext}",
        ["ID"],
    )

# ===========================================================================
# PC-Relate kinship matrices
# ===========================================================================
print("\n--- PC-Relate kinmat long ---")
for ext in ["csv", "tsv"]:
    filter_csv(
        DATA_DIR / f"{NAME}_pcrelate_kinmat_long.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_kinmat_long.{ext}",
        ["ID1", "ID2"],
    )

print("\n--- PC-Relate kinmat wide ---")
for ext in ["csv", "tsv"]:
    filter_kinmat_wide(
        DATA_DIR / f"{NAME}_pcrelate_kinmat_wide.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_kinmat_wide.{ext}",
    )

print("\nDone.")
