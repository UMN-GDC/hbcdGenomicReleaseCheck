#!/usr/bin/env python3

# step 28: filter all derivative outputs to only include release subjects
# Reads keep_list.txt (release IIDs) and subsets:
#   - GCTA binary GRM (.grm.id, .grm.bin, .grm.N.bin)
#   - PLINK GRM (.rel, .rel.id)
#   - PC-AiR scores, unrelated/related IDs
#   - PC-Relate CSV exports (pairs, kinmat)

import pandas as pd
import numpy as np
from pathlib import Path

RELEASE_DIR = Path(
    "/projects/standard/basu_hbcd/shared/HBCD_genomics_release_br_21p2/data/"
)
WORK = Path("/projects/standard/basu_hbcd/shared/HST_HBCD_Transfer_May2026")
PCA_DIR = WORK / "PCA"
PCA_IR = PCA_DIR / "pca_ir"
NAME = "hbcd_rsid_harmonized"

# -- load release IIDs --
keep = pd.read_csv(
    RELEASE_DIR / ".." / "keep_list.txt",
    sep=r"\s+",
    header=None,
    names=["FID", "IID"],
)
release_iids = set(keep["IID"].astype(str))
print(f"Release IIDs: {len(release_iids)}")


def filter_text_by_id(path, id_col, out_path=None):
    """Filter a text file keeping only rows where id_col is in release_iids."""
    if not path.exists():
        print(f"  SKIP (not found): {path}")
        return
    df = pd.read_csv(path, sep="\t", dtype=str)
    before = len(df)
    df = df[df[id_col].isin(release_iids)]
    n = len(df)
    out = out_path or path
    df.to_csv(out, sep="\t", index=False)
    print(f"  {path.name}: {before} → {n} rows ({before - n} removed)")


def filter_matrix_by_ids(rel_path, id_path):
    """Filter a square matrix (.rel) to only keep release subjects."""
    if not rel_path.exists() or not id_path.exists():
        print(f"  SKIP (not found): {rel_path.name} / {id_path.name}")
        return
    ids = pd.read_csv(id_path, sep=r"\s+", header=None, names=["FID", "IID"])
    orig_n = len(ids)
    keep_mask = ids["IID"].astype(str).isin(release_iids)
    keep_idx = [i for i, k in enumerate(keep_mask) if k]
    ids_filtered = ids[keep_mask].reset_index(drop=True)
    ids_filtered.to_csv(id_path, sep=" ", header=False, index=False)

    mat = np.loadtxt(str(rel_path))
    mat_filtered = mat[np.ix_(keep_idx, keep_idx)]
    np.savetxt(str(rel_path), mat_filtered, fmt="%.10g")
    print(f"  {rel_path.name}: {orig_n}x{orig_n} → {len(ids_filtered)}x{len(ids_filtered)}")


def filter_binary_grm(prefix):
    """Filter GCTA binary GRM (.grm.id, .grm.bin, .grm.N.bin) to release subjects."""
    id_file = Path(f"{prefix}.grm.id")
    bin_file = Path(f"{prefix}.grm.bin")
    n_file = Path(f"{prefix}.grm.N.bin")
    if not id_file.exists() or not bin_file.exists():
        print(f"  SKIP (not found): {prefix}.grm.*")
        return

    ids = pd.read_csv(id_file, sep="\t", header=None, names=["FID", "IID"], dtype=str)
    n = len(ids)
    keep_mask = ids["IID"].isin(release_iids)
    keep_idx = [i for i, k in enumerate(keep_mask) if k]
    n_keep = len(keep_idx)
    print(f"  {id_file.name}: {n} → {n_keep} subjects")

    # Read lower triangle → full symmetric → subset → lower triangle
    dt = np.dtype("f4")
    grm_flat = np.fromfile(str(bin_file), dtype=dt)
    grm_full = np.zeros((n, n), dtype=dt)
    grm_full[np.tril_indices(n)] = grm_flat
    grm_full = grm_full + grm_full.T - np.diag(np.diag(grm_full))
    grm_sub = grm_full[np.ix_(keep_idx, keep_idx)]
    grm_sub_flat = grm_sub[np.tril_indices(n_keep)]
    grm_sub_flat.tofile(str(bin_file))
    print(f"  {bin_file.name}: {n}x{n} → {n_keep}x{n_keep}")

    # Same for N.bin
    if n_file.exists():
        n_flat = np.fromfile(str(n_file), dtype=dt)
        n_full = np.zeros((n, n), dtype=dt)
        n_full[np.tril_indices(n)] = n_flat
        n_full = n_full + n_full.T - np.diag(np.diag(n_full))
        n_sub = n_full[np.ix_(keep_idx, keep_idx)]
        n_sub_flat = n_sub[np.tril_indices(n_keep)]
        n_sub_flat.tofile(str(n_file))
        print(f"  {n_file.name}: {n}x{n} → {n_keep}x{n_keep}")

    # Write filtered .grm.id
    ids_filtered = ids[keep_mask].reset_index(drop=True)
    ids_filtered.to_csv(id_file, sep="\t", header=False, index=False)


def filter_csv_by_ids(csv_path, id_cols):
    """Filter a CSV file keeping rows where any id_col matches a release IID."""
    if not csv_path.exists():
        print(f"  SKIP (not found): {csv_path}")
        return
    df = pd.read_csv(csv_path, dtype=str)
    before = len(df)
    mask = df[id_cols].apply(lambda c: c.isin(release_iids)).any(axis=1)
    df = df[mask]
    df.to_csv(csv_path, index=False)
    print(f"  {csv_path.name}: {before} → {len(df)} rows")


# ===========================================================================
# GCTA binary GRM
# ===========================================================================
print("\n--- GCTA binary GRM ---")
filter_binary_grm(str(PCA_DIR / "hbcd_gcta_grm"))

# ===========================================================================
# PLINK GRM
# ===========================================================================
print("\n--- PLINK GRM ---")
filter_matrix_by_ids(
    PCA_DIR / "hbcd_plink_grm.rel",
    PCA_DIR / "hbcd_plink_grm.rel.id",
)

# ===========================================================================
# PC-AiR text exports
# ===========================================================================
print("\n--- PC-AiR ---")
filter_text_by_id(PCA_IR / f"{NAME}_pc_scores.txt", "sample.id")
filter_text_by_id(PCA_IR / f"{NAME}_unrelated_ids.txt", "SampleID")
filter_text_by_id(PCA_IR / f"{NAME}_related_ids.txt", "SampleID")

# ===========================================================================
# PC-Relate CSV/TSV exports
# ===========================================================================
print("\n--- PC-Relate pairs ---")
for ext in ["csv", "tsv"]:
    filter_csv_by_ids(PCA_IR / f"{NAME}_pcrelate_pairs.{ext}", ["ID1", "ID2"])
    filter_csv_by_ids(PCA_IR / f"{NAME}_pcrelate_self.{ext}", ["ID"])

print("\n--- PC-Relate kinmat long ---")
for ext in ["csv", "tsv"]:
    filter_csv_by_ids(PCA_IR / f"{NAME}_pcrelate_kinmat_long.{ext}", ["ID1", "ID2"])

print("\n--- PC-Relate kinmat wide ---")
for ext in ["csv", "tsv"]:
    p = PCA_IR / f"{NAME}_pcrelate_kinmat_wide.{ext}"
    if p.exists():
        df = pd.read_csv(p, dtype=str)
        cols = ["SampleID"] + [c for c in df.columns[1:] if c in release_iids]
        df = df[df["SampleID"].isin(release_iids)][cols]
        df.to_csv(p, index=False)
        print(f"  {p.name}: {len(df)} samples x {len(cols)-1} cols")

# ===========================================================================
# PC-AiR eigenvalues (no IDs to filter, just copy)
# ===========================================================================
print("\nDone.")
