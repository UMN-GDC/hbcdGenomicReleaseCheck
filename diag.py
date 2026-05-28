"""Diagnostic: trace the filter chain step by step."""
import warnings
warnings.filterwarnings("ignore")

from _lib import *
import pandas as pd

par = load_par_visit_candids()
print("par_visit candids:", len(par))
print("  sample:", sorted(par)[:5])

idi = load_identifiers()
print("identifiers rows:", len(idi))
d = idi["candid"].dtype
print("  candid dtype:", d)
uniq = sorted(idi["candid"].dropna().unique())[:5]
print("  candid sample:", uniq)

idi_par = idi[idi["candid"].isin(par)]
print("idi rows with candid in par_visit:", len(idi_par))
valid_rc = set(int(v) for v in idi_par["release_candid"].dropna().unique())
print("  unique release_candids:", len(valid_rc))
print("  sample:", sorted(valid_rc)[:5])
exc = load_excluded_release_candids()
valid_rc -= exc
print("  after removing excluded:", len(valid_rc))

batch = pd.read_csv(DATA_DIR / "batch.info", sep=r"\s+")
print("batch.info rows:", len(batch))
batch["rc"] = pd.to_numeric(batch["IID"].str[:-1])
batch_rc = set(int(v) for v in batch["rc"].dropna().unique())
print("  unique release_candids:", len(batch_rc))
print("  batch_rc & valid_rc:", len(batch_rc & valid_rc))

fam = pd.read_csv(
    str(DATA_DIR / "HBCD.fam"),
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
)
fam["pscid"] = fam["IID"].str[-10:-1]
fam["_orig_rel"] = fam["IID"].str[-1]
print("raw .fam subjects:", len(fam))

m = (
    fam.merge(idi, how="left", on="pscid")
    .merge(
        batch.drop(columns=["IID"]).rename(columns={"rc": "release_candid"}),
        how="left",
        on="release_candid",
    )
)
m["rc_valid"] = m["release_candid"].isin(valid_rc)
print("  subjects with valid_rc:", m["rc_valid"].sum())
print("  + non-NaN visit:", m[m["rc_valid"]]["visit"].notna().sum())
print("  + non-NaN plate:", m[m["rc_valid"]]["plate_number"].notna().sum())
