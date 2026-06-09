import pandas as pd

fam = pd.read_csv("../HST_HBCD_Transfer_July2025/HBCD/HBCD.fam", sep = "\s+", header= None)
fam.columns = ["bad_fid", "bad_iid", "PAT", "MAT", "Sex", "Pheno"]
fam["iid"] = fam['bad_iid'].str[21:]
fam["pscid"] = fam['iid'].str[:9]
fam["rel"] = fam['iid'].str[9]

id = pd.read_csv("../data/release_identifiers_20260526.csv")

fam = fam.merge(
    id[["pscid", "release_candid"]].drop_duplicates(subset="pscid"),
    on="pscid",
    how="left",
)

fam["release_candidMC"] = fam["release_candid"] + fam["rel"]

fam[["release_candid", "release_candidMC", "PAT","MAT","Sex","Pheno"]].to_csv("fullHBCD.fam", sep = " ", index=False, na_rep = "NA")
