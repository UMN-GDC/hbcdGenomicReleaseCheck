import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from _lib import get_release_dir, get_release_base


@pytest.fixture(scope="session")
def release_dir():
    return get_release_dir()


@pytest.fixture(scope="session")
def release_base():
    return get_release_base()


@pytest.fixture(scope="session")
def release_iids():
    """Load the set of release IIDs from keep_list.txt."""
    import pandas as pd
    keep = pd.read_csv(
        get_release_base() / "keep_list.txt",
        sep=r"\s+", header=None, names=["FID", "IID"],
    )
    return set(keep["IID"].astype(str))
