"""Dataset locations shared by the physics provider and provenance checks."""

import os
from pathlib import Path


def higgs_dataset_path(name):
    variables = {"hbdataset": "TRSM_HB_DATASET", "hsdataset": "TRSM_HS_DATASET"}
    variable = variables[name]
    fallback = Path(__file__).resolve().parent.parent / name
    return Path(os.environ.get(variable, fallback)).expanduser().resolve()
