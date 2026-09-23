"""Versioned pole/EW inputs and explicitly separate RG boundary convention."""

import math
import json
from pathlib import Path

_INPUTS = json.loads((Path(__file__).resolve().parent / "config/sm-inputs-v2.json").read_text())

PHYSICS_VERSION = "trsm_constraints_v2"
SCHEMA_VERSION = 2
M1 = _INPUTS["M1_GeV"]
GF = _INPUTS["GF_GeV^-2"]
MW = _INPUTS["MW_GeV"]
MZ = _INPUTS["MZ_GeV"]
VEV = math.sqrt(1.0 / (math.sqrt(2.0) * GF))
SW = math.sqrt(1.0 - (MW / MZ) ** 2)
EE = 2.0 * MW * SW / VEV
RELIC_UPPER_LIMIT = 0.121
ABUNDANCE_REFERENCE = 0.12

# Retained one-loop screening prescription, not precision MS-bar matching.
RG_BOUNDARY = {
    "scale_GeV": 91.0, "maximum_scale_GeV": 1000.0,
    "g3": 1.221, "g2": math.sqrt(0.424),
    "g1": math.sqrt(0.1273), "yt": 0.96738,
    "scalar_matching": "tree masses and mixing at the boundary scale",
    "sm_matching": "inherited fixed running couplings; no precision matching",
    "threshold_matching": False,
}


def sm_inputs():
    return {"scheme": "GF_MW_MZ_tree_v1", "M1_GeV": M1, "GF_GeV^-2": GF,
            "MW_GeV": MW, "MZ_GeV": MZ, "v_GeV": VEV, "SW": SW, "EE": EE}


def micromegas_sm_inputs():
    # Quark running-mass inputs retain their backend scheme; they are echoed
    # by the driver instead of equating them to BSMPT's pole masses.
    return {"Mh": M1, "MW": MW, "SW": SW, "EE": EE, "Mtp": _INPUTS["top_pole_GeV"]}


def nullable_and(values):
    values = tuple(values)
    if any(value is False for value in values):
        return False
    return True if all(value is True for value in values) else None


def json_safe(value):
    """Represent unavailable numbers as JSON null, never nonstandard NaN/Infinity."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value
