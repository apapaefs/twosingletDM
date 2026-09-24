"""Independent, nullable Upsilon lepton-bound assessment for TRSM scans.

The collaborator's prompt, narrow-resonance prescription is kept verbatim.
No limits or low-mass decay inputs are extrapolated.
"""

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from flavour import flavour_observables as ups
from trsm_inputs import M1, json_safe, nullable_and

FLAVOUR_METHOD = "upsilon_leptons_v1"
FLAVOUR_COLUMNS = (
    "flavour", "flavour_method", "flavour_status", "flavour_reason",
    "flavour_max_ratio", "flavour_details",
)


def flavour_configuration():
    """Validate the actual input files and fingerprint the complete prescription."""
    ups._load_limit_curve.cache_clear()
    coverage = {}
    for filename in ups.LIMIT_FILES:
        path = ups.LIMIT_DATA_DIR / filename
        raw = np.genfromtxt(path, delimiter=",", skip_header=1, comments="#")
        if (raw.ndim != 2 or raw.shape[0] < 2 or raw.shape[1] < 2
                or not np.all(np.isfinite(raw[:, :2])) or np.any(raw[:, 1] <= 0)
                or np.any(raw[:, 0] <= 0)
                or len(np.unique(raw[:, 0])) != len(raw)):
            raise ValueError(f"Invalid flavour limit data: {path}")
        masses, _ = ups._load_limit_curve(filename)
        coverage[filename] = [float(masses[0]), float(masses[-1])]
    root = Path(__file__).resolve().parent
    paths = [Path(__file__), root / "generate_trsm_info.py", root / "trsm_inputs.py",
             root / "config/sm-inputs-v2.json", root / "YR/higgsBR_HiggsTools_YR4_lowmass.txt"]
    paths += [ups.LIMIT_DATA_DIR / name for name in (*ups.LIMIT_FILES, "flavour_inputs.py", "flavour_observables.py", "__init__.py")]
    return {
        "method": FLAVOUR_METHOD, "enabled": True,
        "source": "https://arxiv.org/pdf/2112.11852",
        "prescription": "strongest individual Belle/BaBar bound per channel; prompt narrow resonance",
        "search_ranges_GeV": {"mumu": [ups.MUMU_MIN, ups.SEARCH_MAX],
                              "tautau": [ups.TAUTAU_MIN, ups.SEARCH_MAX]},
        "decay_min_mass_GeV": 4.0, "curve_coverage_GeV": coverage,
        "sources_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in paths},
    }


def assess_scalar(name, mass, coupling, base_brs=None, total_width=None, *, stable=False):
    """Assess one scalar using pre-exotic BR arrays and its physical total width.

    Width-independent zero production and absent search coverage are non-vetoes.
    Missing/invalid decay inputs within coverage are unavailable, never passes.
    """
    result = {"name": name, "mass_GeV": None, "mixing_squared": None,
              "base_width_GeV": None, "total_width_GeV": None,
              "passed": None, "status": "unassessed", "reason": "", "channels": {}}
    try:
        mass = float(mass)
        if not math.isfinite(mass) or mass <= 0:
            raise ValueError("scalar mass must be finite and positive")
        result["mass_GeV"] = mass
        if mass < ups.MUMU_MIN or mass > ups.SEARCH_MAX:
            result.update(passed=True, status="outside_coverage", reason="outside both search ranges")
            return result
        coupling = complex(coupling)
        if not (math.isfinite(coupling.real) and math.isfinite(coupling.imag)):
            raise ValueError("doublet projection must be finite")
        result["mixing_squared"] = abs(coupling)**2
        result["coupling_real"] = coupling.real
        result["coupling_imag"] = coupling.imag
        if stable and coupling != 0:
            raise ValueError("stable scalar has a nonzero doublet projection")
        if coupling == 0:
            detail = ups.EvaluateUpsLeptonBounds(mass, coupling, 0.0, 0.0)
            result.update(detail, status="zero_signal", reason="stable scalar" if stable else "zero doublet projection")
            return result
        base_width = float(base_brs[-1])
        width = float(total_width)
        if not math.isfinite(base_width) or base_width <= 0:
            raise ValueError("base decay width must be finite and positive for nonzero mixing")
        if not math.isfinite(width) or width <= 0 or width < base_width:
            raise ValueError("physical total width must be finite, positive, and at least the base width")
        brs = [float(base_brs[index]) for index in (2, 1)]
        if any(not math.isfinite(br) or not 0 <= br <= 1 for br in brs):
            raise ValueError("base lepton branching fractions must lie in [0, 1]")
        result.update(base_width_GeV=base_width, total_width_GeV=width)
        physical = [br * (base_width / width) for br in brs]
        detail = ups.EvaluateUpsLeptonBounds(mass, coupling, *physical)
        covered = [c for c in detail["channels"].values() if c["covered"]]
        status = ("excluded" if not detail["passed"] else "outside_coverage" if not covered
                  else "zero_signal" if all(c["prediction"] == 0 for c in covered) else "passed")
        reason = ", ".join(k for k, c in detail["channels"].items() if c["excluded"])
        result.update(detail, status=status, reason=reason)
    except (ValueError, TypeError, IndexError, KeyError, ArithmeticError, OSError) as error:
        result.update(passed=None, status="unassessed", reason=str(error))
    return result


def combine_flavour(scalars):
    details = {s["name"]: s for s in scalars}
    passed = nullable_and(s["passed"] for s in details.values())
    statuses = {s["status"] for s in details.values()}
    status = ("excluded" if passed is False else "unassessed" if passed is None
              else "passed" if "passed" in statuses else "zero_signal" if "zero_signal" in statuses
              else "outside_coverage")
    ratios = [c["ratio"] for s in details.values() for c in s["channels"].values()
              if c["ratio"] is not None]
    reasons = [f"{name}: {s['reason']}" for name, s in details.items()
               if s["reason"] and (s["passed"] is not True or status != "passed")]
    return {"flavour": passed, "flavour_method": FLAVOUR_METHOD,
            "flavour_status": status, "flavour_reason": "; ".join(reasons),
            "flavour_max_ratio": max(ratios) if ratios else None,
            "flavour_details": json.dumps(json_safe(details), separators=(",", ":"), allow_nan=False)}


def unassessed_flavour(reason):
    return {"flavour": None, "flavour_method": FLAVOUR_METHOD,
            "flavour_status": "unassessed", "flavour_reason": str(reason),
            "flavour_max_ratio": None, "flavour_details": "{}"}


def generated_flavour_updates(m2, m3, couplings, brs, widths, *, vx=0):
    return combine_flavour([
        assess_scalar(f"h{i}", mass, coupling, base, width, stable=(i == 3 and vx == 0))
        for i, (mass, coupling, base, width) in enumerate(
            zip((M1, m2, m3), couplings, brs, widths), start=1)
    ])


def flavour_updates_from_row(row):
    """Reconstruct only the decay inputs relevant to a saved vx=0 scan.

    In the Upsilon mass window h2 -> h1 h1 is closed. The existing SM decay
    and canonical portal helpers therefore give the complete h2 width without
    invoking any external solver or trusting historical cached widths/BRs.
    """
    try:
        def number(key):
            value = float(row[key])
            if not math.isfinite(value):
                raise ValueError(f"{key} is not finite")
            return value
        m2, m3, vx = (number(k) for k in ("M2", "M3", "vx"))
        if vx != 0:
            raise ValueError("flavour-only reevaluation supports vx=0 scans")
        k2, base, width = None, None, None
        if ups.MUMU_MIN <= m2 <= ups.SEARCH_MAX:
            from generate_trsm_info import (Rmatrix, BR_interpolators_SM, calc_h2_BRs,
                                           scalar_to_identical_scalar_width, vxzero_portal_couplings)
            angle = number("a12")
            k2 = Rmatrix(angle, 0.0, 0.0)[1][0]
            if k2 != 0:
                # K112 cannot contribute below 9.2 GeV.
                base = calc_h2_BRs(BR_interpolators_SM, M1, m2, k2, 0.0)
                _, k233 = vxzero_portal_couplings(number("lPhiX"), number("lSX"), number("vs"), angle)
                width = float(base[-1]) + scalar_to_identical_scalar_width(m3, m2, k233)
        return generated_flavour_updates(m2, m3, (None, k2, 0.0), (None, base, None),
                                         (None, width, 0.0), vx=vx)
    except (ValueError, TypeError, KeyError, ArithmeticError, OSError) as error:
        return unassessed_flavour(error)
