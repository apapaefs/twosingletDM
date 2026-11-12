"""Astrophysical bounds on H3 self-interactions for vx=0 TRSM points.

This module is separate from the cross-section code.  ``trsm_sidm.py``
computes sigma_V(v); this module compares it with a table of observational
upper limits on sigma/m (``datafiles/sidm_limits_v1.json``) and writes, for
each point, the model value and the ratio to every selected limit.

How a limit is evaluated
------------------------
Each table entry says how the model cross section is to be reduced to one
number (see the ``conventions`` block of the table):

* ``maxwell``: ``<sigma_V v^n>/<v^n>`` over a Maxwellian whose mean relative
  velocity is the one quoted by the paper, with n = 1 by default (the
  ``<sigma v>/<v>`` definition used for the group and cluster limits);
* ``fixed``: sigma_V(v) over a velocity range, combined with min/max/mean.
  Papers that quote no collision velocity are evaluated over a range and the
  smallest value is used.

Sub-dominant H3
---------------
If H3 is only a fraction ``xi = min(1, Omega h^2 / 0.12)`` of the dark
matter, each model value is multiplied by ``xi**p`` before the comparison,
with p = 2 by default.  The limits constrain the total dark-matter
distribution.  Only a fraction xi of the mass can scatter, and each H3
particle sees a target density reduced by xi, so the number of scatterings
per unit total mass is ``xi^2 rho sigma/m``.  This holds while each H3
particle scatters at most about once.  For larger cross sections the effect
saturates, because at most a fraction xi of the mass is affected.  That is
why no bound is applied when ``xi < min_fraction`` (default 0.1): a
component of 10% or less can have arbitrarily strong self-interactions
(Pollack et al., arXiv:1501.00017; Fan et al., arXiv:1303.1521).  p = 1
rescales only the scattering probability of an H3 particle and overstates
the constraint.  ``--no-rescale`` treats H3 as all of the dark matter.
None of the cited papers gives a rescaling law for cluster limits, so p and
the cut are modelling choices.

Verdict (column ``sidm_excluded``)
----------------------------------
``True``   xi >= min_fraction, the cross section is reliable, and a rescaled
           ratio exceeds 1.
``False``  no selected limit is exceeded, or xi < min_fraction.
``nan``    unavailable: the calculation failed, the abundance is missing, or
           the point is in the ``strong`` regime of trsm_sidm (cross section
           not under control) while the bounds would apply.  Unavailable is
           never to be read as passing.

A ``strong`` point is still given ``False`` when unitarity alone keeps it
below every selected limit: ``sidm_bounds_unitarity_ratio`` is the largest
rescaled ratio that is possible if all even partial waves up to l = 8
saturate ``sigma_l <= 8 pi (2l+1)/k^2`` (and ``sigma_V <= 3/2 sigma``).
This needs no dynamics.  It is used only if the typical momentum is below
the lightest mediator mass, so that partial waves beyond l = 8 are negligible.

Usage::

    python trsm_sidm_bounds.py list --set all
    python trsm_sidm_bounds.py scan scan.dat --output scan_sidm_bounds.dat
    python trsm_sidm_bounds.py point --m2 10.00005 --m3 5 --vs 500 --a12 0 \
        --lx 0.1 --lphix 0 --lsx 0.008 --omega 0.12
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import trsm_sidm as sidm


BOUNDS_METHOD = "trsm_sidm_bounds_v1"
DEFAULT_LIMITS_FILE = Path(__file__).resolve().parent / "datafiles" / "sidm_limits_v1.json"
DEFAULT_SET = "default"
DEFAULT_OMEGA_COLUMN = "dm_omega"
DEFAULT_OMEGA_REFERENCE = 0.12      # trsm_inputs.ABUNDANCE_REFERENCE
DEFAULT_ABUNDANCE_POWER = 2.0
DEFAULT_MIN_FRACTION = 0.1
FIXED_GRID_POINTS = 13
COMBINERS = {"min": np.min, "max": np.max, "mean": np.mean}
# sigma_V <= 3/2 sigma and sigma <= sum over even l <= 8 of 8 pi (2l+1)/k^2.
UNITARITY_WEIGHT = 1.5 * sum(2 * l + 1 for l in range(0, sidm.DEFAULT_LMAX + 1, 2))


@dataclass(frozen=True)
class Limit:
    """One observational upper limit on sigma/m."""

    id: str
    label: str
    limit_cm2_g: float
    confidence_level: str
    kind: str                      # "maxwell" or "fixed"
    reference: str
    sets: tuple
    mean_relative_velocity_kms: float = math.nan
    moment: int = 1
    velocity_range_kms: tuple = ()
    combine: str = "min"
    velocity_assumed: bool = False
    notes: str = ""

    @property
    def dispersion_kms(self):
        return sidm.mean_relative_velocity_to_dispersion(self.mean_relative_velocity_kms)

    @property
    def velocities_kms(self):
        lo, hi = self.velocity_range_kms
        return np.array([lo]) if hi == lo else np.geomspace(lo, hi, FIXED_GRID_POINTS)

    def describe_velocity(self):
        if self.kind == "maxwell":
            return f"<v> = {self.mean_relative_velocity_kms:g} km/s (n = {self.moment})"
        lo, hi = self.velocity_range_kms
        return f"{self.combine} over v = {lo:g}-{hi:g} km/s"


def load_limits(path=None, selected=DEFAULT_SET, moment_override=None):
    """Return (limits, table) for the named set of the limit table."""
    path = Path(path) if path else DEFAULT_LIMITS_FILE
    table = json.loads(path.read_text(encoding="utf-8"))
    if table.get("schema") != "trsm_sidm_limits_v1":
        raise ValueError(f"unsupported limit table schema in {path}: {table.get('schema')!r}")
    if selected not in table["sets"]:
        raise ValueError(f"unknown limit set {selected!r}; available: {sorted(table['sets'])}")
    limits, seen = [], set()
    for entry in table["limits"]:
        if selected not in entry["sets"]:
            continue
        if entry["id"] in seen:
            raise ValueError(f"duplicate limit id {entry['id']!r}")
        seen.add(entry["id"])
        measure = entry["measure"]
        common = dict(
            id=entry["id"], label=entry["label"], limit_cm2_g=float(entry["limit_cm2_g"]),
            confidence_level=entry["confidence_level"], kind=measure["type"],
            reference=entry["reference"], sets=tuple(entry["sets"]),
            velocity_assumed=bool(entry.get("velocity_assumed", False)), notes=entry.get("notes", ""))
        if not common["limit_cm2_g"] > 0.0:
            raise ValueError(f"limit {entry['id']} must be positive")
        if measure["type"] == "maxwell":
            moment = int(measure.get("moment", 1) if moment_override is None else moment_override)
            limits.append(Limit(mean_relative_velocity_kms=float(measure["mean_relative_velocity_kms"]),
                                moment=moment, **common))
        elif measure["type"] == "fixed":
            lo, hi = (float(x) for x in measure["velocity_range_kms"])
            if not 0.0 < lo <= hi:
                raise ValueError(f"bad velocity range in {entry['id']}")
            if measure.get("combine", "min") not in COMBINERS:
                raise ValueError(f"unknown combine rule in {entry['id']}")
            limits.append(Limit(velocity_range_kms=(lo, hi), combine=measure.get("combine", "min"), **common))
        else:
            raise ValueError(f"unknown measure type in {entry['id']}: {measure['type']!r}")
    if not limits:
        raise ValueError(f"limit set {selected!r} is empty")
    return limits, table


def model_values(point, limits, **options):
    """sigma_V/m [cm^2/g] reduced as each limit prescribes; returns (dict, quadrature error)."""
    values, error = {}, 0.0
    maxwell = [limit for limit in limits if limit.kind == "maxwell"]
    if maxwell:
        dispersions = sorted({limit.dispersion_kms for limit in maxwell})
        powers = sorted({limit.moment for limit in maxwell})
        averages, error = sidm.maxwell_averages(point, dispersions, powers, return_error=True, **options)
        for limit in maxwell:
            values[limit.id] = float(sidm.per_mass_cm2_g(averages[(limit.dispersion_kms, limit.moment)], point.M3))
    fixed = [limit for limit in limits if limit.kind == "fixed"]
    if fixed:
        grid = np.unique(np.concatenate([limit.velocities_kms for limit in fixed]))
        _, sigma_v = sidm.cross_sections(point, grid / sidm.C_KMS, **options)
        per_mass = sidm.per_mass_cm2_g(sigma_v, point.M3)
        for limit in fixed:
            pick = np.searchsorted(grid, limit.velocities_kms)
            values[limit.id] = float(COMBINERS[limit.combine](per_mass[pick]))
    return values, error


def unitarity_values(point, limits):
    """Largest sigma_V/m [cm^2/g] allowed by unitarity for each limit, and a range check.

    Uses ``sigma_V <= UNITARITY_WEIGHT * 8 pi/k^2`` with ``k = M3 v/2``.  For a
    Maxwellian, ``<v^(n-2)>/<v^n> = 1/(2 (n+1) sigma_1D^2)``.  The second return
    value is True if the typical momentum of every limit is below the lightest
    mediator mass.
    """
    values, k_top = {}, 0.0
    for limit in limits:
        if limit.kind == "maxwell":
            disp = limit.dispersion_kms / sidm.C_KMS
            sigma = 32.0 * math.pi / point.M3**2 / (2.0 * (limit.moment + 1) * disp**2)
            v_typical = limit.mean_relative_velocity_kms / sidm.C_KMS
        else:
            grid = limit.velocities_kms / sidm.C_KMS
            sigma = float(COMBINERS[limit.combine](32.0 * math.pi / (point.M3 * grid) ** 2))
            v_typical = float(grid.max())
        values[limit.id] = UNITARITY_WEIGHT * float(sidm.per_mass_cm2_g(sigma, point.M3))
        k_top = max(k_top, 0.5 * point.M3 * v_typical)
    return values, k_top <= min(point.M1, point.M2)


def abundance_fraction(omega, reference=DEFAULT_OMEGA_REFERENCE, rescale=True):
    """xi = min(1, omega/reference); 1 without rescaling; nan if omega is unusable."""
    if not rescale:
        return 1.0
    if omega is None or not math.isfinite(omega) or omega < 0.0:
        return math.nan
    return min(1.0, omega / reference)


def bound_columns(limits):
    columns = ["sidm_bounds_status", "sidm_bounds_reason", "sidm_bounds_method", "sidm_bounds_set",
               "sidm_regime", "sidm_reliable", "sidm_abundance_fraction", "sidm_abundance_power",
               "sidm_bounds_applicable", "sidm_bounds_quadrature_relerr"]
    for limit in limits:
        columns += [f"sidm_{limit.id}_sigma_over_m", f"sidm_{limit.id}_ratio"]
    return columns + ["sidm_bounds_max_ratio", "sidm_bounds_strongest",
                      "sidm_bounds_unitarity_ratio", "sidm_excluded"]


def assess_point(point, omega, limits, *, set_name=DEFAULT_SET, rescale=True,
                 omega_reference=DEFAULT_OMEGA_REFERENCE, abundance_power=DEFAULT_ABUNDANCE_POWER,
                 min_fraction=DEFAULT_MIN_FRACTION, alpha_strong=sidm.DEFAULT_ALPHA_STRONG,
                 yukawa_threshold=sidm.DEFAULT_YUKAWA_THRESHOLD, **options):
    """Compare one point with ``limits``; returns the dict of output columns.

    ``sidm_<id>_sigma_over_m`` is the model value before rescaling and
    ``sidm_<id>_ratio`` is ``xi**p * value / limit``.
    """
    options = dict(options, yukawa_threshold=yukawa_threshold)
    regime = point.regime(yukawa_threshold, alpha_strong)
    reliable = regime != "strong"
    xi = abundance_fraction(omega, omega_reference, rescale)
    values, error = model_values(point, limits, **options)
    out = {"sidm_bounds_status": "ok", "sidm_bounds_reason": "", "sidm_bounds_method": BOUNDS_METHOD,
           "sidm_bounds_set": set_name, "sidm_regime": regime, "sidm_reliable": reliable,
           "sidm_abundance_fraction": xi, "sidm_abundance_power": abundance_power if rescale else 0.0,
           "sidm_bounds_applicable": math.nan, "sidm_bounds_quadrature_relerr": error}
    factor = xi**abundance_power if rescale else 1.0
    ratios = {}
    for limit in limits:
        ratios[limit.id] = factor * values[limit.id] / limit.limit_cm2_g
        out[f"sidm_{limit.id}_sigma_over_m"] = values[limit.id]
        out[f"sidm_{limit.id}_ratio"] = ratios[limit.id]
    finite = {key: value for key, value in ratios.items() if math.isfinite(value)}
    strongest = max(finite, key=finite.get) if finite else ""
    out["sidm_bounds_max_ratio"] = finite[strongest] if finite else math.nan
    out["sidm_bounds_strongest"] = strongest
    ceilings, short_range = unitarity_values(point, limits)
    ceiling = max(factor * ceilings[limit.id] / limit.limit_cm2_g for limit in limits)
    out["sidm_bounds_unitarity_ratio"] = ceiling

    if not math.isfinite(xi):
        out.update({"sidm_bounds_status": "unavailable", "sidm_excluded": math.nan,
                    "sidm_bounds_reason": "relic abundance missing; use --no-rescale to ignore it"})
        return out
    applicable = xi >= min_fraction
    out["sidm_bounds_applicable"] = applicable
    if not applicable:
        out.update({"sidm_excluded": False, "sidm_bounds_status": "not_applicable",
                    "sidm_bounds_reason": f"abundance fraction below {min_fraction:g}"})
    elif not reliable and short_range and math.isfinite(ceiling) and ceiling < 1.0:
        out.update({"sidm_excluded": False,
                    "sidm_bounds_reason": "strong coupling, but unitarity keeps it below every limit"})
    elif not reliable:
        out.update({"sidm_excluded": math.nan, "sidm_bounds_status": "unavailable",
                    "sidm_bounds_reason": f"strong coupling (alpha >= {alpha_strong:g}): cross section unreliable"})
    else:
        out["sidm_excluded"] = bool(finite) and finite[strongest] > 1.0
        if out["sidm_excluded"]:
            out["sidm_bounds_reason"] = f"exceeds {strongest}"
    return out


def failed_columns(error, limits, set_name):
    out = {name: math.nan for name in bound_columns(limits)}
    out.update({"sidm_bounds_status": "failed", "sidm_bounds_method": BOUNDS_METHOD,
                "sidm_bounds_set": set_name, "sidm_regime": "", "sidm_reliable": "",
                "sidm_bounds_strongest": "",
                "sidm_bounds_reason": str(error).replace("\t", " ").replace("\n", " ")})
    return out


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def _assess_options(args):
    return {"set_name": args.set, "rescale": args.rescale, "omega_reference": args.omega_reference,
            "abundance_power": args.abundance_power, "min_fraction": args.min_fraction,
            "alpha_strong": args.alpha_strong, "yukawa_threshold": args.yukawa_threshold,
            "mode": args.mode, "lmax": args.lmax, "n_angle": args.n_angle, "yukawa": args.yukawa}


def configuration(args, limits, table):
    return {
        "method": BOUNDS_METHOD, "cross_section_method": sidm.SIDM_METHOD,
        "limits_file": str(Path(args.limits) if args.limits else DEFAULT_LIMITS_FILE),
        "limit_set": args.set, "moment_override": args.moment,
        "limits": [{"id": l.id, "limit_cm2_g": l.limit_cm2_g, "confidence_level": l.confidence_level,
                    "evaluation": l.describe_velocity(), "reference": l.reference,
                    "velocity_assumed": l.velocity_assumed} for l in limits],
        "conventions": table.get("conventions", {}),
        "abundance": {
            "rescale": args.rescale, "omega_column": args.omega_column,
            "omega_reference": args.omega_reference, "power": args.abundance_power,
            "min_fraction": args.min_fraction,
            "rule": "ratio = min(1, omega/reference)**power * value / limit; no bound below min_fraction",
        },
        "cross_section": {"mode": args.mode, "lmax": args.lmax, "n_angle": args.n_angle,
                          "yukawa": args.yukawa, "yukawa_threshold": args.yukawa_threshold,
                          "alpha_strong": args.alpha_strong, "M1_GeV": args.m1},
        "verdict": "sidm_excluded is True/False, or nan when unavailable (never to be read as passing)",
    }


def run_list(args):
    limits, table = load_limits(args.limits, args.set, args.moment)
    print(f"limit set '{args.set}': {table['sets'][args.set]}\n")
    for limit in limits:
        flag = " [velocity assumed]" if limit.velocity_assumed else ""
        print(f"{limit.id}\n    {limit.label}\n    sigma/m < {limit.limit_cm2_g:g} cm^2/g ({limit.confidence_level}); "
              f"{limit.describe_velocity()}{flag}; {limit.reference}")
        if args.notes:
            print(f"    {limit.notes}")


def run_scan(args):
    source, destination = Path(args.input).resolve(), Path(args.output).resolve()
    if destination == source:
        raise SystemExit("--output must differ from the input scan")
    limits, table = load_limits(args.limits, args.set, args.moment)
    columns = bound_columns(limits)
    header, rows = sidm.read_scan(source)
    if args.rescale and args.omega_column not in header:
        raise SystemExit(f"column {args.omega_column!r} not found; give --omega-column or --no-rescale")
    counts = {"excluded": 0, "not_excluded": 0, "unavailable": 0, "not_applicable": 0}

    def evaluate(row):
        try:
            point = sidm.SIDMPoint.from_scan_row(row, M1=args.m1, v=args.vev)
            result = assess_point(point, sidm._row_float(row, args.omega_column), limits, **_assess_options(args))
        except Exception as error:  # record and continue
            counts["unavailable"] += 1
            return failed_columns(error, limits, args.set), False
        verdict = result["sidm_excluded"]
        key = "unavailable" if isinstance(verdict, float) else ("excluded" if verdict else "not_excluded")
        counts[key] += 1
        counts["not_applicable"] += result["sidm_bounds_status"] == "not_applicable"
        return result, True

    n_ok = sidm.write_scan(destination, header, rows, columns, evaluate, progress=args.progress,
                           label="sidm-bounds")
    meta = configuration(args, limits, table)
    meta.update({"input": str(source), "input_sha256": sidm.sha256(source), "output": str(destination),
                 "rows": len(rows), "rows_evaluated": n_ok, "counts": counts, "columns": columns})
    meta_path = destination.with_name(destination.name + ".sidm-bounds.json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {destination} and {meta_path.name}")
    print(f"{len(rows)} rows: {counts['excluded']} excluded, {counts['not_excluded']} not excluded "
          f"({counts['not_applicable']} of them below the abundance threshold), "
          f"{counts['unavailable']} unavailable")


def run_point(args):
    limits, _ = load_limits(args.limits, args.set, args.moment)
    point = sidm.point_from_args(args)
    result = assess_point(point, args.omega, limits, **_assess_options(args))
    if args.json:
        print(json.dumps(result, indent=2, default=str))
        return
    for key, value in result.items():
        print(f"{key:44s} {sidm.format_value(value)}")


def _add_common(parser):
    parser.add_argument("--limits", help=f"limit table (default {DEFAULT_LIMITS_FILE.name} in datafiles/)")
    parser.add_argument("--set", default=DEFAULT_SET, help="named set of limits in the table (see 'list')")
    parser.add_argument("--moment", type=int, default=None,
                        help="override n in <sigma_V v^n>/<v^n> for all Maxwellian limits "
                             "(5 gives the effective cross section of Yang & Yu, arXiv:2205.03392)")


def _add_assessment(parser):
    sidm.add_physics_arguments(parser)
    parser.add_argument("--no-rescale", dest="rescale", action="store_false",
                        help="treat H3 as all of the dark matter")
    parser.add_argument("--omega-reference", type=float, default=DEFAULT_OMEGA_REFERENCE)
    parser.add_argument("--abundance-power", type=float, default=DEFAULT_ABUNDANCE_POWER,
                        help="p in xi**p * sigma/m (default 2)")
    parser.add_argument("--min-fraction", type=float, default=DEFAULT_MIN_FRACTION,
                        help="no bound is applied below this abundance fraction (default 0.1)")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    listing = sub.add_parser("list", help="print the selected limits")
    _add_common(listing)
    listing.add_argument("--notes", action="store_true", help="also print the notes of each entry")
    scan = sub.add_parser("scan", help="add bound columns to a scan TSV (written to a new file)")
    scan.add_argument("input")
    scan.add_argument("--output", required=True)
    scan.add_argument("--omega-column", default=DEFAULT_OMEGA_COLUMN)
    scan.add_argument("--progress", type=int, default=0, help="report every N rows")
    _add_common(scan)
    _add_assessment(scan)
    point = sub.add_parser("point", help="assess one point")
    sidm.add_point_arguments(point)
    point.add_argument("--omega", type=float, default=None, help="relic density Omega h^2 of H3")
    point.add_argument("--json", action="store_true")
    _add_common(point)
    _add_assessment(point)
    args = parser.parse_args(argv)
    {"list": run_list, "scan": run_scan, "point": run_point}[args.command](args)


if __name__ == "__main__":
    main()
