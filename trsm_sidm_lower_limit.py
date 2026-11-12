"""Lower limit on the H3 self-interaction at galactic velocities.

Self-interacting dark matter produces the constant-density cores inferred in
dwarf and low-surface-brightness galaxies if ``sigma/m >~ 1 cm^2/g`` at the
velocities of those halos.  This module tests that requirement for vx=0 TRSM
points.  It is separate from ``trsm_sidm_bounds.py`` (upper limits) and uses
``trsm_sidm.py`` for the cross section.

This is a requirement, not an exclusion.  A point that fails it is not ruled
out; its dark matter simply behaves as collisionless in galaxies.  The
analyses behind the number fix ``sigma/m = 3 cm^2/g`` and report that
values above about 1 give the same results (Ren et al., arXiv:1808.05695;
Kamada et al., arXiv:1611.02716); the explicit ``>~ 1 cm^2/g at
v_rel = 30-200 km/s`` is from Tulin & Yu (arXiv:1705.02358, Table I).

Evaluation
----------
For each entry of ``datafiles/sidm_lower_limits_v1.json`` the model value
is ``<sigma_V v^n>/<v^n>`` (n = 1) over a Maxwellian of mean relative
velocity ``<v> = 4 sigma_v0/sqrt(pi)``, on a logarithmic grid of ``<v>``
across the velocity range of the entry.  The requirement has to hold in all
those halos, so the smallest value over the range is compared with the
minimum.  The largest value is written as well.

Abundance
---------
Cores form only if the self-interacting species is the dark matter.  With
``xi = min(1, Omega h^2/0.12)`` the model value is multiplied by ``xi**p``
(p = 2, as for the upper limits), and the requirement is not met when
``xi < min_fraction`` (default 0.9), whatever the cross section: the
collisionless remainder keeps its cusp.  ``--no-rescale`` treats H3 as all
of the dark matter.

Verdict (column ``sidm_lower_satisfied``)
-----------------------------------------
``True``   every selected entry has rescaled ratio >= 1, xi >= min_fraction
           and the cross section is reliable.
``False``  an entry is not reached, or xi < min_fraction.  A point in the
           ``strong`` regime is also ``False`` if unitarity keeps it below
           an entry (``sidm_lower_unitarity_ratio < 1``).
``nan``    unavailable: the calculation failed, the abundance is missing, or
           the point is in the ``strong`` regime and unitarity does not
           decide.

Usage::

    python trsm_sidm_lower_limit.py list --set all --notes
    python trsm_sidm_lower_limit.py scan scan.dat --output scan_sidm_lower.dat
    python trsm_sidm_lower_limit.py point --m2 39.9999992 --m3 20 --vs 500 --a12 0 \
        --lx 0.1 --lphix 0 --lsx 0.02 --omega 0.12
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import trsm_sidm as sidm
from trsm_sidm_bounds import (
    DEFAULT_ABUNDANCE_POWER, DEFAULT_OMEGA_COLUMN, DEFAULT_OMEGA_REFERENCE, UNITARITY_WEIGHT,
    abundance_fraction,
)


LOWER_METHOD = "trsm_sidm_lower_limit_v1"
DEFAULT_TABLE = Path(__file__).resolve().parent / "datafiles" / "sidm_lower_limits_v1.json"
DEFAULT_SET = "default"
DEFAULT_MIN_FRACTION = 0.9
GRID_POINTS = 7


@dataclass(frozen=True)
class Requirement:
    """A minimum sigma/m over a range of mean relative velocities."""

    id: str
    label: str
    minimum_cm2_g: float
    velocity_range_kms: tuple
    moment: int
    reference: str
    sets: tuple
    velocity_assumed: bool = False
    notes: str = ""

    @property
    def mean_velocities_kms(self):
        lo, hi = self.velocity_range_kms
        return np.array([lo]) if hi == lo else np.geomspace(lo, hi, GRID_POINTS)

    @property
    def dispersions_kms(self):
        return [sidm.mean_relative_velocity_to_dispersion(v) for v in self.mean_velocities_kms]

    def describe_velocity(self):
        lo, hi = self.velocity_range_kms
        span = f"{lo:g}" if hi == lo else f"{lo:g}-{hi:g}"
        return f"<v> = {span} km/s (n = {self.moment})"


def load_requirements(path=None, selected=DEFAULT_SET, moment_override=None):
    """Return (requirements, table) for the named set of the table."""
    path = Path(path) if path else DEFAULT_TABLE
    table = json.loads(path.read_text(encoding="utf-8"))
    if table.get("schema") != "trsm_sidm_lower_limits_v1":
        raise ValueError(f"unsupported table schema in {path}: {table.get('schema')!r}")
    if selected not in table["sets"]:
        raise ValueError(f"unknown set {selected!r}; available: {sorted(table['sets'])}")
    requirements, seen = [], set()
    for entry in table["requirements"]:
        if selected not in entry["sets"]:
            continue
        if entry["id"] in seen:
            raise ValueError(f"duplicate id {entry['id']!r}")
        seen.add(entry["id"])
        lo, hi = (float(x) for x in entry["mean_relative_velocity_range_kms"])
        if not 0.0 < lo <= hi:
            raise ValueError(f"bad velocity range in {entry['id']}")
        minimum = float(entry["minimum_cm2_g"])
        if not minimum > 0.0:
            raise ValueError(f"minimum of {entry['id']} must be positive")
        requirements.append(Requirement(
            id=entry["id"], label=entry["label"], minimum_cm2_g=minimum, velocity_range_kms=(lo, hi),
            moment=int(entry.get("moment", 1) if moment_override is None else moment_override),
            reference=entry["reference"], sets=tuple(entry["sets"]),
            velocity_assumed=bool(entry.get("velocity_assumed", False)), notes=entry.get("notes", "")))
    if not requirements:
        raise ValueError(f"set {selected!r} is empty")
    return requirements, table


def model_values(point, requirements, **options):
    """(min, max) of the averaged sigma_V/m [cm^2/g] over each velocity range, and the quadrature error."""
    dispersions = sorted({d for requirement in requirements for d in requirement.dispersions_kms})
    powers = sorted({requirement.moment for requirement in requirements})
    averages, error = sidm.maxwell_averages(point, dispersions, powers, return_error=True, **options)
    values = {}
    for requirement in requirements:
        series = sidm.per_mass_cm2_g(
            np.array([averages[(d, requirement.moment)] for d in requirement.dispersions_kms]), point.M3)
        values[requirement.id] = (float(series.min()), float(series.max()))
    return values, error


def unitarity_values(point, requirements):
    """Largest value of the range minimum allowed by unitarity, and a range check.

    ``sigma_V <= UNITARITY_WEIGHT * 8 pi/k^2`` (even partial waves up to l = 8)
    falls like 1/v^2, so the minimum over the range sits at its upper end.
    The second return value is True if the typical momentum there is below
    the lightest mediator mass.
    """
    values, k_top = {}, 0.0
    for requirement in requirements:
        top = requirement.velocity_range_kms[1]
        disp = sidm.mean_relative_velocity_to_dispersion(top) / sidm.C_KMS
        sigma = 32.0 * math.pi / point.M3**2 / (2.0 * (requirement.moment + 1) * disp**2)
        values[requirement.id] = UNITARITY_WEIGHT * float(sidm.per_mass_cm2_g(sigma, point.M3))
        k_top = max(k_top, 0.5 * point.M3 * top / sidm.C_KMS)
    return values, k_top <= min(point.M1, point.M2)


def lower_columns(requirements):
    columns = ["sidm_lower_status", "sidm_lower_reason", "sidm_lower_method", "sidm_lower_set",
               "sidm_regime", "sidm_reliable", "sidm_lower_abundance_fraction",
               "sidm_lower_abundance_power", "sidm_lower_quadrature_relerr"]
    for requirement in requirements:
        columns += [f"sidm_lower_{requirement.id}_min", f"sidm_lower_{requirement.id}_max",
                    f"sidm_lower_{requirement.id}_ratio"]
    return columns + ["sidm_lower_min_ratio", "sidm_lower_weakest", "sidm_lower_unitarity_ratio",
                      "sidm_lower_satisfied"]


def assess_point(point, omega, requirements, *, set_name=DEFAULT_SET, rescale=True,
                 omega_reference=DEFAULT_OMEGA_REFERENCE, abundance_power=DEFAULT_ABUNDANCE_POWER,
                 min_fraction=DEFAULT_MIN_FRACTION, alpha_strong=sidm.DEFAULT_ALPHA_STRONG,
                 yukawa_threshold=sidm.DEFAULT_YUKAWA_THRESHOLD, **options):
    """Test one point against ``requirements``; returns the dict of output columns.

    ``sidm_lower_<id>_min`` and ``_max`` are the model values before
    rescaling; ``sidm_lower_<id>_ratio`` is ``xi**p * min / minimum``.
    """
    options = dict(options, yukawa_threshold=yukawa_threshold)
    regime = point.regime(yukawa_threshold, alpha_strong)
    reliable = regime != "strong"
    xi = abundance_fraction(omega, omega_reference, rescale)
    values, error = model_values(point, requirements, **options)
    factor = xi**abundance_power if rescale else 1.0
    out = {"sidm_lower_status": "ok", "sidm_lower_reason": "", "sidm_lower_method": LOWER_METHOD,
           "sidm_lower_set": set_name, "sidm_regime": regime, "sidm_reliable": reliable,
           "sidm_lower_abundance_fraction": xi,
           "sidm_lower_abundance_power": abundance_power if rescale else 0.0,
           "sidm_lower_quadrature_relerr": error}
    ratios = {}
    for requirement in requirements:
        low, high = values[requirement.id]
        ratios[requirement.id] = factor * low / requirement.minimum_cm2_g
        out[f"sidm_lower_{requirement.id}_min"] = low
        out[f"sidm_lower_{requirement.id}_max"] = high
        out[f"sidm_lower_{requirement.id}_ratio"] = ratios[requirement.id]
    finite = {key: value for key, value in ratios.items() if math.isfinite(value)}
    weakest = min(finite, key=finite.get) if finite else ""
    out["sidm_lower_min_ratio"] = finite[weakest] if finite else math.nan
    out["sidm_lower_weakest"] = weakest
    ceilings, short_range = unitarity_values(point, requirements)
    ceiling = min(factor * ceilings[r.id] / r.minimum_cm2_g for r in requirements)
    out["sidm_lower_unitarity_ratio"] = ceiling

    if not math.isfinite(xi):
        out.update({"sidm_lower_status": "unavailable", "sidm_lower_satisfied": math.nan,
                    "sidm_lower_reason": "relic abundance missing; use --no-rescale to ignore it"})
    elif xi < min_fraction:
        out.update({"sidm_lower_satisfied": False, "sidm_lower_status": "subdominant",
                    "sidm_lower_reason": f"abundance fraction below {min_fraction:g}"})
    elif reliable:
        satisfied = len(finite) == len(requirements) and finite[weakest] >= 1.0
        out["sidm_lower_satisfied"] = satisfied
        if not satisfied:
            out["sidm_lower_reason"] = f"below {weakest}"
    elif short_range and math.isfinite(ceiling) and ceiling < 1.0:
        out.update({"sidm_lower_satisfied": False,
                    "sidm_lower_reason": "strong coupling, but out of reach by unitarity"})
    else:
        out.update({"sidm_lower_satisfied": math.nan, "sidm_lower_status": "unavailable",
                    "sidm_lower_reason": f"strong coupling (alpha >= {alpha_strong:g}): cross section unreliable"})
    return out


def failed_columns(error, requirements, set_name):
    out = {name: math.nan for name in lower_columns(requirements)}
    out.update({"sidm_lower_status": "failed", "sidm_lower_method": LOWER_METHOD,
                "sidm_lower_set": set_name, "sidm_regime": "", "sidm_reliable": "",
                "sidm_lower_weakest": "",
                "sidm_lower_reason": str(error).replace("\t", " ").replace("\n", " ")})
    return out


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def _assess_options(args):
    return {"set_name": args.set, "rescale": args.rescale, "omega_reference": args.omega_reference,
            "abundance_power": args.abundance_power, "min_fraction": args.min_fraction,
            "alpha_strong": args.alpha_strong, "yukawa_threshold": args.yukawa_threshold,
            "mode": args.mode, "lmax": args.lmax, "n_angle": args.n_angle, "yukawa": args.yukawa}


def configuration(args, requirements, table):
    return {
        "method": LOWER_METHOD, "cross_section_method": sidm.SIDM_METHOD,
        "table": str(Path(args.table) if args.table else DEFAULT_TABLE),
        "set": args.set, "moment_override": args.moment,
        "requirements": [{"id": r.id, "minimum_cm2_g": r.minimum_cm2_g, "evaluation": r.describe_velocity(),
                          "reference": r.reference, "velocity_assumed": r.velocity_assumed}
                         for r in requirements],
        "conventions": table.get("conventions", {}),
        "abundance": {
            "rescale": args.rescale, "omega_column": args.omega_column,
            "omega_reference": args.omega_reference, "power": args.abundance_power,
            "min_fraction": args.min_fraction,
            "rule": "ratio = min(1, omega/reference)**power * (smallest value in range) / minimum; "
                    "not satisfied below min_fraction",
        },
        "cross_section": {"mode": args.mode, "lmax": args.lmax, "n_angle": args.n_angle,
                          "yukawa": args.yukawa, "yukawa_threshold": args.yukawa_threshold,
                          "alpha_strong": args.alpha_strong, "M1_GeV": args.m1},
        "verdict": "sidm_lower_satisfied is True/False, or nan when unavailable; False is not an exclusion",
    }


def run_list(args):
    requirements, table = load_requirements(args.table, args.set, args.moment)
    print(f"set '{args.set}': {table['sets'][args.set]}\n")
    for requirement in requirements:
        flag = " [velocity assumed]" if requirement.velocity_assumed else ""
        print(f"{requirement.id}\n    {requirement.label}\n    sigma/m >= {requirement.minimum_cm2_g:g} cm^2/g; "
              f"{requirement.describe_velocity()}{flag}; {requirement.reference}")
        if args.notes:
            print(f"    {requirement.notes}")


def run_scan(args):
    source, destination = Path(args.input).resolve(), Path(args.output).resolve()
    if destination == source:
        raise SystemExit("--output must differ from the input scan")
    requirements, table = load_requirements(args.table, args.set, args.moment)
    columns = lower_columns(requirements)
    header, rows = sidm.read_scan(source)
    if args.rescale and args.omega_column not in header:
        raise SystemExit(f"column {args.omega_column!r} not found; give --omega-column or --no-rescale")
    counts = {"satisfied": 0, "not_satisfied": 0, "unavailable": 0, "subdominant": 0}

    def evaluate(row):
        try:
            point = sidm.SIDMPoint.from_scan_row(row, M1=args.m1, v=args.vev)
            result = assess_point(point, sidm._row_float(row, args.omega_column), requirements,
                                  **_assess_options(args))
        except Exception as error:  # record and continue
            counts["unavailable"] += 1
            return failed_columns(error, requirements, args.set), False
        verdict = result["sidm_lower_satisfied"]
        key = "unavailable" if isinstance(verdict, float) else ("satisfied" if verdict else "not_satisfied")
        counts[key] += 1
        counts["subdominant"] += result["sidm_lower_status"] == "subdominant"
        return result, True

    n_ok = sidm.write_scan(destination, header, rows, columns, evaluate, progress=args.progress,
                           label="sidm-lower")
    meta = configuration(args, requirements, table)
    meta.update({"input": str(source), "input_sha256": sidm.sha256(source), "output": str(destination),
                 "rows": len(rows), "rows_evaluated": n_ok, "counts": counts, "columns": columns})
    meta_path = destination.with_name(destination.name + ".sidm-lower.json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {destination} and {meta_path.name}")
    print(f"{len(rows)} rows: {counts['satisfied']} satisfy the lower limit, {counts['not_satisfied']} do not "
          f"({counts['subdominant']} of them below the abundance threshold), {counts['unavailable']} unavailable")


def run_point(args):
    requirements, _ = load_requirements(args.table, args.set, args.moment)
    point = sidm.point_from_args(args)
    result = assess_point(point, args.omega, requirements, **_assess_options(args))
    if args.json:
        print(json.dumps(result, indent=2, default=str))
        return
    for key, value in result.items():
        print(f"{key:44s} {sidm.format_value(value)}")


def _add_common(parser):
    parser.add_argument("--table", help=f"requirement table (default {DEFAULT_TABLE.name} in datafiles/)")
    parser.add_argument("--set", default=DEFAULT_SET, help="named set of requirements in the table (see 'list')")
    parser.add_argument("--moment", type=int, default=None, help="override n in <sigma_V v^n>/<v^n>")


def _add_assessment(parser):
    sidm.add_physics_arguments(parser)
    parser.add_argument("--no-rescale", dest="rescale", action="store_false",
                        help="treat H3 as all of the dark matter")
    parser.add_argument("--omega-reference", type=float, default=DEFAULT_OMEGA_REFERENCE)
    parser.add_argument("--abundance-power", type=float, default=DEFAULT_ABUNDANCE_POWER,
                        help="p in xi**p * sigma/m (default 2)")
    parser.add_argument("--min-fraction", type=float, default=DEFAULT_MIN_FRACTION,
                        help="abundance fraction below which the requirement cannot be met (default 0.9)")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    listing = sub.add_parser("list", help="print the selected requirements")
    _add_common(listing)
    listing.add_argument("--notes", action="store_true", help="also print the notes of each entry")
    scan = sub.add_parser("scan", help="add lower-limit columns to a scan TSV (written to a new file)")
    scan.add_argument("input")
    scan.add_argument("--output", required=True)
    scan.add_argument("--omega-column", default=DEFAULT_OMEGA_COLUMN)
    scan.add_argument("--progress", type=int, default=0, help="report every N rows")
    _add_common(scan)
    _add_assessment(scan)
    point = sub.add_parser("point", help="test one point")
    sidm.add_point_arguments(point)
    point.add_argument("--omega", type=float, default=None, help="relic density Omega h^2 of H3")
    point.add_argument("--json", action="store_true")
    _add_common(point)
    _add_assessment(point)
    args = parser.parse_args(argv)
    {"list": run_list, "scan": run_scan, "point": run_point}[args.command](args)


if __name__ == "__main__":
    main()
