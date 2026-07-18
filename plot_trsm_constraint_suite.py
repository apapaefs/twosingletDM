#!/usr/bin/env python3

"""Create a comprehensive constraint and diagnostic plot suite for TRSM scans.

The script intentionally lives beside, rather than inside, the simple
``plot_trsm_observables.py`` workflow.  It uses the stored scan flags as the
authoritative selections and provides fixed, reproducible styles for comparing
dark-matter and experimental constraints.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize, SymLogNorm, TwoSlopeNorm
from matplotlib.lines import Line2D


NOMINAL_MASS_GAP_GEV = 125.0
NOMINAL_M3_MAX_GEV = 1000.0
SM_LIKE_HIGGS_MASS_GEV = 125.09

BOOLEAN_COLUMNS = (
    "evo",
    "thc",
    "hb",
    "hs",
    "ewpo",
    "wmass",
    "dm",
    "dm_relic_excluded",
    "dm_direct_detection_excluded",
    "dm_indirect_available",
    "dm_indirect_detection_excluded",
)

OPTIONAL_BOOLEAN_COLUMNS = (
    "higgs_invisible_widths_included",
)

NUMERIC_COLUMNS = (
    "M2",
    "M3",
    "a12",
    "K133",
    "K233",
    "lSX",
    "dm_omega",
    "dm_relic_upper_limit",
    "dm_dir_det",
    "dm_dir_det_limit",
    "dm_indirect_ratio",
    "higgstools_hs_delta_chi2",
    "ewpt_ew_true_over_T",
)


class PlotUnavailable(RuntimeError):
    """Raised when a requested observable has no finite plottable values."""


@dataclass(frozen=True)
class ScanData:
    source: Path
    columns: tuple[str, ...]
    floats: dict[str, np.ndarray]
    bools: dict[str, np.ndarray]
    derived: dict[str, np.ndarray]

    def __len__(self) -> int:
        return len(self.floats["M2"])

    def f(self, name: str) -> np.ndarray:
        if name in self.derived:
            return self.derived[name]
        return self.floats[name]

    def b(self, name: str) -> np.ndarray:
        if name in self.derived:
            return self.derived[name]
        return self.bools[name]


@dataclass(frozen=True)
class CategoryStyle:
    display: str
    color: str
    marker: str
    size: float
    alpha: float
    zorder: float
    edgecolor: str = "none"
    linewidth: float = 0.0


@dataclass(frozen=True)
class PlotSpec:
    stem: str
    title: str
    kind: str
    scheme: str | None = None
    x: str | None = None
    y: str | None = None
    value: str | None = None
    norm_kind: str | None = None
    cmap: str | None = None
    xlabel: str | None = None
    ylabel: str | None = None
    colorbar_label: str | None = None


@dataclass(frozen=True)
class SummaryRow:
    metric: str
    count: int
    denominator: int
    note: str = ""

    @property
    def percent(self) -> float:
        if self.denominator == 0:
            return math.nan
        return 100.0 * self.count / self.denominator


FOURWAY_STYLES = OrderedDict(
    [
        (
            "neither",
            CategoryStyle("Neither", "#BDBDBD", "o", 8.0, 0.18, 1.0),
        ),
        (
            "DM only",
            CategoryStyle("DM only", "#0072B2", "^", 21.0, 0.72, 2.0),
        ),
        (
            "experimental only",
            CategoryStyle(
                "Experimental only", "#E69F00", "s", 21.0, 0.72, 3.0
            ),
        ),
        (
            "both",
            CategoryStyle(
                "Both", "#009E73", "*", 52.0, 0.95, 4.0, "#202020", 0.35
            ),
        ),
    ]
)

BINARY_COLORS = {
    "dm": "#0072B2",
    "experimental": "#E69F00",
    "full_viability": "#009E73",
    "hb": "#E69F00",
    "hs": "#E69F00",
    "wmass": "#E69F00",
    "relic_pass": "#0072B2",
    "direct_pass": "#0072B2",
}


PLOT_SPECS = (
    PlotSpec(
        "01_dm_experimental_fourway_m2_m3",
        "Dark-matter and experimental constraints",
        "categorical_mass",
        scheme="fourway",
    ),
    PlotSpec(
        "02_dm_status_m2_m3",
        "Aggregate dark-matter constraint",
        "categorical_mass",
        scheme="dm",
    ),
    PlotSpec(
        "03_experimental_status_m2_m3",
        "Combined experimental constraint",
        "categorical_mass",
        scheme="experimental",
    ),
    PlotSpec(
        "04_full_viability_status_m2_m3",
        "Full viability",
        "categorical_mass",
        scheme="full_viability",
    ),
    PlotSpec(
        "05_higgsbounds_status_m2_m3",
        "HiggsBounds constraint",
        "categorical_mass",
        scheme="hb",
    ),
    PlotSpec(
        "06_higgssignals_status_m2_m3",
        "HiggsSignals constraint",
        "categorical_mass",
        scheme="hs",
    ),
    PlotSpec(
        "07_wmass_status_m2_m3",
        r"$W$-mass constraint",
        "categorical_mass",
        scheme="wmass",
    ),
    PlotSpec(
        "08_relic_density_status_m2_m3",
        "Relic-density constraint",
        "categorical_mass",
        scheme="relic_pass",
    ),
    PlotSpec(
        "09_direct_detection_status_m2_m3",
        "Direct-detection constraint",
        "categorical_mass",
        scheme="direct_pass",
    ),
    PlotSpec(
        "10_indirect_detection_status_m2_m3",
        "Indirect-detection coverage and exclusion",
        "categorical_mass",
        scheme="indirect",
    ),
    PlotSpec(
        "11_dm_failure_modes_m2_m3",
        "Dark-matter failure modes",
        "categorical_mass",
        scheme="dm_failure",
    ),
    PlotSpec(
        "12_log10_relic_ratio_m2_m3",
        "Relic-density ratio",
        "continuous_mass",
        value="log10_relic_ratio",
        norm_kind="threshold0",
        cmap="RdBu_r",
        colorbar_label=r"$\log_{10}(\Omega/\Omega_{\max})$",
    ),
    PlotSpec(
        "13_log10_direct_detection_ratio_m2_m3",
        "Direct-detection ratio",
        "continuous_mass",
        value="log10_direct_ratio",
        norm_kind="threshold0",
        cmap="RdBu_r",
        colorbar_label=r"$\log_{10}(\sigma_{\rm SI}/\sigma_{\rm limit})$",
    ),
    PlotSpec(
        "14_indirect_ratio_m2_m3",
        "Indirect-detection ratio (available points)",
        "continuous_mass",
        value="indirect_ratio_available",
        norm_kind="log_to_one",
        cmap="viridis",
        colorbar_label=r"$\Phi/\Phi_{\rm limit}$",
    ),
    PlotSpec(
        "15_abs_a12_m2_m3",
        "Scalar mixing magnitude",
        "continuous_mass",
        value="abs_a12",
        norm_kind="positive",
        cmap="viridis",
        colorbar_label=r"$|a_{12}|$",
    ),
    PlotSpec(
        "16_k233_m2_m3",
        r"Portal trilinear $K_{233}$",
        "continuous_mass",
        value="K233",
        norm_kind="signed",
        cmap="coolwarm",
        colorbar_label=r"$K_{233}$ [GeV]",
    ),
    PlotSpec(
        "17_lsx_m2_m3",
        r"Portal coupling $\lambda_{SX}$",
        "continuous_mass",
        value="lSX",
        norm_kind="signed",
        cmap="coolwarm",
        colorbar_label=r"$\lambda_{SX}$",
    ),
    PlotSpec(
        "18_higgssignals_delta_chi2_m2_m3",
        r"HiggsSignals $\Delta\chi^2$",
        "continuous_mass",
        value="higgstools_hs_delta_chi2",
        norm_kind="threshold4",
        cmap="RdBu_r",
        colorbar_label=r"HiggsSignals $\Delta\chi^2$ (pass $<4$)",
    ),
    PlotSpec(
        "19_m2_vs_abs_a12",
        "Mixing-angle selection",
        "categorical_xy",
        scheme="fourway",
        x="M2",
        y="abs_a12",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$|a_{12}|$",
    ),
    PlotSpec(
        "20_k133_vs_k233",
        "Dark-sector trilinears",
        "categorical_xy",
        scheme="fourway",
        x="K133",
        y="K233",
        xlabel=r"$K_{133}$ [GeV]",
        ylabel=r"$K_{233}$ [GeV]",
    ),
    PlotSpec(
        "21_relic_vs_direct_ratio",
        "Relic-density versus direct-detection ratios",
        "ratio_plane",
    ),
    PlotSpec(
        "22_constraint_counts",
        "Constraint and coverage counts",
        "bars",
    ),
)

PLOT_BY_STEM = {spec.stem: spec for spec in PLOT_SPECS}

DASHBOARDS = OrderedDict(
    [
        (
            "dashboard_status_summary",
            (
                "01_dm_experimental_fourway_m2_m3",
                "04_full_viability_status_m2_m3",
                "03_experimental_status_m2_m3",
                "05_higgsbounds_status_m2_m3",
                "06_higgssignals_status_m2_m3",
                "07_wmass_status_m2_m3",
            ),
        ),
        (
            "dashboard_dm_summary",
            (
                "02_dm_status_m2_m3",
                "08_relic_density_status_m2_m3",
                "09_direct_detection_status_m2_m3",
                "10_indirect_detection_status_m2_m3",
            ),
        ),
        (
            "dashboard_diagnostic_summary",
            (
                "12_log10_relic_ratio_m2_m3",
                "13_log10_direct_detection_ratio_m2_m3",
                "15_abs_a12_m2_m3",
                "16_k233_m2_m3",
                "17_lsx_m2_m3",
                "18_higgssignals_delta_chi2_m2_m3",
            ),
        ),
    ]
)

DASHBOARD_TITLES = {
    "dashboard_status_summary": "TRSM constraint-status summary",
    "dashboard_dm_summary": "TRSM dark-matter constraint summary",
    "dashboard_diagnostic_summary": "TRSM diagnostic maps",
}


def strict_bool(value: str, column: str = "value", row_number: int | None = None) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    location = f" on row {row_number}" if row_number is not None else ""
    raise ValueError(
        f"Invalid boolean for {column}{location}: {value!r}; expected 'True' or 'False'."
    )


def strict_float(value: str, column: str, row_number: int) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid numeric value for {column} on row {row_number}: {value!r}."
        ) from exc


def safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    result = np.full(np.broadcast_shapes(numerator.shape, denominator.shape), np.nan)
    numerator, denominator = np.broadcast_arrays(numerator, denominator)
    valid = (
        np.isfinite(numerator)
        & np.isfinite(denominator)
        & (numerator >= 0.0)
        & (denominator > 0.0)
    )
    np.divide(numerator, denominator, out=result, where=valid)
    return result


def positive_log10(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.isfinite(values) & (values > 0.0)
    result[valid] = np.log10(values[valid])
    return result


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    if not arrays:
        return np.array([], dtype=bool)
    result = np.ones(np.asarray(arrays[0]).shape, dtype=bool)
    for array in arrays:
        result &= np.isfinite(np.asarray(array, dtype=float))
    return result


def four_way_categories(dm: np.ndarray, experimental: np.ndarray) -> np.ndarray:
    dm = np.asarray(dm, dtype=bool)
    experimental = np.asarray(experimental, dtype=bool)
    categories = np.full(dm.shape, "neither", dtype=object)
    categories[dm & ~experimental] = "DM only"
    categories[~dm & experimental] = "experimental only"
    categories[dm & experimental] = "both"
    return categories


def dm_failure_categories(
    dm: np.ndarray, relic_excluded: np.ndarray, direct_excluded: np.ndarray
) -> np.ndarray:
    dm = np.asarray(dm, dtype=bool)
    relic_excluded = np.asarray(relic_excluded, dtype=bool)
    direct_excluded = np.asarray(direct_excluded, dtype=bool)
    categories = np.full(dm.shape, "other DM failure", dtype=object)
    categories[dm] = "pass"
    failed = ~dm
    categories[failed & relic_excluded & ~direct_excluded] = "relic only"
    categories[failed & ~relic_excluded & direct_excluded] = "direct only"
    categories[failed & relic_excluded & direct_excluded] = "relic + direct"
    return categories


def indirect_categories(available: np.ndarray, excluded: np.ndarray) -> np.ndarray:
    available = np.asarray(available, dtype=bool)
    excluded = np.asarray(excluded, dtype=bool)
    categories = np.full(available.shape, "unavailable", dtype=object)
    categories[available & ~excluded] = "allowed"
    categories[available & excluded] = "excluded"
    return categories


def load_scan(path: Path | str) -> ScanData:
    path = Path(path)
    with path.open(encoding="ascii", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Input file is empty: {path}") from exc

        if len(header) != len(set(header)):
            duplicates = sorted({name for name in header if header.count(name) > 1})
            raise ValueError(f"Duplicate input columns: {', '.join(duplicates)}")

        required = set(BOOLEAN_COLUMNS) | set(NUMERIC_COLUMNS)
        missing = sorted(required - set(header))
        if missing:
            raise ValueError(f"Missing required input columns: {', '.join(missing)}")

        index = {column: position for position, column in enumerate(header)}
        bool_buffers = {
            column: []
            for column in BOOLEAN_COLUMNS + OPTIONAL_BOOLEAN_COLUMNS
            if column in index
        }
        float_buffers = {column: [] for column in NUMERIC_COLUMNS}

        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"Row {row_number} has {len(row)} fields; expected {len(header)}."
                )
            for column in bool_buffers:
                bool_buffers[column].append(
                    strict_bool(row[index[column]], column, row_number)
                )
            for column in NUMERIC_COLUMNS:
                float_buffers[column].append(
                    strict_float(row[index[column]], column, row_number)
                )

    bools = {
        column: np.asarray(values, dtype=bool)
        for column, values in bool_buffers.items()
    }
    floats = {
        column: np.asarray(values, dtype=float)
        for column, values in float_buffers.items()
    }
    if len(floats["M2"]) == 0:
        raise ValueError(f"Input file has a header but no data rows: {path}")
    for column in OPTIONAL_BOOLEAN_COLUMNS:
        if column not in bools:
            # Legacy scans predate explicit invisible-width provenance.  Such
            # rows are intentionally treated as unmodelled, not as passing.
            bools[column] = np.zeros(len(floats["M2"]), dtype=bool)
    if not np.all(finite_mask(floats["M2"], floats["M3"])):
        raise ValueError("M2 and M3 must be finite for every scan row.")

    theory = bools["evo"] & bools["thc"]
    experimental = bools["hb"] & bools["hs"] & bools["ewpo"] & bools["wmass"]
    full_viability = theory & experimental & bools["dm"]
    relic_pass = ~bools["dm_relic_excluded"]
    direct_pass = ~bools["dm_direct_detection_excluded"]
    relic_ratio = safe_ratio(floats["dm_omega"], floats["dm_relic_upper_limit"])
    direct_ratio = safe_ratio(floats["dm_dir_det"], floats["dm_dir_det_limit"])
    indirect_ratio = np.where(
        bools["dm_indirect_available"]
        & np.isfinite(floats["dm_indirect_ratio"])
        & (floats["dm_indirect_ratio"] > 0.0),
        floats["dm_indirect_ratio"],
        np.nan,
    )

    derived = {
        "theory": theory,
        "experimental": experimental,
        "full_viability": full_viability,
        "relic_pass": relic_pass,
        "direct_pass": direct_pass,
        "fourway": four_way_categories(bools["dm"], experimental),
        "dm_failure": dm_failure_categories(
            bools["dm"],
            bools["dm_relic_excluded"],
            bools["dm_direct_detection_excluded"],
        ),
        "indirect": indirect_categories(
            bools["dm_indirect_available"],
            bools["dm_indirect_detection_excluded"],
        ),
        "relic_ratio": relic_ratio,
        "direct_ratio": direct_ratio,
        "indirect_ratio_available": indirect_ratio,
        "log10_relic_ratio": positive_log10(relic_ratio),
        "log10_direct_ratio": positive_log10(direct_ratio),
        "log10_indirect_ratio": positive_log10(indirect_ratio),
        "abs_a12": np.abs(floats["a12"]),
    }

    return ScanData(
        source=path,
        columns=tuple(header),
        floats=floats,
        bools=bools,
        derived=derived,
    )


def robust_threshold_norm(
    values: np.ndarray, center: float = 0.0, lower: float = 1.0, upper: float = 99.0
) -> TwoSlopeNorm:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise PlotUnavailable("no finite values for threshold-centered normalization")
    vmin, vmax = np.percentile(finite, [lower, upper])
    scale = max(abs(vmin - center), abs(vmax - center), 1.0e-9)
    if not vmin < center:
        vmin = center - scale
    if not vmax > center:
        vmax = center + scale
    if math.isclose(vmin, center):
        vmin = center - scale
    if math.isclose(vmax, center):
        vmax = center + scale
    return TwoSlopeNorm(vmin=float(vmin), vcenter=float(center), vmax=float(vmax))


def robust_linear_norm(values: np.ndarray, include_zero: bool = False) -> Normalize:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise PlotUnavailable("no finite values for linear normalization")
    vmin, vmax = np.percentile(finite, [1.0, 99.0])
    if include_zero:
        vmin = min(0.0, float(vmin))
    if math.isclose(float(vmin), float(vmax)):
        span = max(abs(float(vmin)), 1.0) * 0.05
        vmin -= span
        vmax += span
    return Normalize(vmin=float(vmin), vmax=float(vmax), clip=True)


def robust_symlog_parameters(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise PlotUnavailable("no finite values for symmetric-log normalization")
    absolute = np.abs(finite)
    vmax = float(np.percentile(absolute, 99.0))
    if vmax <= 0.0:
        vmax = 1.0
    nonzero = absolute[absolute > 0.0]
    if nonzero.size == 0:
        linthresh = vmax * 0.01
    else:
        linthresh = float(np.percentile(nonzero, 5.0))
        linthresh = min(max(linthresh, vmax * 1.0e-6), vmax * 0.25)
    return linthresh, vmax


def robust_symlog_norm(values: np.ndarray) -> SymLogNorm:
    linthresh, vmax = robust_symlog_parameters(values)
    return SymLogNorm(
        linthresh=linthresh,
        vmin=-vmax,
        vmax=vmax,
        base=10,
        clip=True,
    )


def sparse_symlog_ticks(norm: SymLogNorm, max_per_side: int = 3) -> list[float]:
    maximum = max(abs(float(norm.vmin)), abs(float(norm.vmax)))
    if maximum <= 0.0:
        return [0.0]
    first_exponent = math.ceil(math.log10(float(norm.linthresh))) + 1
    last_exponent = math.floor(math.log10(maximum))
    exponents = list(range(first_exponent, last_exponent + 1))
    if not exponents:
        positive = [maximum]
    elif len(exponents) <= max_per_side:
        positive = [10.0**exponent for exponent in exponents]
    else:
        indices = np.linspace(0, len(exponents) - 1, max_per_side)
        selected = sorted({exponents[int(round(index))] for index in indices})
        positive = [10.0**exponent for exponent in selected]
    return [-value for value in reversed(positive)] + [0.0] + positive


def robust_log_norm_to_one(values: np.ndarray) -> LogNorm:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite) & (finite > 0.0)]
    if finite.size == 0:
        raise PlotUnavailable("no finite positive values for logarithmic normalization")
    vmin = float(np.percentile(finite, 1.0))
    vmax = max(1.0, float(np.percentile(finite, 99.0)))
    if not vmin < vmax:
        vmin = max(float(np.min(finite)), vmax * 0.1)
    if not vmin < vmax:
        vmin = vmax * 0.1
    return LogNorm(vmin=vmin, vmax=vmax, clip=True)


def norm_for(spec: PlotSpec, values: np.ndarray):
    if spec.norm_kind == "threshold0":
        return robust_threshold_norm(values, center=0.0)
    if spec.norm_kind == "threshold4":
        return robust_threshold_norm(values, center=4.0)
    if spec.norm_kind == "positive":
        return robust_linear_norm(values, include_zero=True)
    if spec.norm_kind == "signed":
        return robust_symlog_norm(values)
    if spec.norm_kind == "log_to_one":
        return robust_log_norm_to_one(values)
    return robust_linear_norm(values)


def category_styles(data: ScanData, scheme: str):
    if scheme == "fourway":
        return data.derived["fourway"], FOURWAY_STYLES
    if scheme == "indirect":
        styles = OrderedDict(
            [
                (
                    "unavailable",
                    CategoryStyle("Unavailable", "#BDBDBD", "o", 8.0, 0.18, 1.0),
                ),
                (
                    "allowed",
                    CategoryStyle("Available / allowed", "#0072B2", "^", 22.0, 0.75, 2.0),
                ),
                (
                    "excluded",
                    CategoryStyle("Available / excluded", "#D55E00", "X", 38.0, 0.9, 3.0),
                ),
            ]
        )
        return data.derived["indirect"], styles
    if scheme == "dm_failure":
        styles = OrderedDict(
            [
                (
                    "relic + direct",
                    CategoryStyle("Relic + direct fail", "#6F6F6F", "X", 17.0, 0.3, 1.0),
                ),
                (
                    "relic only",
                    CategoryStyle("Relic-only fail", "#CC79A7", "^", 22.0, 0.72, 2.0),
                ),
                (
                    "direct only",
                    CategoryStyle("Direct-only fail", "#D55E00", "s", 22.0, 0.72, 3.0),
                ),
                (
                    "pass",
                    CategoryStyle("DM pass", "#0072B2", "*", 48.0, 0.92, 4.0, "#202020", 0.3),
                ),
            ]
        )
        if np.any(data.derived["dm_failure"] == "other DM failure"):
            styles["other DM failure"] = CategoryStyle(
                "Other DM failure", "#F0E442", "D", 28.0, 0.9, 3.5, "#202020", 0.3
            )
        return data.derived["dm_failure"], styles

    mask = data.b(scheme)
    color = BINARY_COLORS[scheme]
    categories = np.where(mask, "pass", "fail")
    styles = OrderedDict(
        [
            (
                "fail",
                CategoryStyle("Fail", "#BDBDBD", "o", 8.0, 0.18, 1.0),
            ),
            (
                "pass",
                CategoryStyle("Pass", color, "D", 24.0, 0.78, 2.0, "#202020", 0.25),
            ),
        ]
    )
    return categories, styles


def mass_limits(data: ScanData) -> tuple[tuple[float, float], tuple[float, float]]:
    m2 = data.f("M2")
    m3 = data.f("M3")
    xspan = max(float(np.max(m2) - np.min(m2)), 1.0)
    yspan = max(float(np.max(m3) - np.min(m3)), 1.0)
    return (
        (float(np.min(m2) - 0.035 * xspan), float(np.max(m2) + 0.035 * xspan)),
        (float(np.min(m3) - 0.035 * yspan), float(np.max(m3) + 0.035 * yspan)),
    )


def draw_mass_guides(ax, data: ScanData, annotate: bool = False) -> None:
    (xmin, xmax), (ymin, ymax) = mass_limits(data)
    guide_min = max(xmin, ymin - NOMINAL_MASS_GAP_GEV)
    guide_max = min(xmax, ymax - NOMINAL_MASS_GAP_GEV)
    if guide_min < guide_max:
        guide_x = np.linspace(guide_min, guide_max, 200)
        ax.plot(
            guide_x,
            guide_x + NOMINAL_MASS_GAP_GEV,
            color="#666666",
            linestyle="--",
            linewidth=0.8,
            alpha=0.65,
            zorder=0.1,
        )
    ax.axhline(
        NOMINAL_M3_MAX_GEV,
        color="#666666",
        linestyle=":",
        linewidth=0.9,
        alpha=0.7,
        zorder=0.1,
    )
    if annotate:
        tail = int(np.count_nonzero(data.f("M3") > NOMINAL_M3_MAX_GEV))
        ax.text(
            0.985,
            0.018,
            (
                r"References: legacy nominal $M_3=M_2+125$ GeV (dashed), "
                r"nominal $M_3=1000$ GeV maximum (dotted)"
                f"\n{tail:,} points above the nominal $M_3$ maximum"
            ),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7.2,
            color="#444444",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 2.0},
            zorder=10,
        )


def style_mass_axis(ax, data: ScanData, annotate_guides: bool = False) -> None:
    draw_mass_guides(ax, data, annotate=annotate_guides)
    (xmin, xmax), (ymin, ymax) = mass_limits(data)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"$M_2$ [GeV]")
    ax.set_ylabel(r"$M_3$ [GeV]")
    ax.grid(True, alpha=0.18, linewidth=0.6)


def category_legend_handles(
    categories: np.ndarray,
    styles: OrderedDict[str, CategoryStyle],
    denominator: int,
    neutral_colors: bool = False,
) -> list[Line2D]:
    handles = []
    for key, style in styles.items():
        count = int(np.count_nonzero(categories == key))
        percent = 100.0 * count / denominator if denominator else 0.0
        facecolor = "#D9D9D9" if neutral_colors else style.color
        handles.append(
            Line2D(
                [0],
                [0],
                linestyle="None",
                marker=style.marker,
                markersize=max(4.0, math.sqrt(style.size)),
                markerfacecolor=facecolor,
                markeredgecolor=style.edgecolor if style.edgecolor != "none" else facecolor,
                markeredgewidth=max(style.linewidth, 0.35),
                color="none",
                label=f"{style.display}: {count:,} ({percent:.2f}%)",
            )
        )
    return handles


def render_categorical_mass(
    ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    categories, styles = category_styles(data, spec.scheme)
    m2 = data.f("M2")
    m3 = data.f("M3")
    valid = finite_mask(m2, m3)
    for key, style in styles.items():
        mask = valid & (categories == key)
        if not np.any(mask):
            continue
        ax.scatter(
            m2[mask],
            m3[mask],
            s=style.size * (0.82 if compact else 1.0),
            c=style.color,
            marker=style.marker,
            alpha=style.alpha,
            edgecolors=style.edgecolor,
            linewidths=style.linewidth,
            rasterized=True,
            zorder=style.zorder,
        )
    style_mass_axis(
        ax,
        data,
        annotate_guides=(spec.stem == "01_dm_experimental_fourway_m2_m3" and not compact),
    )
    if spec.scheme == "fourway":
        ax.set_title(
            spec.title
            + "\n"
            + r"Experimental = HB $\wedge$ HS $\wedge$ EWPO $\wedge$ $M_W$",
            fontsize=9.5 if compact else 12.0,
        )
    else:
        ax.set_title(spec.title, fontsize=9.5 if compact else 12.0)
    handles = category_legend_handles(categories[valid], styles, int(np.count_nonzero(valid)))
    ax.legend(
        handles=handles,
        loc="best",
        frameon=True,
        framealpha=0.82,
        edgecolor="none",
        fontsize=6.7 if compact else 8.0,
        handletextpad=0.5,
        borderpad=0.45,
    )


def render_continuous_mass(
    fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    values = data.f(spec.value)
    m2 = data.f("M2")
    m3 = data.f("M3")
    valid = finite_mask(m2, m3, values)
    if not np.any(valid):
        raise PlotUnavailable(f"{spec.value} has no finite values")

    norm = norm_for(spec, values[valid])
    cmap = plt.get_cmap(spec.cmap)
    categories = data.derived["fourway"]
    for key, style in FOURWAY_STYLES.items():
        mask = valid & (categories == key)
        if not np.any(mask):
            continue
        group_values = values[mask]
        if spec.norm_kind in {"threshold0", "threshold4"}:
            center = 0.0 if spec.norm_kind == "threshold0" else 4.0
            order = np.argsort(np.abs(group_values - center), kind="stable")
        elif spec.norm_kind == "signed":
            order = np.argsort(np.abs(group_values), kind="stable")
        else:
            order = np.argsort(group_values, kind="stable")
        indices = np.flatnonzero(mask)[order]
        ax.scatter(
            m2[indices],
            m3[indices],
            c=values[indices],
            norm=norm,
            cmap=cmap,
            s=max(7.0, style.size * (0.62 if compact else 0.72)),
            marker=style.marker,
            alpha=0.58 if key == "neither" else min(style.alpha + 0.05, 0.95),
            edgecolors="#202020" if key == "both" else "none",
            linewidths=0.25 if key == "both" else 0.0,
            rasterized=True,
            zorder=style.zorder,
        )

    style_mass_axis(ax, data)
    ax.set_title(spec.title, fontsize=9.5 if compact else 12.0)
    ax.text(
        0.985,
        0.018,
        f"finite N = {int(np.count_nonzero(valid)):,}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.7 if compact else 7.6,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
        zorder=10,
    )
    scalar_mappable = ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(
        scalar_mappable,
        ax=ax,
        extend="both",
        fraction=0.048 if compact else 0.046,
        pad=0.025,
    )
    colorbar.set_label(spec.colorbar_label, fontsize=8.0 if compact else 9.5)
    colorbar.ax.tick_params(labelsize=7.0 if compact else 8.0)
    if spec.norm_kind == "signed":
        colorbar.set_ticks(sparse_symlog_ticks(norm))
    if spec.norm_kind == "threshold4":
        ticks = [tick for tick in colorbar.get_ticks() if norm.vmin <= tick <= norm.vmax]
        colorbar.set_ticks(sorted(set(ticks + [4.0])))
        colorbar.ax.axhline(4.0, color="#222222", linewidth=0.8, alpha=0.8)
    handles = category_legend_handles(
        categories[valid], FOURWAY_STYLES, int(np.count_nonzero(valid)), neutral_colors=True
    )
    ax.legend(
        handles=handles,
        loc="best",
        frameon=True,
        framealpha=0.78,
        edgecolor="none",
        fontsize=5.9 if compact else 7.0,
        handletextpad=0.4,
        borderpad=0.4,
    )


def render_categorical_xy(
    ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    x = data.f(spec.x)
    y = data.f(spec.y)
    categories, styles = category_styles(data, spec.scheme)
    valid = finite_mask(x, y)
    if not np.any(valid):
        raise PlotUnavailable(f"{spec.x} and {spec.y} have no finite pairs")
    for key, style in styles.items():
        mask = valid & (categories == key)
        if not np.any(mask):
            continue
        ax.scatter(
            x[mask],
            y[mask],
            s=style.size,
            c=style.color,
            marker=style.marker,
            alpha=style.alpha,
            edgecolors=style.edgecolor,
            linewidths=style.linewidth,
            rasterized=True,
            zorder=style.zorder,
        )
    if spec.stem == "20_k133_vs_k233":
        x_linthresh, _ = robust_symlog_parameters(x[valid])
        y_linthresh, _ = robust_symlog_parameters(y[valid])
        ax.set_xscale("symlog", linthresh=x_linthresh)
        ax.set_yscale("symlog", linthresh=y_linthresh)
        x_max = float(np.max(np.abs(x[valid])))
        y_max = float(np.max(np.abs(y[valid])))
        ax.set_xlim(-1.08 * x_max, 1.08 * x_max)
        ax.set_ylim(-1.08 * y_max, 1.08 * y_max)
        ax.axhline(0.0, color="#666666", linewidth=0.7, alpha=0.55, zorder=0.1)
        ax.axvline(0.0, color="#666666", linewidth=0.7, alpha=0.55, zorder=0.1)
    ax.set_xlabel(spec.xlabel)
    ax.set_ylabel(spec.ylabel)
    ax.set_title(spec.title, fontsize=9.5 if compact else 12.0)
    ax.grid(True, alpha=0.18, linewidth=0.6)
    handles = category_legend_handles(categories[valid], styles, int(np.count_nonzero(valid)))
    ax.legend(
        handles=handles,
        loc="best",
        frameon=True,
        framealpha=0.82,
        edgecolor="none",
        fontsize=6.8 if compact else 8.0,
    )


def render_ratio_plane(ax, data: ScanData, spec: PlotSpec, compact: bool = False) -> None:
    relic = data.f("relic_ratio")
    direct = data.f("direct_ratio")
    valid = finite_mask(relic, direct) & (relic > 0.0) & (direct > 0.0)
    if not np.any(valid):
        raise PlotUnavailable("relic/direct ratios have no finite positive pairs")

    ax.scatter(
        relic[valid],
        direct[valid],
        s=8.0,
        c="#AFAFAF",
        marker="o",
        alpha=0.24,
        edgecolors="none",
        rasterized=True,
        zorder=1.0,
        label=f"All scan points: {int(np.count_nonzero(valid)):,} (100.00%)",
    )
    experimental = valid & data.b("experimental")
    ax.scatter(
        relic[experimental],
        direct[experimental],
        s=22.0,
        facecolors="none",
        edgecolors="#E69F00",
        linewidths=0.55,
        marker="o",
        alpha=0.76,
        rasterized=True,
        zorder=2.0,
        label=(
            f"Experimental pass: {int(np.count_nonzero(experimental)):,} "
            f"({100.0 * np.count_nonzero(experimental) / np.count_nonzero(valid):.2f}%)"
        ),
    )
    full = valid & data.b("full_viability")
    ax.scatter(
        relic[full],
        direct[full],
        s=52.0,
        c="#009E73",
        edgecolors="#202020",
        linewidths=0.3,
        marker="*",
        alpha=0.95,
        rasterized=True,
        zorder=3.0,
        label=(
            f"Full viable: {int(np.count_nonzero(full)):,} "
            f"({100.0 * np.count_nonzero(full) / np.count_nonzero(valid):.2f}%)"
        ),
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axvline(1.0, color="#222222", linestyle="--", linewidth=1.0, zorder=0.5)
    ax.axhline(1.0, color="#222222", linestyle="--", linewidth=1.0, zorder=0.5)
    ax.set_xlabel(r"$\Omega/\Omega_{\max}$")
    ax.set_ylabel(r"$\sigma_{\rm SI}/\sigma_{\rm limit}$")
    ax.set_title(spec.title, fontsize=9.5 if compact else 12.0)
    ax.grid(True, which="both", alpha=0.16, linewidth=0.55)
    ax.legend(frameon=True, framealpha=0.82, edgecolor="none", fontsize=7.5)

    failure = data.derived["dm_failure"]
    annotations = (
        ("pass", 0.02, 0.02, "Relic + direct pass"),
        ("relic only", 0.60, 0.02, "Relic-only fail"),
        ("direct only", 0.02, 0.91, "Direct-only fail"),
        ("relic + direct", 0.60, 0.91, "Both fail"),
    )
    for key, xpos, ypos, label in annotations:
        count = int(np.count_nonzero(failure == key))
        ax.text(
            xpos,
            ypos,
            f"{label}: {count:,}",
            transform=ax.transAxes,
            fontsize=7.2 if compact else 8.2,
            va="bottom" if ypos < 0.5 else "top",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.76, "pad": 1.6},
            zorder=10,
        )


def constraint_bar_metrics(data: ScanData):
    return (
        ("HiggsBounds pass", int(np.count_nonzero(data.b("hb"))), "#E69F00", ""),
        ("HiggsSignals pass", int(np.count_nonzero(data.b("hs"))), "#E69F00", ""),
        ("EWPO pass", int(np.count_nonzero(data.b("ewpo"))), "#E69F00", ""),
        (r"$W$ mass pass", int(np.count_nonzero(data.b("wmass"))), "#E69F00", ""),
        ("Experimental pass", int(np.count_nonzero(data.b("experimental"))), "#E69F00", ""),
        ("Relic-density pass", int(np.count_nonzero(data.b("relic_pass"))), "#0072B2", ""),
        ("Direct-detection pass", int(np.count_nonzero(data.b("direct_pass"))), "#0072B2", ""),
        (
            "Indirect result available",
            int(np.count_nonzero(data.b("dm_indirect_available"))),
            "#56B4E9",
            "//",
        ),
        ("Aggregate DM pass", int(np.count_nonzero(data.b("dm"))), "#0072B2", ""),
        ("Full viable", int(np.count_nonzero(data.b("full_viability"))), "#009E73", ""),
    )


def render_bars(ax, data: ScanData, spec: PlotSpec, compact: bool = False) -> None:
    metrics = constraint_bar_metrics(data)
    labels = [metric[0] for metric in metrics]
    counts = [metric[1] for metric in metrics]
    colors = [metric[2] for metric in metrics]
    hatches = [metric[3] for metric in metrics]
    positions = np.arange(len(labels))
    bars = ax.barh(positions, counts, color=colors, alpha=0.82, edgecolor="#333333", linewidth=0.25)
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ax.set_yticks(positions, labels)
    ax.invert_yaxis()
    ax.set_xlim(0.0, len(data) * 1.13)
    ax.set_xlabel(f"Points (total N = {len(data):,})")
    ax.set_title(spec.title, fontsize=9.5 if compact else 12.0)
    ax.grid(True, axis="x", alpha=0.18, linewidth=0.6)
    for position, count in zip(positions, counts):
        ax.text(
            count + len(data) * 0.012,
            position,
            f"{count:,} ({100.0 * count / len(data):.2f}%)",
            va="center",
            fontsize=7.2 if compact else 8.3,
        )
    ax.text(
        0.99,
        0.01,
        "Hatched bar denotes availability, not exclusion passing.",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.2,
        color="#444444",
    )


def render_spec(fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False) -> None:
    if spec.kind == "categorical_mass":
        render_categorical_mass(ax, data, spec, compact=compact)
    elif spec.kind == "continuous_mass":
        render_continuous_mass(fig, ax, data, spec, compact=compact)
    elif spec.kind == "categorical_xy":
        render_categorical_xy(ax, data, spec, compact=compact)
    elif spec.kind == "ratio_plane":
        render_ratio_plane(ax, data, spec, compact=compact)
    elif spec.kind == "bars":
        render_bars(ax, data, spec, compact=compact)
    else:
        raise ValueError(f"Unknown plot kind: {spec.kind}")


def extensions_for_format(plot_format: str) -> tuple[str, ...]:
    if plot_format == "both":
        return ("png", "pdf")
    return (plot_format,)


def all_figure_stems() -> tuple[str, ...]:
    return tuple(spec.stem for spec in PLOT_SPECS) + tuple(DASHBOARDS)


def expected_figure_paths(output_dir: Path | str, plot_format: str) -> list[Path]:
    output_dir = Path(output_dir)
    return [
        output_dir / f"{stem}.{extension}"
        for stem in all_figure_stems()
        for extension in extensions_for_format(plot_format)
    ]


def save_figure(
    fig,
    output_dir: Path,
    stem: str,
    plot_format: str,
    dpi: int,
) -> list[Path]:
    paths = []
    for extension in extensions_for_format(plot_format):
        path = output_dir / f"{stem}.{extension}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append(path)
        print(f"Saved {path}")
    return paths


def render_standalone(
    data: ScanData,
    spec: PlotSpec,
    output_dir: Path,
    plot_format: str,
    dpi: int,
) -> list[Path]:
    if spec.kind == "bars":
        figsize = (9.2, 6.3)
    else:
        figsize = (8.2, 6.2)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    try:
        render_spec(fig, ax, data, spec, compact=False)
        return save_figure(fig, output_dir, spec.stem, plot_format, dpi)
    finally:
        plt.close(fig)


def render_dashboard(
    data: ScanData,
    stem: str,
    plot_stems: Sequence[str],
    output_dir: Path,
    plot_format: str,
    dpi: int,
) -> list[Path]:
    if len(plot_stems) == 4:
        rows, columns = 2, 2
        figsize = (14.0, 10.0)
    else:
        rows, columns = 2, 3
        figsize = (18.0, 10.5)
    fig, axes = plt.subplots(rows, columns, figsize=figsize, constrained_layout=True)
    axes = np.asarray(axes).ravel()
    try:
        for panel_index, (ax, plot_stem) in enumerate(zip(axes, plot_stems)):
            spec = PLOT_BY_STEM[plot_stem]
            try:
                render_spec(fig, ax, data, spec, compact=True)
            except PlotUnavailable as exc:
                ax.axis("off")
                ax.text(0.5, 0.5, f"Unavailable\n{exc}", ha="center", va="center")
            ax.text(
                0.012,
                0.985,
                f"({chr(ord('a') + panel_index)})",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.0,
                fontweight="bold",
                zorder=20,
            )
        fig.suptitle(DASHBOARD_TITLES[stem], fontsize=16.0)
        return save_figure(fig, output_dir, stem, plot_format, dpi)
    finally:
        plt.close(fig)


def build_summary(data: ScanData, skipped_figures: Iterable[tuple[str, str]] = ()) -> list[SummaryRow]:
    n = len(data)
    rows = [SummaryRow("total_rows", n, n, "Input scan rows")]

    for name, note in (
        ("evo", "Evolution check pass"),
        ("thc", "Theory-constraint pass"),
        ("hb", "HiggsBounds pass"),
        ("hs", "HiggsSignals pass"),
        ("ewpo", "Electroweak-precision pass"),
        ("wmass", "W-mass pass"),
        ("theory", "evo & thc"),
        ("experimental", "hb & hs & ewpo & wmass"),
        ("dm", "Stored aggregate DM pass"),
        ("full_viability", "theory & experimental & dm"),
        ("relic_pass", "Not relic-density excluded"),
        ("direct_pass", "Not direct-detection excluded"),
    ):
        rows.append(SummaryRow(name, int(np.count_nonzero(data.b(name))), n, note))

    for category in FOURWAY_STYLES:
        rows.append(
            SummaryRow(
                f"fourway_{category.replace(' ', '_')}",
                int(np.count_nonzero(data.derived["fourway"] == category)),
                n,
                "DM/experimental four-way category",
            )
        )

    dm_failure_names = ["pass", "relic only", "direct only", "relic + direct"]
    if np.any(data.derived["dm_failure"] == "other DM failure"):
        dm_failure_names.append("other DM failure")
    for category in dm_failure_names:
        rows.append(
            SummaryRow(
                f"dm_failure_{category.replace(' ', '_').replace('+', 'and')}",
                int(np.count_nonzero(data.derived["dm_failure"] == category)),
                n,
                "Stored DM result split by relic/direct flags",
            )
        )

    for category in ("unavailable", "allowed", "excluded"):
        rows.append(
            SummaryRow(
                f"indirect_{category}",
                int(np.count_nonzero(data.derived["indirect"] == category)),
                n,
                "Indirect-detection availability/status",
            )
        )

    ewpt_finite = int(np.count_nonzero(np.isfinite(data.f("ewpt_ew_true_over_T"))))
    rows.append(
        SummaryRow(
            "ewpt_finite",
            ewpt_finite,
            n,
            "EWPT plots omitted because no finite values are present" if ewpt_finite == 0 else "Finite EWPT strengths",
        )
    )
    rows.append(
        SummaryRow(
            "m3_above_nominal_max",
            int(np.count_nonzero(data.f("M3") > NOMINAL_M3_MAX_GEV)),
            n,
            f"M3 > {NOMINAL_M3_MAX_GEV:g} GeV",
        )
    )
    rows.append(
        SummaryRow(
            "m2_in_reversed_sampling_domain",
            int(
                np.count_nonzero(
                    data.f("M2") > NOMINAL_M3_MAX_GEV - NOMINAL_MASS_GAP_GEV
                )
            ),
            n,
            "M2 > M3_max - 125 GeV, where the legacy conditional sampler reverses its endpoints",
        )
    )
    invisible_width_open = (
        (2.0 * data.f("M3") < SM_LIKE_HIGGS_MASS_GEV)
        | (2.0 * data.f("M3") < data.f("M2"))
    )
    rows.append(
        SummaryRow(
            "higgs_invisible_decay_open",
            int(np.count_nonzero(invisible_width_open)),
            n,
            "h1 or h2 -> h3 h3 is kinematically open",
        )
    )
    invisible_width_unmodelled = (
        invisible_width_open & ~data.b("higgs_invisible_widths_included")
    )
    provenance_note = (
        "Open invisible decay with widths not included in stored hb/hs"
        if "higgs_invisible_widths_included" in data.columns
        else "Open invisible decay in legacy input without width-model provenance"
    )
    rows.append(
        SummaryRow(
            "higgs_invisible_decay_open_but_unmodelled",
            int(np.count_nonzero(invisible_width_unmodelled)),
            n,
            provenance_note,
        )
    )

    component_dm = (
        ~data.bools["dm_relic_excluded"]
        & ~data.bools["dm_direct_detection_excluded"]
        & ~data.bools["dm_indirect_detection_excluded"]
    )
    rows.append(
        SummaryRow(
            "dm_component_mismatch",
            int(np.count_nonzero(component_dm != data.bools["dm"])),
            n,
            "Stored dm differs from conjunction of stored component exclusions",
        )
    )

    for name in ("evo", "thc", "ewpo"):
        values = data.b(name)
        if np.all(values == values[0]):
            rows.append(
                SummaryRow(
                    f"omitted_plot_{name}",
                    int(np.count_nonzero(values)),
                    n,
                    f"Standalone {name} plot omitted: flag is uniform {bool(values[0])}",
                )
            )
    for stem, reason in skipped_figures:
        rows.append(SummaryRow(f"skipped_figure_{stem}", 0, n, reason))
    return rows


def write_summary(path: Path, rows: Sequence[SummaryRow]) -> None:
    with path.open("w", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "count", "denominator", "percent", "note"])
        for row in rows:
            percent = "nan" if not math.isfinite(row.percent) else f"{row.percent:.8g}"
            writer.writerow([row.metric, row.count, row.denominator, percent, row.note])


def write_plot_index(
    path: Path,
    data: ScanData,
    figure_paths: Sequence[Path],
    summary_rows: Sequence[SummaryRow],
    skipped_figures: Iterable[tuple[str, str]] = (),
) -> None:
    """Write a self-contained HTML inventory for the generated plot suite."""
    files_by_stem: dict[str, dict[str, str]] = {}
    for figure_path in figure_paths:
        files_by_stem.setdefault(figure_path.stem, {})[
            figure_path.suffix.lstrip(".").lower()
        ] = figure_path.name
    skipped = dict(skipped_figures)

    def escaped(value: object) -> str:
        return html.escape(str(value), quote=True)

    def plot_card(stem: str, title: str) -> str:
        files = files_by_stem.get(stem, {})
        reason = skipped.get(stem)
        classes = "plot-card unavailable" if reason else "plot-card"
        parts = [
            f'<article class="{classes}" id="{escaped(stem)}">',
            f"<h3>{escaped(title)}</h3>",
            f'<p class="stem"><code>{escaped(stem)}</code></p>',
        ]
        if reason:
            parts.append(
                f'<div class="placeholder"><strong>Unavailable</strong><br>{escaped(reason)}</div>'
            )
        elif "png" in files:
            png = escaped(files["png"])
            parts.append(
                f'<a class="preview" href="{png}" target="_blank">'
                f'<img loading="lazy" src="{png}" alt="{escaped(title)}"></a>'
            )
        elif files:
            parts.append(
                '<div class="placeholder">Preview unavailable for PDF-only output.</div>'
            )
        else:
            parts.append('<div class="placeholder">No file generated for this plot.</div>')

        links = []
        for extension in ("png", "pdf"):
            if extension in files:
                filename = escaped(files[extension])
                links.append(
                    f'<a href="{filename}" target="_blank">{extension.upper()}</a>'
                )
        if links:
            parts.append(f'<p class="file-links">{" ".join(links)}</p>')
        parts.append("</article>")
        return "\n".join(parts)

    dashboard_cards = "\n".join(
        plot_card(stem, DASHBOARD_TITLES[stem]) for stem in DASHBOARDS
    )
    standalone_cards = "\n".join(
        plot_card(spec.stem, spec.title) for spec in PLOT_SPECS
    )
    summary_table_rows = []
    for row in summary_rows:
        percent = "&mdash;" if not math.isfinite(row.percent) else f"{row.percent:.6g}%"
        summary_table_rows.append(
            "<tr>"
            f"<td><code>{escaped(row.metric)}</code></td>"
            f"<td>{row.count:,}</td>"
            f"<td>{row.denominator:,}</td>"
            f"<td>{percent}</td>"
            f"<td>{escaped(row.note)}</td>"
            "</tr>"
        )

    experimental = int(np.count_nonzero(data.b("experimental")))
    dm_pass = int(np.count_nonzero(data.b("dm")))
    full = int(np.count_nonzero(data.b("full_viability")))
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TRSM constraint plots - {escaped(data.source.name)}</title>
<style>
:root {{ color-scheme: light; --ink: #17202a; --muted: #5f6b76; --line: #d9e0e6; --panel: #f7f9fb; --accent: #0072b2; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; color: var(--ink); background: #fff; font: 15px/1.5 system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }}
main {{ width: min(1500px, calc(100% - 32px)); margin: 0 auto 64px; }}
header {{ padding: 34px 0 20px; border-bottom: 1px solid var(--line); }}
h1 {{ margin: 0 0 8px; font-size: clamp(1.65rem, 3vw, 2.45rem); }}
h2 {{ margin-top: 38px; }}
h3 {{ margin: 0; font-size: 1rem; }}
p {{ margin: .45rem 0; }}
.meta, .stem {{ color: var(--muted); }}
nav {{ display: flex; flex-wrap: wrap; gap: 8px; margin-top: 18px; }}
nav a, .file-links a {{ color: var(--accent); text-decoration: none; border: 1px solid #b9d9ea; border-radius: 999px; padding: 5px 10px; background: #f3faff; }}
nav a:hover, .file-links a:hover {{ background: #e3f4fc; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 360px), 1fr)); gap: 20px; }}
.plot-card {{ min-width: 0; padding: 15px; border: 1px solid var(--line); border-radius: 10px; background: var(--panel); }}
.plot-card img {{ display: block; width: 100%; height: auto; margin-top: 12px; border: 1px solid var(--line); background: #fff; }}
.placeholder {{ display: grid; place-items: center; min-height: 210px; margin-top: 12px; padding: 24px; color: var(--muted); text-align: center; border: 1px dashed #aeb8c1; background: #fff; }}
.unavailable {{ opacity: .78; }}
.file-links {{ display: flex; gap: 8px; margin-top: 13px; }}
table {{ width: 100%; border-collapse: collapse; font-variant-numeric: tabular-nums; }}
th, td {{ padding: 8px 10px; border-bottom: 1px solid var(--line); text-align: left; vertical-align: top; }}
th {{ position: sticky; top: 0; background: #eef3f7; }}
.table-wrap {{ overflow-x: auto; border: 1px solid var(--line); border-radius: 8px; }}
footer {{ margin-top: 36px; color: var(--muted); }}
</style>
</head>
<body>
<main>
<header>
<h1>TRSM constraint plot suite</h1>
<p class="meta">Input: <code>{escaped(data.source)}</code> &middot; {len(data):,} rows &middot; experimental {experimental:,} &middot; DM {dm_pass:,} &middot; full viability {full:,}</p>
<nav><a href="#dashboards">Dashboards</a><a href="#standalone">Individual plots</a><a href="#summary">Constraint summary</a><a href="constraint_summary.tsv">Download TSV</a></nav>
</header>
<section id="dashboards"><h2>Dashboards</h2><div class="grid">{dashboard_cards}</div></section>
<section id="standalone"><h2>Individual plots</h2><div class="grid">{standalone_cards}</div></section>
<section id="summary"><h2>Constraint summary</h2><p><a href="constraint_summary.tsv">Download constraint_summary.tsv</a></p>
<div class="table-wrap"><table><thead><tr><th>Metric</th><th>Count</th><th>Denominator</th><th>Percent</th><th>Note</th></tr></thead>
<tbody>{''.join(summary_table_rows)}</tbody></table></div></section>
<footer>Generated by <code>plot_trsm_constraint_suite.py</code>.</footer>
</main>
</body>
</html>
"""
    path.write_text(document, encoding="utf-8")


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.0,
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def default_output_dir(input_path: Path) -> Path:
    return Path(__file__).resolve().parent / "plots" / f"{input_path.stem}_constraints"


def parse_args(argv: Sequence[str] | None = None):
    parser = argparse.ArgumentParser(
        description="Create the comprehensive TRSM constraint and diagnostic plot suite."
    )
    parser.add_argument("input", type=Path, help="TRSM scan TSV file.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory. Defaults to plots/<input-stem>_constraints beside this script.",
    )
    parser.add_argument(
        "--format", choices=("png", "pdf", "both"), default="both"
    )
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args(argv)
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    if args.output_dir is None:
        args.output_dir = default_output_dir(args.input)
    return args


def run(argv: Sequence[str] | None = None) -> list[Path]:
    args = parse_args(argv)
    configure_style()
    data = load_scan(args.input)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loaded {len(data):,} rows from {data.source}")
    print(
        "Selections: "
        f"experimental={np.count_nonzero(data.b('experimental')):,} "
        f"dm={np.count_nonzero(data.b('dm')):,} "
        f"full={np.count_nonzero(data.b('full_viability')):,}"
    )

    paths: list[Path] = []
    skipped: list[tuple[str, str]] = []
    for spec in PLOT_SPECS:
        try:
            paths.extend(
                render_standalone(data, spec, output_dir, args.format, args.dpi)
            )
        except PlotUnavailable as exc:
            reason = str(exc)
            skipped.append((spec.stem, reason))
            print(f"Skipped {spec.stem}: {reason}")

    for dashboard_stem, plot_stems in DASHBOARDS.items():
        paths.extend(
            render_dashboard(
                data,
                dashboard_stem,
                plot_stems,
                output_dir,
                args.format,
                args.dpi,
            )
        )

    summary_rows = build_summary(data, skipped)
    summary_path = output_dir / "constraint_summary.tsv"
    write_summary(summary_path, summary_rows)
    print(f"Saved {summary_path}")
    index_path = output_dir / "index.html"
    write_plot_index(index_path, data, paths, summary_rows, skipped)
    print(f"Saved {index_path}")
    print(f"Created {len(paths)} figure files in {output_dir}")
    return paths


if __name__ == "__main__":
    run()
