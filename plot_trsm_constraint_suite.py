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
import json
import math
import shlex
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import (
    LinearSegmentedColormap,
    LogNorm,
    Normalize,
    SymLogNorm,
    TwoSlopeNorm,
)
from matplotlib.lines import Line2D


NOMINAL_MASS_GAP_GEV = 125.0
NOMINAL_M3_MAX_GEV = 1000.0
SM_LIKE_HIGGS_MASS_GEV = 125.09
BSMPT_STRONG_EWPT_THRESHOLD = 1.0
M2_EQ_2M3_GUIDE_COLOR = "#CC79A7"
M3_EQ_2M2_GUIDE_COLOR = "#0072B2"
RESONANCE_CMAP = LinearSegmentedColormap.from_list(
    "trsm_resonance",
    ("#0072B2", "#171717", "#D55E00"),
)
SCAN_METADATA_SCHEMA = "trsm_scan_metadata_v1"
HL_LHC_LUMINOSITY_FB = 3000.0
YR4_SIGNAL_TABLE_NAME = "lhchxswg_yr4_bsm_13p6tev_ggf_vbf.tsv"
YR4_SIGNAL_SOURCE_URL = (
    "https://gitlab.cern.ch/LHCHIGGSXS/LHCHXSWG1/crosssections"
)
YR4_SIGNAL_CITATION_URL = "https://arxiv.org/abs/1610.07922"
YR4_SIGNAL_REPOSITORY_COMMIT = "aad67de39778537fa36fa0692abf66fd43f660a4"
YR4_SIGNAL_COLUMNS = (
    "mass_gev",
    "ggf_pb",
    "ggf_scale_up_pct",
    "ggf_scale_down_pct",
    "ggf_pdfalphas_pct",
    "vbf_pb",
    "vbf_scale_up_pct",
    "vbf_scale_down_pct",
    "vbf_pdfalphas_pct",
)
SIGNAL_REQUIRED_COLUMNS = ("k2", "w2", "h2_h3h3_br")

BOOLEAN_COLUMNS = (
    "evo",
    "thc",
    "hb",
    "hs",
    "ewpo",
    "wmass",
    "dm",
)

NULLABLE_BOOLEAN_COLUMNS = (
    "dm_relic_excluded",
    "dm_direct_detection_excluded",
    "dm_indirect_available",
    "dm_indirect_detection_excluded",
)

OPTIONAL_BOOLEAN_COLUMNS = (
    "higgs_invisible_widths_included",
)

OPTIONAL_NULLABLE_BOOLEAN_COLUMNS = (
    "ewpt_has_x_broken",
)

OPTIONAL_TEXT_COLUMNS = (
    "ewpt_status",
    "ewpt_error",
    "ewpt_global_phase_path",
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

OPTIONAL_NUMERIC_COLUMNS = (
    "w2",
    "ewpt_ew_step_index",
    "vs",
    "vx",
    "a13",
    "a23",
    "lX",
    "lPhiX",
    "k1",
    "k2",
    "k3",
    "h1_h2h2_width",
    "h1_h2h2_br",
    "h2_h1h1_width",
    "h2_h1h1_br",
    "h1_h3h3_br",
    "h2_h3h3_br",
    "xs136_lo_h1_pb",
    "xs136_lo_h2_pb",
    "xsec_h2_h1h1_one_h1_invisible_pb",
    "xsec_h1_h2h2_one_h2_invisible_pb",
    "mg5_xsec_gg_heta0_pb",
    "mg5_xsec_pp_eta0Z_pb",
    "mono_higgs_xsec_pb",
    "mono_z_xsec_pb",
)

OBSERVED_PARAMETER_COLUMNS = (
    ("M2", "GeV"),
    ("M3", "GeV"),
    ("vs", "GeV"),
    ("vx", "GeV"),
    ("a12", "rad"),
    ("a13", "rad"),
    ("a23", "rad"),
    ("k1", ""),
    ("k2", ""),
    ("k3", ""),
    ("lX", ""),
    ("lPhiX", ""),
    ("lSX", ""),
    ("K133", "GeV"),
    ("K233", "GeV"),
)


class PlotUnavailable(RuntimeError):
    """Raised when a requested observable has no finite plottable values."""


@dataclass(frozen=True)
class ScanData:
    source: Path
    columns: tuple[str, ...]
    floats: dict[str, np.ndarray]
    bools: dict[str, np.ndarray]
    strings: dict[str, np.ndarray]
    derived: dict[str, np.ndarray]
    metadata: dict[str, object] | None = None
    metadata_source: Path | None = None
    metadata_error: str | None = None
    signal_error: str | None = None

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

    def s(self, name: str) -> np.ndarray:
        if name in self.derived:
            return self.derived[name]
        return self.strings[name]


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
    selection: str | None = None
    requires_bsmpt: bool = False
    required_columns: tuple[str, ...] = ()
    requires_signal: bool = False


@dataclass(frozen=True)
class YR4CrossSectionGrid:
    mass_gev: np.ndarray
    ggf_pb: np.ndarray
    ggf_scale_up_pct: np.ndarray
    ggf_scale_down_pct: np.ndarray
    ggf_pdfalphas_pct: np.ndarray
    vbf_pb: np.ndarray
    vbf_scale_up_pct: np.ndarray
    vbf_scale_down_pct: np.ndarray
    vbf_pdfalphas_pct: np.ndarray


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

CUMULATIVE_CONSTRAINT_STYLES = OrderedDict(
    [
        (
            "all",
            CategoryStyle("All stored", "#BDBDBD", "o", 7.0, 0.16, 1.0),
        ),
        (
            "hb",
            CategoryStyle("HB", "#CC79A7", "^", 18.0, 0.54, 2.0),
        ),
        (
            "hb_hs",
            CategoryStyle("HB + HS", "#E69F00", "s", 22.0, 0.66, 3.0),
        ),
        (
            "hb_hs_wmass",
            CategoryStyle(
                r"HB + HS + $M_W$", "#0072B2", "D", 26.0, 0.78, 4.0
            ),
        ),
        (
            "hb_hs_wmass_dm",
            CategoryStyle(
                r"HB + HS + $M_W$ + DM",
                "#009E73",
                "*",
                54.0,
                0.96,
                5.0,
                "#202020",
                0.35,
            ),
        ),
    ]
)

BSMPT_STATUS_STYLES = OrderedDict(
    [
        (
            "not run",
            CategoryStyle("Not run", "#BDBDBD", "o", 7.0, 0.14, 1.0),
        ),
        (
            "failed",
            CategoryStyle("BSMPT failed", "#D55E00", "X", 34.0, 0.9, 4.0),
        ),
        (
            "success / no selected FOPT",
            CategoryStyle(
                "Success / no selected FOPT",
                "#0072B2",
                "D",
                24.0,
                0.78,
                2.0,
                "#202020",
                0.25,
            ),
        ),
        (
            "selected weak FOPT",
            CategoryStyle(
                r"Selected FOPT: $v_{\rm EW,true}/T<1$",
                "#E69F00",
                "^",
                30.0,
                0.86,
                3.0,
                "#202020",
                0.25,
            ),
        ),
        (
            "selected strong FOPT",
            CategoryStyle(
                r"Selected FOPT: $v_{\rm EW,true}/T\geq1$",
                "#009E73",
                "*",
                54.0,
                0.96,
                5.0,
                "#202020",
                0.35,
            ),
        ),
    ]
)

BSMPT_PHASE_STYLES = OrderedDict(
    [
        (
            "not run",
            CategoryStyle("Not run", "#BDBDBD", "o", 7.0, 0.12, 1.0),
        ),
        (
            "failed",
            CategoryStyle("BSMPT failed", "#D55E00", "X", 34.0, 0.9, 6.0),
        ),
        (
            "phase unavailable",
            CategoryStyle(
                "Phase history unavailable", "#4D4D4D", "P", 30.0, 0.82, 2.0
            ),
        ),
        (
            "direct EW",
            CategoryStyle(
                "Direct EW entry", "#0072B2", "D", 25.0, 0.82, 3.0
            ),
        ),
        (
            "singlet-assisted",
            CategoryStyle(
                r"Path through $S$-broken phase",
                "#E69F00",
                "s",
                28.0,
                0.86,
                4.0,
            ),
        ),
        (
            "X-broken",
            CategoryStyle(
                r"Path through $X$-broken phase",
                "#CC79A7",
                "X",
                38.0,
                0.9,
                5.0,
                "#202020",
                0.25,
            ),
        ),
        (
            "other / multistep",
            CategoryStyle(
                "Other / multistep path",
                "#009E73",
                "^",
                30.0,
                0.86,
                4.5,
            ),
        ),
    ]
)

BSMPT_STEP_STYLES = OrderedDict(
    [
        (
            "not run",
            CategoryStyle("Not run", "#BDBDBD", "o", 7.0, 0.12, 1.0),
        ),
        (
            "failed",
            CategoryStyle("BSMPT failed", "#D55E00", "X", 34.0, 0.9, 6.0),
        ),
        (
            "unavailable",
            CategoryStyle("EW entry unavailable", "#4D4D4D", "P", 30.0, 0.82, 2.0),
        ),
        (
            "step 0",
            CategoryStyle("Direct EW entry (step 0)", "#0072B2", "D", 25.0, 0.82, 3.0),
        ),
        (
            "step 1",
            CategoryStyle("One prior phase (step 1)", "#E69F00", "s", 28.0, 0.86, 4.0),
        ),
        (
            "step 2+",
            CategoryStyle(
                "Two or more prior phases", "#009E73", "*", 50.0, 0.94, 5.0
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

SIGNAL_WIDTH_STYLES = OrderedDict(
    [
        (
            "narrow",
            CategoryStyle(
                r"$\Gamma_2/M_2<1\%$",
                "#009E73",
                "o",
                28.0,
                0.88,
                2.0,
                "#202020",
                0.25,
            ),
        ),
        (
            "intermediate",
            CategoryStyle(
                r"$1\%\leq\Gamma_2/M_2<10\%$",
                "#E69F00",
                "^",
                38.0,
                0.9,
                3.0,
                "#202020",
                0.3,
            ),
        ),
        (
            "broad",
            CategoryStyle(
                r"$\Gamma_2/M_2\geq10\%$ (NWA unreliable)",
                "#D55E00",
                "X",
                46.0,
                0.94,
                4.0,
                "#202020",
                0.35,
            ),
        ),
    ]
)


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
        r"Portal trilinear magnitude $|K_{233}|$",
        "continuous_mass",
        value="abs_K233",
        norm_kind="log_to_one",
        cmap="viridis",
        colorbar_label=r"$|K_{233}|$ [GeV]",
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
    PlotSpec(
        "23_h2_width_dm_status_m2_m3",
        r"Dark-matter status and $h_2$ total width",
        "continuous_mass",
        scheme="dm",
        value="w2",
        norm_kind="log",
        cmap="viridis",
        colorbar_label=r"$\Gamma(h_2)$ [GeV]",
    ),
    PlotSpec(
        "24_cumulative_constraints_m2_m3",
        r"Cumulative constraint survival: $M_2$ versus $M_3$",
        "cumulative_xy",
        x="M2",
        y="M3",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$M_3$ [GeV]",
    ),
    PlotSpec(
        "25_cumulative_constraints_m2_vs",
        r"Cumulative constraint survival: $M_2$ versus $v_s$",
        "cumulative_xy",
        x="M2",
        y="vs",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$v_s$ [GeV]",
    ),
    PlotSpec(
        "26_cumulative_constraints_m2_a12",
        r"Cumulative constraint survival: $M_2$ versus $a_{12}$",
        "cumulative_xy",
        x="M2",
        y="a12",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$a_{12}$",
    ),
    PlotSpec(
        "27_cumulative_constraints_m3_lx",
        r"Cumulative constraint survival: $M_3$ versus $\lambda_X$",
        "cumulative_xy",
        x="M3",
        y="lX",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\lambda_X$",
    ),
    PlotSpec(
        "28_cumulative_constraints_m3_lphix",
        r"Cumulative constraint survival: $M_3$ versus $\lambda_{\Phi X}$",
        "cumulative_xy",
        x="M3",
        y="lPhiX",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\lambda_{\Phi X}$",
    ),
    PlotSpec(
        "29_cumulative_constraints_m3_lsx",
        r"Cumulative constraint survival: $M_3$ versus $\lambda_{SX}$",
        "cumulative_xy",
        x="M3",
        y="lSX",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\lambda_{SX}$",
    ),
    PlotSpec(
        "30_bsmpt_status_m2_m3",
        "BSMPT evaluation and selected FOPT status",
        "categorical_mass",
        scheme="bsmpt_status",
        requires_bsmpt=True,
    ),
    PlotSpec(
        "31_bsmpt_ew_true_over_t_m2_m3",
        r"Selected BSMPT electroweak order parameter",
        "continuous_mass",
        scheme="bsmpt_phase",
        value="ewpt_ew_true_over_T",
        norm_kind="threshold1",
        cmap="RdBu_r",
        colorbar_label=r"$v_{\rm EW,true}(T_*)/T_*$",
        requires_bsmpt=True,
    ),
    PlotSpec(
        "32_bsmpt_phase_history_m2_m3",
        "BSMPT global phase history",
        "categorical_mass",
        scheme="bsmpt_phase",
        requires_bsmpt=True,
    ),
    PlotSpec(
        "33_bsmpt_ew_entry_step_m2_m3",
        "BSMPT electroweak-entry step",
        "categorical_mass",
        scheme="bsmpt_step",
        requires_bsmpt=True,
    ),
    PlotSpec(
        "34_bsmpt_strength_vs_m2",
        r"Selected BSMPT order parameter versus $M_2$",
        "bsmpt_strength_xy",
        scheme="bsmpt_phase",
        x="M2",
        y="ewpt_ew_true_over_T",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$v_{\rm EW,true}(T_*)/T_*$",
        requires_bsmpt=True,
    ),
    PlotSpec(
        "35_bsmpt_strength_vs_m3",
        r"Selected BSMPT order parameter versus $M_3$",
        "bsmpt_strength_xy",
        scheme="bsmpt_phase",
        x="M3",
        y="ewpt_ew_true_over_T",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$v_{\rm EW,true}(T_*)/T_*$",
        requires_bsmpt=True,
    ),
    PlotSpec(
        "36_bsmpt_counts",
        "BSMPT evaluation and phase-transition counts",
        "bsmpt_bars",
        requires_bsmpt=True,
    ),
    PlotSpec(
        "37_k133_vs_m3_all_resonance",
        r"$K_{133}$ versus $M_3$: all stored points",
        "resonance_xy",
        x="M3",
        y="K133",
        value="m2_minus_2m3",
        norm_kind="signed",
        cmap="trsm_resonance",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$K_{133}$ [GeV]",
        colorbar_label=r"$M_2-2M_3$ [GeV]",
        selection="all",
    ),
    PlotSpec(
        "38_k133_vs_m3_experimental_resonance",
        r"$K_{133}$ versus $M_3$: experimental pass",
        "resonance_xy",
        x="M3",
        y="K133",
        value="m2_minus_2m3",
        norm_kind="signed",
        cmap="trsm_resonance",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$K_{133}$ [GeV]",
        colorbar_label=r"$M_2-2M_3$ [GeV]",
        selection="experimental",
    ),
    PlotSpec(
        "39_k133_vs_m3_relic_pass_resonance",
        r"$K_{133}$ versus $M_3$: relic-density pass",
        "resonance_xy",
        x="M3",
        y="K133",
        value="m2_minus_2m3",
        norm_kind="signed",
        cmap="trsm_resonance",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$K_{133}$ [GeV]",
        colorbar_label=r"$M_2-2M_3$ [GeV]",
        selection="relic_pass",
    ),
    PlotSpec(
        "40_k233_vs_m3_all_resonance",
        r"$K_{233}$ versus $M_3$: all stored points",
        "resonance_xy",
        x="M3",
        y="K233",
        value="m2_minus_2m3",
        norm_kind="signed",
        cmap="trsm_resonance",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$K_{233}$ [GeV]",
        colorbar_label=r"$M_2-2M_3$ [GeV]",
        selection="all",
    ),
    PlotSpec(
        "41_k233_vs_m3_experimental_resonance",
        r"$K_{233}$ versus $M_3$: experimental pass",
        "resonance_xy",
        x="M3",
        y="K233",
        value="m2_minus_2m3",
        norm_kind="signed",
        cmap="trsm_resonance",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$K_{233}$ [GeV]",
        colorbar_label=r"$M_2-2M_3$ [GeV]",
        selection="experimental",
    ),
    PlotSpec(
        "42_k233_vs_m3_relic_pass_resonance",
        r"$K_{233}$ versus $M_3$: relic-density pass",
        "resonance_xy",
        x="M3",
        y="K233",
        value="m2_minus_2m3",
        norm_kind="signed",
        cmap="trsm_resonance",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$K_{233}$ [GeV]",
        colorbar_label=r"$M_2-2M_3$ [GeV]",
        selection="relic_pass",
    ),
    PlotSpec(
        "43_k133_experimental_m2_m3",
        r"$K_{133}$ for experimentally allowed points",
        "selected_continuous_mass",
        value="K133",
        norm_kind="signed",
        cmap="coolwarm",
        colorbar_label=r"$K_{133}$ [GeV]",
        selection="experimental",
    ),
    PlotSpec(
        "44_k133_dm_m2_m3",
        r"$K_{133}$ for aggregate-DM-passing points",
        "selected_continuous_mass",
        value="K133",
        norm_kind="signed",
        cmap="coolwarm",
        colorbar_label=r"$K_{133}$ [GeV]",
        selection="dm",
    ),
    PlotSpec(
        "45_k233_experimental_m2_m3",
        r"$K_{233}$ for experimentally allowed points",
        "selected_continuous_mass",
        value="K233",
        norm_kind="signed",
        cmap="coolwarm",
        colorbar_label=r"$K_{233}$ [GeV]",
        selection="experimental",
    ),
    PlotSpec(
        "46_k233_dm_m2_m3",
        r"$K_{233}$ for aggregate-DM-passing points",
        "selected_continuous_mass",
        value="K233",
        norm_kind="signed",
        cmap="coolwarm",
        colorbar_label=r"$K_{233}$ [GeV]",
        selection="dm",
    ),
    PlotSpec(
        "47_m2_vs_a12",
        r"Signed $H_1$--$H_2$ mixing angle versus $M_2$",
        "categorical_xy",
        scheme="fourway",
        x="M2",
        y="a12",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$a_{12}$ [rad]",
    ),
    PlotSpec(
        "48_h2_h1h1_one_h1_invisible_xsec_vs_m2",
        r"$H_2\to H_1H_1$ with exactly one invisible $H_1$",
        "rate_xy",
        x="M2",
        y="xsec_h2_h1h1_one_h1_invisible_pb",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$\sigma_{\rm excl}$ [pb]",
        selection="full_viability",
        required_columns=("xsec_h2_h1h1_one_h1_invisible_pb",),
    ),
    PlotSpec(
        "49_h2_h1h1_one_h1_invisible_xsec_vs_m3",
        r"$H_2\to H_1H_1$ with exactly one invisible $H_1$",
        "rate_xy",
        x="M3",
        y="xsec_h2_h1h1_one_h1_invisible_pb",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\sigma_{\rm excl}$ [pb]",
        selection="full_viability",
        required_columns=("xsec_h2_h1h1_one_h1_invisible_pb",),
    ),
    PlotSpec(
        "50_h1_h2h2_one_h2_invisible_xsec_vs_m2",
        r"$H_1\to H_2H_2$ with exactly one invisible $H_2$",
        "rate_xy",
        x="M2",
        y="xsec_h1_h2h2_one_h2_invisible_pb",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$\sigma_{\rm excl}$ [pb]",
        selection="full_viability",
        required_columns=("xsec_h1_h2h2_one_h2_invisible_pb",),
    ),
    PlotSpec(
        "51_h1_h2h2_one_h2_invisible_xsec_vs_m3",
        r"$H_1\to H_2H_2$ with exactly one invisible $H_2$",
        "rate_xy",
        x="M3",
        y="xsec_h1_h2h2_one_h2_invisible_pb",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\sigma_{\rm excl}$ [pb]",
        selection="full_viability",
        required_columns=("xsec_h1_h2h2_one_h2_invisible_pb",),
    ),
    PlotSpec(
        "52_mono_higgs_xsec_vs_m2",
        r"Mono-Higgs: $gg\to H_1H_2$, $H_2\to H_3H_3$",
        "rate_xy",
        x="M2",
        y="mono_higgs_xsec_pb",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$\sigma(gg\to H_1H_2)\,\mathrm{BR}(H_2\to H_3H_3)$ [pb]",
        selection="full_viability",
        required_columns=("mono_higgs_xsec_pb",),
    ),
    PlotSpec(
        "53_mono_higgs_xsec_vs_m3",
        r"Mono-Higgs: $gg\to H_1H_2$, $H_2\to H_3H_3$",
        "rate_xy",
        x="M3",
        y="mono_higgs_xsec_pb",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\sigma(gg\to H_1H_2)\,\mathrm{BR}(H_2\to H_3H_3)$ [pb]",
        selection="full_viability",
        required_columns=("mono_higgs_xsec_pb",),
    ),
    PlotSpec(
        "54_mono_z_xsec_vs_m2",
        r"Mono-$Z$: $pp\to H_2Z$, $H_2\to H_3H_3$",
        "rate_xy",
        x="M2",
        y="mono_z_xsec_pb",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$\sigma(pp\to H_2Z)\,\mathrm{BR}(H_2\to H_3H_3)$ [pb]",
        selection="full_viability",
        required_columns=("mono_z_xsec_pb",),
    ),
    PlotSpec(
        "55_mono_z_xsec_vs_m3",
        r"Mono-$Z$: $pp\to H_2Z$, $H_2\to H_3H_3$",
        "rate_xy",
        x="M3",
        y="mono_z_xsec_pb",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\sigma(pp\to H_2Z)\,\mathrm{BR}(H_2\to H_3H_3)$ [pb]",
        selection="full_viability",
        required_columns=("mono_z_xsec_pb",),
    ),
    PlotSpec(
        "56_h2_h1h1_one_h1_invisible_xsec_no_dm_vs_m2",
        r"$H_2\to H_1H_1$ with exactly one invisible $H_1$",
        "rate_xy",
        x="M2",
        y="xsec_h2_h1h1_one_h1_invisible_pb",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$\sigma_{\rm excl}$ [pb]",
        selection="non_dm_viability",
        required_columns=("xsec_h2_h1h1_one_h1_invisible_pb",),
    ),
    PlotSpec(
        "57_h2_h1h1_one_h1_invisible_xsec_no_dm_vs_m3",
        r"$H_2\to H_1H_1$ with exactly one invisible $H_1$",
        "rate_xy",
        x="M3",
        y="xsec_h2_h1h1_one_h1_invisible_pb",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\sigma_{\rm excl}$ [pb]",
        selection="non_dm_viability",
        required_columns=("xsec_h2_h1h1_one_h1_invisible_pb",),
    ),
    PlotSpec(
        "58_h1_h2h2_one_h2_invisible_xsec_no_dm_vs_m2",
        r"$H_1\to H_2H_2$ with exactly one invisible $H_2$",
        "rate_xy",
        x="M2",
        y="xsec_h1_h2h2_one_h2_invisible_pb",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$\sigma_{\rm excl}$ [pb]",
        selection="non_dm_viability",
        required_columns=("xsec_h1_h2h2_one_h2_invisible_pb",),
    ),
    PlotSpec(
        "59_h1_h2h2_one_h2_invisible_xsec_no_dm_vs_m3",
        r"$H_1\to H_2H_2$ with exactly one invisible $H_2$",
        "rate_xy",
        x="M3",
        y="xsec_h1_h2h2_one_h2_invisible_pb",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\sigma_{\rm excl}$ [pb]",
        selection="non_dm_viability",
        required_columns=("xsec_h1_h2h2_one_h2_invisible_pb",),
    ),
    PlotSpec(
        "60_mono_higgs_xsec_no_dm_vs_m2",
        r"Mono-Higgs: $gg\to H_1H_2$, $H_2\to H_3H_3$",
        "rate_xy",
        x="M2",
        y="mono_higgs_xsec_pb",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$\sigma(gg\to H_1H_2)\,\mathrm{BR}(H_2\to H_3H_3)$ [pb]",
        selection="non_dm_viability",
        required_columns=("mono_higgs_xsec_pb",),
    ),
    PlotSpec(
        "61_mono_higgs_xsec_no_dm_vs_m3",
        r"Mono-Higgs: $gg\to H_1H_2$, $H_2\to H_3H_3$",
        "rate_xy",
        x="M3",
        y="mono_higgs_xsec_pb",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\sigma(gg\to H_1H_2)\,\mathrm{BR}(H_2\to H_3H_3)$ [pb]",
        selection="non_dm_viability",
        required_columns=("mono_higgs_xsec_pb",),
    ),
    PlotSpec(
        "62_mono_z_xsec_no_dm_vs_m2",
        r"Mono-$Z$: $pp\to H_2Z$, $H_2\to H_3H_3$",
        "rate_xy",
        x="M2",
        y="mono_z_xsec_pb",
        xlabel=r"$M_2$ [GeV]",
        ylabel=r"$\sigma(pp\to H_2Z)\,\mathrm{BR}(H_2\to H_3H_3)$ [pb]",
        selection="non_dm_viability",
        required_columns=("mono_z_xsec_pb",),
    ),
    PlotSpec(
        "63_mono_z_xsec_no_dm_vs_m3",
        r"Mono-$Z$: $pp\to H_2Z$, $H_2\to H_3H_3$",
        "rate_xy",
        x="M3",
        y="mono_z_xsec_pb",
        xlabel=r"$M_3$ [GeV]",
        ylabel=r"$\sigma(pp\to H_2Z)\,\mathrm{BR}(H_2\to H_3H_3)$ [pb]",
        selection="non_dm_viability",
        required_columns=("mono_z_xsec_pb",),
    ),
    PlotSpec(
        "64_signal_rate_m2_m3",
        r"Full-viable $h_2\to h_3h_3$ signal rate",
        "signal_mass",
        value="signal_dominant_rate_fb",
        norm_kind="log",
        cmap="viridis",
        colorbar_label=(
            r"$[\sigma_{\rm ggF}+\sigma_{\rm VBF}]"
            r"\,\mathrm{BR}(h_2\to h_3h_3)$ [fb]"
        ),
        requires_signal=True,
    ),
    PlotSpec(
        "65_signal_rates_vs_m2",
        r"$h_2\to h_3h_3$ production rates versus $M_2$",
        "signal_rates_m2",
        requires_signal=True,
    ),
    PlotSpec(
        "66_signal_rate_vs_m3",
        r"Dominant invisible rate versus dark-matter mass",
        "signal_rate_m3",
        requires_signal=True,
    ),
    PlotSpec(
        "67_signal_k2sq_vs_br",
        r"Production coupling versus invisible branching fraction",
        "signal_k2sq_br",
        requires_signal=True,
    ),
    PlotSpec(
        "68_signal_width_fraction_vs_rate",
        r"Finite-width diagnostic for the invisible signal",
        "signal_width_rate",
        requires_signal=True,
    ),
    PlotSpec(
        "69_signal_rate_vs_direct_detection",
        r"Collider and dark-matter complementarity",
        "signal_dm_complementarity",
        requires_signal=True,
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
        (
            "dashboard_cumulative_constraint_summary",
            (
                "24_cumulative_constraints_m2_m3",
                "25_cumulative_constraints_m2_vs",
                "26_cumulative_constraints_m2_a12",
                "27_cumulative_constraints_m3_lx",
                "28_cumulative_constraints_m3_lphix",
                "29_cumulative_constraints_m3_lsx",
            ),
        ),
        (
            "dashboard_bsmpt_summary",
            (
                "30_bsmpt_status_m2_m3",
                "31_bsmpt_ew_true_over_t_m2_m3",
                "32_bsmpt_phase_history_m2_m3",
                "33_bsmpt_ew_entry_step_m2_m3",
                "34_bsmpt_strength_vs_m2",
                "36_bsmpt_counts",
            ),
        ),
        (
            "dashboard_portal_resonance_summary",
            (
                "37_k133_vs_m3_all_resonance",
                "38_k133_vs_m3_experimental_resonance",
                "39_k133_vs_m3_relic_pass_resonance",
                "40_k233_vs_m3_all_resonance",
                "41_k233_vs_m3_experimental_resonance",
                "42_k233_vs_m3_relic_pass_resonance",
            ),
        ),
        (
            "dashboard_portal_mass_plane_summary",
            (
                "43_k133_experimental_m2_m3",
                "44_k133_dm_m2_m3",
                "45_k233_experimental_m2_m3",
                "46_k233_dm_m2_m3",
            ),
        ),
        (
            "dashboard_scalar_cascade_rates",
            (
                "48_h2_h1h1_one_h1_invisible_xsec_vs_m2",
                "49_h2_h1h1_one_h1_invisible_xsec_vs_m3",
                "50_h1_h2h2_one_h2_invisible_xsec_vs_m2",
                "51_h1_h2h2_one_h2_invisible_xsec_vs_m3",
            ),
        ),
        (
            "dashboard_mg5_mono_rates",
            (
                "52_mono_higgs_xsec_vs_m2",
                "53_mono_higgs_xsec_vs_m3",
                "54_mono_z_xsec_vs_m2",
                "55_mono_z_xsec_vs_m3",
            ),
        ),
        (
            "dashboard_scalar_cascade_rates_no_dm",
            (
                "56_h2_h1h1_one_h1_invisible_xsec_no_dm_vs_m2",
                "57_h2_h1h1_one_h1_invisible_xsec_no_dm_vs_m3",
                "58_h1_h2h2_one_h2_invisible_xsec_no_dm_vs_m2",
                "59_h1_h2h2_one_h2_invisible_xsec_no_dm_vs_m3",
            ),
        ),
        (
            "dashboard_mg5_mono_rates_no_dm",
            (
                "60_mono_higgs_xsec_no_dm_vs_m2",
                "61_mono_higgs_xsec_no_dm_vs_m3",
                "62_mono_z_xsec_no_dm_vs_m2",
                "63_mono_z_xsec_no_dm_vs_m3",
            ),
        ),
        (
            "dashboard_signal_summary",
            (
                "64_signal_rate_m2_m3",
                "65_signal_rates_vs_m2",
                "66_signal_rate_vs_m3",
                "67_signal_k2sq_vs_br",
                "68_signal_width_fraction_vs_rate",
                "69_signal_rate_vs_direct_detection",
            ),
        ),
    ]
)

DASHBOARD_TITLES = {
    "dashboard_status_summary": "TRSM constraint-status summary",
    "dashboard_dm_summary": "TRSM dark-matter constraint summary",
    "dashboard_diagnostic_summary": "TRSM diagnostic maps",
    "dashboard_cumulative_constraint_summary": (
        "TRSM cumulative constraint-survival diagnostics"
    ),
    "dashboard_bsmpt_summary": "TRSM BSMPT electroweak phase-transition summary",
    "dashboard_portal_resonance_summary": (
        r"TRSM portal trilinears and the $M_2=2M_3$ resonance"
    ),
    "dashboard_portal_mass_plane_summary": (
        "TRSM portal trilinears on the mass plane"
    ),
    "dashboard_scalar_cascade_rates": (
        "TRSM exclusive one-invisible scalar-cascade rates at 13.6 TeV"
    ),
    "dashboard_mg5_mono_rates": (
        "TRSM MadGraph mono-Higgs and mono-Z rates at 13.6 TeV"
    ),
    "dashboard_scalar_cascade_rates_no_dm": (
        "TRSM scalar-cascade rates without the DM selection at 13.6 TeV"
    ),
    "dashboard_mg5_mono_rates_no_dm": (
        "TRSM MadGraph mono-Higgs and mono-Z rates without the DM selection"
    ),
    "dashboard_signal_summary": (
        r"Full-viable $h_2\to h_3h_3$ collider-signal summary"
    ),
}

BSMPT_DASHBOARDS = frozenset({"dashboard_bsmpt_summary"})
SIGNAL_DASHBOARDS = frozenset({"dashboard_signal_summary"})


def strict_bool(value: str, column: str = "value", row_number: int | None = None) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    location = f" on row {row_number}" if row_number is not None else ""
    raise ValueError(
        f"Invalid boolean for {column}{location}: {value!r}; expected 'True' or 'False'."
    )


def strict_nullable_bool(
    value: str,
    column: str = "value",
    row_number: int | None = None,
) -> bool | None:
    """Parse a component result that can be absent after a provider failure."""
    if value == "nan":
        return None
    return strict_bool(value, column, row_number)


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


def normalize_optional_text(value: object) -> str:
    """Normalize the scan writer's ``nan`` sentinel for optional text fields."""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def default_yr4_signal_table_path() -> Path:
    return Path(__file__).resolve().parent / "datafiles" / YR4_SIGNAL_TABLE_NAME


def load_yr4_cross_section_grid(
    path: Path | str | None = None,
) -> YR4CrossSectionGrid:
    """Load and strictly validate the tracked YR4 13.6 TeV ggF/VBF grid."""
    path = Path(path) if path is not None else default_yr4_signal_table_path()
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(
            (line for line in stream if not line.startswith("#")),
            delimiter="\t",
        )
        if tuple(reader.fieldnames or ()) != YR4_SIGNAL_COLUMNS:
            raise ValueError(
                f"Unexpected YR4 signal-table columns in {path}: "
                f"{reader.fieldnames!r}; expected {YR4_SIGNAL_COLUMNS!r}."
            )
        buffers = {column: [] for column in YR4_SIGNAL_COLUMNS}
        for row_number, row in enumerate(reader, start=2):
            for column in YR4_SIGNAL_COLUMNS:
                try:
                    value = float(row[column])
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Invalid YR4 value for {column} on data row "
                        f"{row_number}: {row.get(column)!r}."
                    ) from exc
                if not math.isfinite(value):
                    raise ValueError(
                        f"Non-finite YR4 value for {column} on data row "
                        f"{row_number}: {value!r}."
                    )
                buffers[column].append(value)

    if not buffers["mass_gev"]:
        raise ValueError(f"YR4 signal table has no data rows: {path}")
    arrays = {
        column: np.asarray(values, dtype=float)
        for column, values in buffers.items()
    }
    mass = arrays["mass_gev"]
    if np.any(np.diff(mass) <= 0.0):
        raise ValueError(f"YR4 mass grid is not strictly increasing: {path}")
    for column in ("ggf_pb", "vbf_pb"):
        if np.any(arrays[column] <= 0.0):
            raise ValueError(f"YR4 {column} values must all be positive: {path}")
    for prefix in ("ggf", "vbf"):
        if np.any(arrays[f"{prefix}_scale_up_pct"] < 0.0):
            raise ValueError(f"YR4 {prefix} scale-up uncertainties must be nonnegative")
        if np.any(arrays[f"{prefix}_scale_down_pct"] > 0.0):
            raise ValueError(f"YR4 {prefix} scale-down uncertainties must be nonpositive")
        if np.any(arrays[f"{prefix}_pdfalphas_pct"] < 0.0):
            raise ValueError(f"YR4 {prefix} PDF+alpha_s uncertainties must be nonnegative")
    return YR4CrossSectionGrid(**arrays)


def interpolate_no_extrapolation(
    x: np.ndarray,
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    *,
    log_y: bool = False,
) -> np.ndarray:
    """Interpolate a one-dimensional grid and return NaN outside its support."""
    x = np.asarray(x, dtype=float)
    grid_x = np.asarray(grid_x, dtype=float)
    grid_y = np.asarray(grid_y, dtype=float)
    result = np.full(x.shape, np.nan, dtype=float)
    valid = (
        np.isfinite(x)
        & (x >= grid_x[0])
        & (x <= grid_x[-1])
    )
    if not np.any(valid):
        return result
    if log_y:
        if np.any(grid_y <= 0.0):
            raise ValueError("Log interpolation requires a positive grid")
        interpolated = np.exp(
            np.interp(x[valid], grid_x, np.log(grid_y))
        )
    else:
        interpolated = np.interp(x[valid], grid_x, grid_y)
    result[valid] = interpolated
    return result


def interpolate_yr4_cross_sections(
    mass_gev: np.ndarray,
    grid: YR4CrossSectionGrid,
) -> dict[str, np.ndarray]:
    """Interpolate YR4 central rates logarithmically and uncertainties linearly."""
    result = {}
    for column in YR4_SIGNAL_COLUMNS[1:]:
        result[column] = interpolate_no_extrapolation(
            mass_gev,
            grid.mass_gev,
            getattr(grid, column),
            log_y=column.endswith("_pb"),
        )
    return result


def derive_signal_observables(
    m2: np.ndarray,
    m3: np.ndarray,
    k2: np.ndarray,
    w2: np.ndarray,
    h2_h3h3_br: np.ndarray,
    full_viability: np.ndarray,
    grid: YR4CrossSectionGrid,
) -> dict[str, np.ndarray]:
    """Build NWA signal-rate proxies from YR4 production cross sections."""
    m2, m3, k2, w2, h2_h3h3_br = np.broadcast_arrays(
        np.asarray(m2, dtype=float),
        np.asarray(m3, dtype=float),
        np.asarray(k2, dtype=float),
        np.asarray(w2, dtype=float),
        np.asarray(h2_h3h3_br, dtype=float),
    )
    full_viability = np.asarray(full_viability, dtype=bool)
    if full_viability.shape != m2.shape:
        raise ValueError("Signal input arrays must have identical shapes")

    interpolated = interpolate_yr4_cross_sections(m2, grid)
    k2_sq = np.square(k2)
    width_fraction = safe_ratio(w2, m2)
    physical_inputs = (
        finite_mask(m2, m3, k2, w2, h2_h3h3_br)
        & (m2 > 0.0)
        & (m3 >= 0.0)
        & (w2 >= 0.0)
        & (h2_h3h3_br >= 0.0)
        & (h2_h3h3_br <= 1.0)
    )
    decay_open = m2 > 2.0 * m3
    grid_available = finite_mask(
        interpolated["ggf_pb"],
        interpolated["vbf_pb"],
    )

    rates = {}
    for mode in ("ggf", "vbf"):
        central = (
            1000.0
            * interpolated[f"{mode}_pb"]
            * k2_sq
            * h2_h3h3_br
        )
        scale_up = interpolated[f"{mode}_scale_up_pct"]
        scale_down = np.abs(interpolated[f"{mode}_scale_down_pct"])
        pdfalphas = interpolated[f"{mode}_pdfalphas_pct"]
        rel_up = np.hypot(scale_up, pdfalphas) / 100.0
        rel_down = np.hypot(scale_down, pdfalphas) / 100.0
        rates[f"signal_{mode}_rate_fb"] = central
        rates[f"signal_{mode}_rate_low_fb"] = np.maximum(
            central * (1.0 - rel_down), 0.0
        )
        rates[f"signal_{mode}_rate_high_fb"] = central * (1.0 + rel_up)
        rates[f"signal_{mode}_relative_up"] = rel_up
        rates[f"signal_{mode}_relative_down"] = rel_down

    dominant = rates["signal_ggf_rate_fb"] + rates["signal_vbf_rate_fb"]
    dominant_low = (
        rates["signal_ggf_rate_low_fb"]
        + rates["signal_vbf_rate_low_fb"]
    )
    dominant_high = (
        rates["signal_ggf_rate_high_fb"]
        + rates["signal_vbf_rate_high_fb"]
    )
    signal_available = (
        physical_inputs
        & decay_open
        & (h2_h3h3_br > 0.0)
        & grid_available
        & np.isfinite(dominant)
        & (dominant > 0.0)
    )
    viable_open = full_viability & signal_available

    width_category = np.full(m2.shape, "unavailable", dtype=object)
    width_valid = np.isfinite(width_fraction) & (width_fraction >= 0.0)
    width_category[width_valid & (width_fraction < 0.01)] = "narrow"
    width_category[
        width_valid
        & (width_fraction >= 0.01)
        & (width_fraction < 0.10)
    ] = "intermediate"
    width_category[width_valid & (width_fraction >= 0.10)] = "broad"

    return {
        **interpolated,
        **rates,
        "k2_sq": k2_sq,
        "h2_width_fraction": width_fraction,
        "h2_h3h3_kinematically_open": decay_open,
        "signal_yr4_grid_available": grid_available,
        "signal_inputs_physical": physical_inputs,
        "signal_rate_available": signal_available,
        "signal_viable_open": viable_open,
        "signal_width_category": width_category,
        "signal_dominant_rate_fb": dominant,
        "signal_dominant_rate_low_fb": dominant_low,
        "signal_dominant_rate_high_fb": dominant_high,
        "signal_raw_hllhc_events": dominant * HL_LHC_LUMINOSITY_FB,
    }


def derive_bsmpt_results(
    status: np.ndarray,
    strength: np.ndarray,
    phase_path: np.ndarray,
    ew_step_index: np.ndarray,
    has_x_broken: np.ndarray,
    has_x_broken_available: np.ndarray,
) -> dict[str, np.ndarray]:
    """Derive stable BSMPT status and phase-history categories.

    ``ewpt_ew_true_over_T`` is the generator's selected FOPT diagnostic, with
    temperature priority nucleation, percolation, completion, then critical.
    A finite value therefore indicates a selected first-order transition
    result, while a successful BSMPT run can legitimately have no such value.
    """
    status = np.asarray(
        [normalize_optional_text(value).lower() for value in status],
        dtype=object,
    )
    strength = np.asarray(strength, dtype=float)
    phase_path = np.asarray(
        [normalize_optional_text(value) for value in phase_path],
        dtype=object,
    )
    ew_step_index = np.asarray(ew_step_index, dtype=float)
    has_x_broken = np.asarray(has_x_broken, dtype=bool)
    has_x_broken_available = np.asarray(has_x_broken_available, dtype=bool)

    shape = strength.shape
    for name, values in (
        ("status", status),
        ("phase_path", phase_path),
        ("ew_step_index", ew_step_index),
        ("has_x_broken", has_x_broken),
        ("has_x_broken_available", has_x_broken_available),
    ):
        if values.shape != shape:
            raise ValueError(f"BSMPT {name} array has shape {values.shape}, expected {shape}")

    status_present = status != ""
    strength_available = np.isfinite(strength)
    phase_available = phase_path != ""
    step_available = np.isfinite(ew_step_index) & (ew_step_index >= 0.0)
    attempted = (
        status_present
        | strength_available
        | phase_available
        | step_available
        | has_x_broken_available
    )
    failed = status_present & (status != "success")
    success = attempted & ~failed
    selected_fopt = success & strength_available
    strong_fopt = selected_fopt & (strength >= BSMPT_STRONG_EWPT_THRESHOLD)
    weak_fopt = selected_fopt & ~strong_fopt

    status_categories = np.full(shape, "not run", dtype=object)
    status_categories[failed] = "failed"
    status_categories[success] = "success / no selected FOPT"
    status_categories[weak_fopt] = "selected weak FOPT"
    status_categories[strong_fopt] = "selected strong FOPT"

    normalized_paths = np.asarray(
        [" -> ".join(part.strip().upper() for part in path.split("->")) for path in phase_path],
        dtype=object,
    )
    path_has_x = np.asarray(
        ["X_BROKEN" in path for path in normalized_paths], dtype=bool
    )
    path_has_s = np.asarray(
        ["SINGLET_S" in path for path in normalized_paths], dtype=bool
    )
    x_broken_path = success & phase_available & (
        path_has_x | (has_x_broken_available & has_x_broken)
    )
    singlet_path = success & phase_available & path_has_s & ~x_broken_path
    direct_ew = success & step_available & (ew_step_index < 0.5)
    other_path = (
        success
        & phase_available
        & ~x_broken_path
        & ~singlet_path
        & ~direct_ew
    )

    phase_categories = np.full(shape, "not run", dtype=object)
    phase_categories[failed] = "failed"
    phase_categories[success] = "phase unavailable"
    phase_categories[direct_ew & phase_available] = "direct EW"
    phase_categories[other_path] = "other / multistep"
    phase_categories[singlet_path] = "singlet-assisted"
    phase_categories[x_broken_path] = "X-broken"

    step_categories = np.full(shape, "not run", dtype=object)
    step_categories[failed] = "failed"
    step_categories[success] = "unavailable"
    step_categories[success & step_available & (ew_step_index < 0.5)] = "step 0"
    step_categories[
        success
        & step_available
        & (ew_step_index >= 0.5)
        & (ew_step_index < 1.5)
    ] = "step 1"
    step_categories[success & step_available & (ew_step_index >= 1.5)] = "step 2+"

    return {
        "bsmpt_attempted": attempted,
        "bsmpt_success": success,
        "bsmpt_failed": failed,
        "bsmpt_selected_fopt": selected_fopt,
        "bsmpt_strong_fopt": strong_fopt,
        "bsmpt_weak_fopt": weak_fopt,
        "bsmpt_phase_available": success & phase_available,
        "bsmpt_x_broken_path": x_broken_path,
        "bsmpt_direct_ew_entry": direct_ew,
        "bsmpt_multistep_ew_entry": success & step_available & (ew_step_index >= 0.5),
        "bsmpt_status": status_categories,
        "bsmpt_phase": phase_categories,
        "bsmpt_step": step_categories,
    }


def has_bsmpt_results(data: ScanData) -> bool:
    return bool(np.any(data.b("bsmpt_attempted")))


def signal_availability_reason(data: ScanData) -> str | None:
    if data.signal_error is not None:
        return data.signal_error
    missing = [column for column in SIGNAL_REQUIRED_COLUMNS if column not in data.floats]
    if missing:
        return (
            "Signal plots require stored "
            + ", ".join(missing)
            + " columns."
        )
    if "signal_viable_open" not in data.derived:
        return "Signal observables could not be derived."
    if not np.any(data.b("signal_viable_open")):
        return (
            "No full-viable point has a positive, YR4-supported "
            "h2 -> h3 h3 signal rate."
        )
    return None


def has_signal_results(data: ScanData) -> bool:
    return signal_availability_reason(data) is None


def four_way_categories(dm: np.ndarray, experimental: np.ndarray) -> np.ndarray:
    dm = np.asarray(dm, dtype=bool)
    experimental = np.asarray(experimental, dtype=bool)
    categories = np.full(dm.shape, "neither", dtype=object)
    categories[dm & ~experimental] = "DM only"
    categories[~dm & experimental] = "experimental only"
    categories[dm & experimental] = "both"
    return categories


def cumulative_constraint_masks(
    data: ScanData,
) -> OrderedDict[str, np.ndarray]:
    """Return the collaborator diagnostic's nested survivor selections.

    This sequence intentionally mirrors the supplied plots and therefore does
    not apply the EWPO flag.  It is kept separate from ``experimental``, whose
    definition remains HB & HS & EWPO & W-mass.
    """
    hb = data.b("hb")
    hb_hs = hb & data.b("hs")
    hb_hs_wmass = hb_hs & data.b("wmass")
    return OrderedDict(
        [
            ("all", np.ones(len(data), dtype=bool)),
            ("hb", hb),
            ("hb_hs", hb_hs),
            ("hb_hs_wmass", hb_hs_wmass),
            ("hb_hs_wmass_dm", hb_hs_wmass & data.b("dm")),
        ]
    )


def dm_failure_categories(
    dm: np.ndarray,
    relic_excluded: np.ndarray,
    direct_excluded: np.ndarray,
    available: np.ndarray | None = None,
) -> np.ndarray:
    dm = np.asarray(dm, dtype=bool)
    relic_excluded = np.asarray(relic_excluded, dtype=bool)
    direct_excluded = np.asarray(direct_excluded, dtype=bool)
    if available is None:
        available = np.ones(dm.shape, dtype=bool)
    available = np.asarray(available, dtype=bool)
    if not (dm.shape == relic_excluded.shape == direct_excluded.shape == available.shape):
        raise ValueError("DM failure-category arrays must have identical shapes")
    categories = np.full(dm.shape, "other DM failure", dtype=object)
    categories[~available] = "DM unavailable"
    categories[available & dm] = "pass"
    failed = available & ~dm
    categories[failed & relic_excluded & ~direct_excluded] = "relic only"
    categories[failed & ~relic_excluded & direct_excluded] = "direct only"
    categories[failed & relic_excluded & direct_excluded] = "relic + direct"
    return categories


def indirect_categories(
    available: np.ndarray,
    excluded: np.ndarray,
    dm_result_available: np.ndarray | None = None,
) -> np.ndarray:
    available = np.asarray(available, dtype=bool)
    excluded = np.asarray(excluded, dtype=bool)
    if dm_result_available is None:
        dm_result_available = np.ones(available.shape, dtype=bool)
    dm_result_available = np.asarray(dm_result_available, dtype=bool)
    if not (available.shape == excluded.shape == dm_result_available.shape):
        raise ValueError("Indirect-category arrays must have identical shapes")
    categories = np.full(available.shape, "unavailable", dtype=object)
    categories[~dm_result_available] = "DM unavailable"
    categories[dm_result_available & available & ~excluded] = "allowed"
    categories[dm_result_available & available & excluded] = "excluded"
    return categories


def default_scan_metadata_path(scan_path: Path | str) -> Path:
    return Path(scan_path).with_suffix(".metadata.json")


def load_scan_metadata(
    scan_path: Path | str,
    metadata_path: Path | str | None = None,
) -> tuple[dict[str, object] | None, Path | None, str | None]:
    """Load optional generator metadata without making legacy scans unusable."""
    explicit = metadata_path is not None
    candidate = (
        Path(metadata_path)
        if explicit
        else default_scan_metadata_path(scan_path)
    )
    if not candidate.exists():
        error = f"Scan metadata file does not exist: {candidate}" if explicit else None
        return None, candidate if explicit else None, error

    try:
        with candidate.open(encoding="utf-8") as stream:
            payload = json.load(stream)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return None, candidate, f"Could not read scan metadata {candidate}: {exc}"

    if not isinstance(payload, dict):
        return None, candidate, f"Scan metadata {candidate} must contain a JSON object."
    if payload.get("schema") != SCAN_METADATA_SCHEMA:
        return (
            None,
            candidate,
            f"Unsupported scan metadata schema in {candidate}: {payload.get('schema')!r}",
        )
    if not isinstance(payload.get("variable_ranges"), list):
        return None, candidate, f"Scan metadata {candidate} has no variable_ranges list."
    return payload, candidate, None


def load_scan(
    path: Path | str,
    metadata_path: Path | str | None = None,
) -> ScanData:
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

        required = (
            set(BOOLEAN_COLUMNS)
            | set(NULLABLE_BOOLEAN_COLUMNS)
            | set(NUMERIC_COLUMNS)
        )
        missing = sorted(required - set(header))
        if missing:
            raise ValueError(f"Missing required input columns: {', '.join(missing)}")

        index = {column: position for position, column in enumerate(header)}
        bool_buffers = {
            column: []
            for column in BOOLEAN_COLUMNS + OPTIONAL_BOOLEAN_COLUMNS
            if column in index
        }
        nullable_columns = list(NULLABLE_BOOLEAN_COLUMNS) + [
            column
            for column in OPTIONAL_NULLABLE_BOOLEAN_COLUMNS
            if column in index
        ]
        nullable_bool_buffers = {column: [] for column in nullable_columns}
        nullable_available_buffers = {
            column: [] for column in nullable_columns
        }
        float_buffers = {
            column: []
            for column in NUMERIC_COLUMNS + OPTIONAL_NUMERIC_COLUMNS
            if column in index
        }
        string_buffers = {
            column: [] for column in OPTIONAL_TEXT_COLUMNS if column in index
        }

        for row_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(
                    f"Row {row_number} has {len(row)} fields; expected {len(header)}."
                )
            for column in bool_buffers:
                bool_buffers[column].append(
                    strict_bool(row[index[column]], column, row_number)
                )
            for column in nullable_bool_buffers:
                value = strict_nullable_bool(
                    row[index[column]], column, row_number
                )
                nullable_bool_buffers[column].append(value is True)
                nullable_available_buffers[column].append(value is not None)
            for column in float_buffers:
                float_buffers[column].append(
                    strict_float(row[index[column]], column, row_number)
                )
            for column in string_buffers:
                string_buffers[column].append(
                    normalize_optional_text(row[index[column]])
                )

    bools = {
        column: np.asarray(values, dtype=bool)
        for column, values in bool_buffers.items()
    }
    bools.update(
        {
            column: np.asarray(values, dtype=bool)
            for column, values in nullable_bool_buffers.items()
        }
    )
    nullable_available = {
        column: np.asarray(values, dtype=bool)
        for column, values in nullable_available_buffers.items()
    }
    floats = {
        column: np.asarray(values, dtype=float)
        for column, values in float_buffers.items()
    }
    strings = {
        column: np.asarray(values, dtype=object)
        for column, values in string_buffers.items()
    }
    if len(floats["M2"]) == 0:
        raise ValueError(f"Input file has a header but no data rows: {path}")
    for column in OPTIONAL_BOOLEAN_COLUMNS:
        if column not in bools:
            # Legacy scans predate explicit invisible-width provenance.  Such
            # rows are intentionally treated as unmodelled, not as passing.
            bools[column] = np.zeros(len(floats["M2"]), dtype=bool)
    for column in OPTIONAL_NULLABLE_BOOLEAN_COLUMNS:
        if column not in bools:
            bools[column] = np.zeros(len(floats["M2"]), dtype=bool)
            nullable_available[column] = np.zeros(
                len(floats["M2"]), dtype=bool
            )
    for column in OPTIONAL_TEXT_COLUMNS:
        if column not in strings:
            strings[column] = np.full(len(floats["M2"]), "", dtype=object)
    if not np.all(finite_mask(floats["M2"], floats["M3"])):
        raise ValueError("M2 and M3 must be finite for every scan row.")

    theory = bools["evo"] & bools["thc"]
    experimental = bools["hb"] & bools["hs"] & bools["ewpo"] & bools["wmass"]
    non_dm_viability = theory & experimental
    full_viability = theory & experimental & bools["dm"]
    relic_available = nullable_available["dm_relic_excluded"]
    direct_available = nullable_available["dm_direct_detection_excluded"]
    dm_result_available = np.logical_and.reduce(
        [nullable_available[column] for column in NULLABLE_BOOLEAN_COLUMNS]
    )
    relic_pass = relic_available & ~bools["dm_relic_excluded"]
    direct_pass = direct_available & ~bools["dm_direct_detection_excluded"]
    relic_ratio = safe_ratio(floats["dm_omega"], floats["dm_relic_upper_limit"])
    direct_ratio = safe_ratio(floats["dm_dir_det"], floats["dm_dir_det_limit"])
    indirect_ratio = np.where(
        dm_result_available
        & bools["dm_indirect_available"]
        & np.isfinite(floats["dm_indirect_ratio"])
        & (floats["dm_indirect_ratio"] > 0.0),
        floats["dm_indirect_ratio"],
        np.nan,
    )
    ewpt_step_index = floats.get(
        "ewpt_ew_step_index",
        np.full(len(floats["M2"]), np.nan, dtype=float),
    )
    bsmpt_results = derive_bsmpt_results(
        strings["ewpt_status"],
        floats["ewpt_ew_true_over_T"],
        strings["ewpt_global_phase_path"],
        ewpt_step_index,
        bools["ewpt_has_x_broken"],
        nullable_available["ewpt_has_x_broken"],
    )

    derived = {
        "theory": theory,
        "experimental": experimental,
        "non_dm_viability": non_dm_viability,
        "full_viability": full_viability,
        "relic_pass": relic_pass,
        "direct_pass": direct_pass,
        "relic_available": relic_available,
        "direct_available": direct_available,
        "dm_result_available": dm_result_available,
        "fourway": four_way_categories(bools["dm"], experimental),
        "dm_failure": dm_failure_categories(
            bools["dm"],
            bools["dm_relic_excluded"],
            bools["dm_direct_detection_excluded"],
            dm_result_available,
        ),
        "indirect": indirect_categories(
            bools["dm_indirect_available"],
            bools["dm_indirect_detection_excluded"],
            dm_result_available,
        ),
        "relic_ratio": relic_ratio,
        "direct_ratio": direct_ratio,
        "indirect_ratio_available": indirect_ratio,
        "log10_relic_ratio": positive_log10(relic_ratio),
        "log10_direct_ratio": positive_log10(direct_ratio),
        "log10_indirect_ratio": positive_log10(indirect_ratio),
        "abs_a12": np.abs(floats["a12"]),
        "abs_K233": np.abs(floats["K233"]),
        "m2_minus_2m3": floats["M2"] - 2.0 * floats["M3"],
    }
    if (
        "mono_higgs_xsec_pb" not in floats
        and "mg5_xsec_gg_heta0_pb" in floats
        and "h2_h3h3_br" in floats
    ):
        derived["mono_higgs_xsec_pb"] = (
            floats["mg5_xsec_gg_heta0_pb"] * floats["h2_h3h3_br"]
        )
    if (
        "mono_z_xsec_pb" not in floats
        and "mg5_xsec_pp_eta0Z_pb" in floats
        and "h2_h3h3_br" in floats
    ):
        derived["mono_z_xsec_pb"] = (
            floats["mg5_xsec_pp_eta0Z_pb"] * floats["h2_h3h3_br"]
        )
    if (
        "xsec_h2_h1h1_one_h1_invisible_pb" not in floats
        and all(
            column in floats
            for column in ("xs136_lo_h2_pb", "h2_h1h1_br", "h1_h3h3_br")
        )
    ):
        br = floats["h1_h3h3_br"]
        derived["xsec_h2_h1h1_one_h1_invisible_pb"] = (
            floats["xs136_lo_h2_pb"]
            * floats["h2_h1h1_br"]
            * 2.0
            * br
            * (1.0 - br)
        )
    if (
        "xsec_h1_h2h2_one_h2_invisible_pb" not in floats
        and all(
            column in floats
            for column in ("xs136_lo_h1_pb", "h1_h2h2_br", "h2_h3h3_br")
        )
    ):
        br = floats["h2_h3h3_br"]
        derived["xsec_h1_h2h2_one_h2_invisible_pb"] = (
            floats["xs136_lo_h1_pb"]
            * floats["h1_h2h2_br"]
            * 2.0
            * br
            * (1.0 - br)
        )
    derived.update(bsmpt_results)
    signal_error = None
    missing_signal_columns = [
        column for column in SIGNAL_REQUIRED_COLUMNS if column not in floats
    ]
    if missing_signal_columns:
        signal_error = (
            "Signal plots require stored "
            + ", ".join(missing_signal_columns)
            + " columns."
        )
    else:
        try:
            signal_grid = load_yr4_cross_section_grid()
            derived.update(
                derive_signal_observables(
                    floats["M2"],
                    floats["M3"],
                    floats["k2"],
                    floats["w2"],
                    floats["h2_h3h3_br"],
                    full_viability,
                    signal_grid,
                )
            )
        except (OSError, ValueError) as exc:
            signal_error = f"Could not load/derive YR4 signal rates: {exc}"

    metadata, metadata_source, metadata_error = load_scan_metadata(
        path, metadata_path
    )
    if metadata_path is not None and metadata_error is not None:
        raise ValueError(metadata_error)

    return ScanData(
        source=path,
        columns=tuple(header),
        floats=floats,
        bools=bools,
        strings=strings,
        derived=derived,
        metadata=metadata,
        metadata_source=metadata_source,
        metadata_error=metadata_error,
        signal_error=signal_error,
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


def robust_log_norm(values: np.ndarray, include_one: bool = False) -> LogNorm:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite) & (finite > 0.0)]
    if finite.size == 0:
        raise PlotUnavailable("no finite positive values for logarithmic normalization")
    vmin = float(np.percentile(finite, 1.0))
    vmax = float(np.percentile(finite, 99.0))
    if include_one:
        vmax = max(1.0, vmax)
    if not vmin < vmax:
        vmin = max(float(np.min(finite)), vmax * 0.1)
    if not vmin < vmax:
        vmin = vmax * 0.1
    return LogNorm(vmin=vmin, vmax=vmax, clip=True)


def robust_log_norm_to_one(values: np.ndarray) -> LogNorm:
    return robust_log_norm(values, include_one=True)


def norm_for(spec: PlotSpec, values: np.ndarray):
    if spec.norm_kind == "threshold0":
        return robust_threshold_norm(values, center=0.0)
    if spec.norm_kind == "threshold1":
        return robust_threshold_norm(
            values, center=BSMPT_STRONG_EWPT_THRESHOLD
        )
    if spec.norm_kind == "threshold4":
        return robust_threshold_norm(values, center=4.0)
    if spec.norm_kind == "positive":
        return robust_linear_norm(values, include_zero=True)
    if spec.norm_kind == "signed":
        return robust_symlog_norm(values)
    if spec.norm_kind == "log":
        return robust_log_norm(values)
    if spec.norm_kind == "log_to_one":
        return robust_log_norm_to_one(values)
    return robust_linear_norm(values)


def plot_colormap(name: str):
    if name == "trsm_resonance":
        return RESONANCE_CMAP
    return plt.get_cmap(name)


SELECTION_LABELS = {
    "all": "All stored points",
    "experimental": (
        r"Experimental pass: HB $\wedge$ HS $\wedge$ EWPO $\wedge$ $M_W$"
    ),
    "relic_pass": "Relic-density pass",
    "dm": "Aggregate DM pass",
    "non_dm_viability": (
        r"All non-DM constraints: evolution $\wedge$ theory $\wedge$ experimental"
    ),
    "full_viability": (
        r"Full viability: evolution $\wedge$ theory $\wedge$ experimental $\wedge$ DM"
    ),
}


def selection_mask(data: ScanData, selection: str | None) -> np.ndarray:
    if selection in {None, "all"}:
        return np.ones(len(data), dtype=bool)
    try:
        return np.asarray(data.b(selection), dtype=bool)
    except KeyError as exc:
        raise ValueError(f"Unknown point selection: {selection}") from exc


def selection_label(selection: str | None) -> str:
    try:
        return SELECTION_LABELS["all" if selection is None else selection]
    except KeyError as exc:
        raise ValueError(f"Unknown point selection: {selection}") from exc


def category_styles(data: ScanData, scheme: str):
    if scheme == "fourway":
        return data.derived["fourway"], FOURWAY_STYLES
    if scheme == "bsmpt_status":
        return data.derived["bsmpt_status"], BSMPT_STATUS_STYLES
    if scheme == "bsmpt_phase":
        return data.derived["bsmpt_phase"], BSMPT_PHASE_STYLES
    if scheme == "bsmpt_step":
        return data.derived["bsmpt_step"], BSMPT_STEP_STYLES
    if scheme == "indirect":
        styles = OrderedDict(
            [
                (
                    "unavailable",
                    CategoryStyle("Indirect unavailable", "#BDBDBD", "o", 8.0, 0.18, 1.0),
                ),
                (
                    "allowed",
                    CategoryStyle("Available / allowed", "#0072B2", "^", 22.0, 0.75, 2.0),
                ),
                (
                    "excluded",
                    CategoryStyle("Available / excluded", "#D55E00", "X", 38.0, 0.9, 3.0),
                ),
                (
                    "DM unavailable",
                    CategoryStyle("DM evaluation unavailable", "#4D4D4D", "P", 34.0, 0.9, 4.0),
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
        if np.any(data.derived["dm_failure"] == "DM unavailable"):
            styles["DM unavailable"] = CategoryStyle(
                "DM evaluation unavailable", "#4D4D4D", "P", 34.0, 0.9, 4.5
            )
        return data.derived["dm_failure"], styles

    if scheme in {"relic_pass", "direct_pass"}:
        available_name = (
            "relic_available" if scheme == "relic_pass" else "direct_available"
        )
        available = data.b(available_name)
        passed = data.b(scheme)
        categories = np.full(passed.shape, "unavailable", dtype=object)
        categories[available] = "fail"
        categories[passed] = "pass"
        styles = OrderedDict(
            [
                (
                    "unavailable",
                    CategoryStyle(
                        "DM evaluation unavailable",
                        "#4D4D4D",
                        "P",
                        34.0,
                        0.9,
                        3.0,
                    ),
                ),
                (
                    "fail",
                    CategoryStyle("Fail", "#BDBDBD", "o", 8.0, 0.18, 1.0),
                ),
                (
                    "pass",
                    CategoryStyle(
                        "Pass", "#0072B2", "D", 24.0, 0.78, 2.0, "#202020", 0.25
                    ),
                ),
            ]
        )
        return categories, styles

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

    m2_eq_2m3_min = max(xmin, 2.0 * ymin)
    m2_eq_2m3_max = min(xmax, 2.0 * ymax)
    if m2_eq_2m3_min < m2_eq_2m3_max:
        guide_x = np.linspace(m2_eq_2m3_min, m2_eq_2m3_max, 200)
        ax.plot(
            guide_x,
            0.5 * guide_x,
            color=M2_EQ_2M3_GUIDE_COLOR,
            linestyle="-.",
            linewidth=1.0,
            alpha=0.82,
            zorder=1.5,
        )

    m3_eq_2m2_min = max(xmin, 0.5 * ymin)
    m3_eq_2m2_max = min(xmax, 0.5 * ymax)
    if m3_eq_2m2_min < m3_eq_2m2_max:
        guide_x = np.linspace(m3_eq_2m2_min, m3_eq_2m2_max, 200)
        ax.plot(
            guide_x,
            2.0 * guide_x,
            color=M3_EQ_2M2_GUIDE_COLOR,
            linestyle="--",
            linewidth=1.0,
            alpha=0.82,
            zorder=1.5,
        )

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
            0.015,
            0.018,
            (
                r"Guides: $M_2=2M_3$ (magenta dash-dot), "
                r"$M_3=2M_2$ (blue dashed)"
                "\n"
                r"Gray: $M_3=M_2+125$ GeV (dashed), "
                r"$M_3=1000$ GeV (dotted); "
                f"{tail:,} points above the latter"
            ),
            transform=ax.transAxes,
            ha="left",
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
    include_empty: bool = True,
) -> list[Line2D]:
    handles = []
    for key, style in styles.items():
        count = int(np.count_nonzero(categories == key))
        if count == 0 and not include_empty:
            continue
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
    elif spec.scheme == "bsmpt_status":
        ax.set_title(
            spec.title
            + "\n"
            + r"Conventional strong-FOPT diagnostic: "
            + r"$v_{\rm EW,true}(T_*)/T_*\geq1$",
            fontsize=8.8 if compact else 11.2,
        )
    elif spec.scheme == "bsmpt_phase":
        ax.set_title(
            spec.title + "\nGlobal-minimum route on cooling",
            fontsize=8.8 if compact else 11.2,
        )
    elif spec.scheme == "bsmpt_step":
        ax.set_title(
            spec.title + "\nStep 0 is the first transition into an EW phase",
            fontsize=8.8 if compact else 11.2,
        )
    else:
        ax.set_title(spec.title, fontsize=9.5 if compact else 12.0)
    handles = category_legend_handles(
        categories[valid],
        styles,
        int(np.count_nonzero(valid)),
        include_empty=not str(spec.scheme).startswith("bsmpt"),
    )
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
    try:
        values = data.f(spec.value)
    except KeyError as exc:
        raise PlotUnavailable(f"{spec.value} column is unavailable") from exc
    m2 = data.f("M2")
    m3 = data.f("M3")
    valid = finite_mask(m2, m3, values)
    if not np.any(valid):
        raise PlotUnavailable(f"{spec.value} has no finite values")

    norm = norm_for(spec, values[valid])
    cmap = plot_colormap(spec.cmap)
    marker_scheme = spec.scheme or "fourway"
    categories, marker_styles = category_styles(data, marker_scheme)
    for key, style in marker_styles.items():
        mask = valid & (categories == key)
        if not np.any(mask):
            continue
        group_values = values[mask]
        if spec.norm_kind in {"threshold0", "threshold1", "threshold4"}:
            center = {
                "threshold0": 0.0,
                "threshold1": BSMPT_STRONG_EWPT_THRESHOLD,
                "threshold4": 4.0,
            }[spec.norm_kind]
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
            edgecolors=style.edgecolor,
            linewidths=style.linewidth,
            rasterized=True,
            zorder=style.zorder,
        )

    style_mass_axis(ax, data)
    title = spec.title
    if spec.norm_kind == "threshold1":
        title += (
            "\n"
            + r"$T_*$ priority: nucleation, percolation, completion, critical"
        )
    ax.set_title(title, fontsize=8.8 if compact else 12.0)
    if spec.scheme is None:
        ax.text(
            0.985,
            0.018,
            f"finite N = {int(np.count_nonzero(valid)):,}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=6.7 if compact else 7.6,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.72,
                "pad": 1.5,
            },
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
    if spec.norm_kind in {"threshold1", "threshold4"}:
        threshold = (
            BSMPT_STRONG_EWPT_THRESHOLD
            if spec.norm_kind == "threshold1"
            else 4.0
        )
        ticks = [tick for tick in colorbar.get_ticks() if norm.vmin <= tick <= norm.vmax]
        colorbar.set_ticks(sorted(set(ticks + [threshold])))
        colorbar.ax.axhline(
            threshold, color="#222222", linewidth=0.8, alpha=0.8
        )
    handles = category_legend_handles(
        categories[valid],
        marker_styles,
        int(np.count_nonzero(valid)),
        neutral_colors=True,
        include_empty=not str(marker_scheme).startswith("bsmpt"),
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


def add_continuous_colorbar(
    fig,
    ax,
    spec: PlotSpec,
    norm,
    cmap,
    compact: bool = False,
) -> None:
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


def render_resonance_xy(
    fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    x = data.f(spec.x)
    y = data.f(spec.y)
    values = data.f(spec.value)
    selected = selection_mask(data, spec.selection)
    reference = finite_mask(x, y, values)
    valid = selected & reference
    if not np.any(valid):
        raise PlotUnavailable(
            f"{selection_label(spec.selection)} has no finite "
            f"{spec.x}, {spec.y}, and {spec.value} values"
        )

    # Fix the axes and color normalization from the complete finite scan so
    # the all/experimental/relic variants remain directly comparable.
    norm = norm_for(spec, values[reference])
    cmap = plot_colormap(spec.cmap)
    # Draw points farthest from resonance first so the M2 = 2*M3 region is
    # not hidden by a denser off-resonance population.
    order = np.argsort(np.abs(values[valid]), kind="stable")[::-1]
    indices = np.flatnonzero(valid)[order]
    ax.scatter(
        x[indices],
        y[indices],
        c=values[indices],
        norm=norm,
        cmap=cmap,
        s=9.0 if compact else 14.0,
        marker="o",
        alpha=0.82,
        edgecolors="none",
        rasterized=True,
        zorder=2.0,
    )

    xspan = max(float(np.max(x[reference]) - np.min(x[reference])), 1.0)
    ax.set_xlim(
        float(np.min(x[reference]) - 0.035 * xspan),
        float(np.max(x[reference]) + 0.035 * xspan),
    )
    y_linthresh, _ = robust_symlog_parameters(y[reference])
    y_max = max(float(np.max(np.abs(y[reference]))), 1.0e-12)
    ax.set_yscale("symlog", linthresh=y_linthresh)
    ax.set_ylim(-1.08 * y_max, 1.08 * y_max)
    ax.axhline(0.0, color="#666666", linewidth=0.7, alpha=0.55, zorder=0.1)
    ax.set_xlabel(spec.xlabel)
    ax.set_ylabel(spec.ylabel)
    ax.set_title(
        spec.title
        + "\n"
        + selection_label(spec.selection)
        + rf"; $M_2-2M_3=0$ is resonant; $N={int(np.count_nonzero(valid)):,}$",
        fontsize=8.4 if compact else 11.2,
    )
    ax.grid(True, alpha=0.18, linewidth=0.6)
    add_continuous_colorbar(fig, ax, spec, norm, cmap, compact=compact)


def render_selected_continuous_mass(
    fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    values = data.f(spec.value)
    m2 = data.f("M2")
    m3 = data.f("M3")
    mass_valid = finite_mask(m2, m3)
    selected = selection_mask(data, spec.selection)
    colored = mass_valid & selected & np.isfinite(values)
    if not np.any(colored):
        raise PlotUnavailable(
            f"{selection_label(spec.selection)} has no finite {spec.value} values"
        )

    ax.scatter(
        m2[mass_valid],
        m3[mass_valid],
        s=5.0 if compact else 8.0,
        c="#BDBDBD",
        marker="o",
        alpha=0.18,
        edgecolors="none",
        rasterized=True,
        zorder=0.7,
    )
    value_reference = mass_valid & np.isfinite(values)
    # Experimental and DM variants of a given coupling use the same scale,
    # fixed by all finite stored rows rather than by the selected subset.
    norm = norm_for(spec, values[value_reference])
    cmap = plot_colormap(spec.cmap)
    # Keep small-magnitude portal couplings visible in dense regions.
    order = np.argsort(np.abs(values[colored]), kind="stable")[::-1]
    indices = np.flatnonzero(colored)[order]
    ax.scatter(
        m2[indices],
        m3[indices],
        c=values[indices],
        norm=norm,
        cmap=cmap,
        s=10.0 if compact else 16.0,
        marker="o",
        alpha=0.9,
        edgecolors="#4D4D4D",
        linewidths=0.08 if compact else 0.12,
        rasterized=True,
        zorder=3.0,
    )

    style_mass_axis(ax, data)
    ax.set_title(
        spec.title
        + "\nColored: "
        + selection_label(spec.selection)
        + "; gray: all stored points",
        fontsize=8.4 if compact else 11.2,
    )
    ax.text(
        0.985,
        0.018,
        (
            f"colored N = {int(np.count_nonzero(colored)):,}; "
            f"stored N = {int(np.count_nonzero(mass_valid)):,}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.7 if compact else 7.6,
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.72,
            "pad": 1.5,
        },
        zorder=10,
    )
    add_continuous_colorbar(fig, ax, spec, norm, cmap, compact=compact)


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


def render_rate_xy(
    ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    """Render a positive cross section versus one scalar mass on a log axis."""
    try:
        x = data.f(spec.x)
        rate = data.f(spec.y)
    except KeyError as exc:
        raise PlotUnavailable(f"{exc.args[0]} column is unavailable") from exc
    selected = selection_mask(data, spec.selection)
    finite = finite_mask(x, rate)
    valid = selected & finite & (rate > 0.0)
    if not np.any(valid):
        raise PlotUnavailable(
            f"{selection_label(spec.selection)} has no finite positive {spec.y} values"
        )

    indices = np.flatnonzero(valid)[np.argsort(rate[valid], kind="stable")]
    color = "#009E73" if spec.selection == "full_viability" else "#0072B2"
    ax.scatter(
        x[indices],
        rate[indices],
        s=13.0 if compact else 21.0,
        c=color,
        marker="o",
        alpha=0.78,
        edgecolors="#202020",
        linewidths=0.18,
        rasterized=True,
        zorder=2.0,
    )
    xspan = max(float(np.max(x[valid]) - np.min(x[valid])), 1.0)
    ax.set_xlim(
        float(np.min(x[valid]) - 0.035 * xspan),
        float(np.max(x[valid]) + 0.035 * xspan),
    )
    ymin = float(np.min(rate[valid]))
    ymax = float(np.max(rate[valid]))
    if math.isclose(ymin, ymax):
        ymin *= 0.5
        ymax *= 2.0
    else:
        ymin /= 1.35
        ymax *= 1.35
    ax.set_yscale("log")
    ax.set_ylim(max(ymin, np.finfo(float).tiny), ymax)
    ax.set_xlabel(spec.xlabel)
    ax.set_ylabel(spec.ylabel)
    ax.set_title(
        spec.title + "\n" + selection_label(spec.selection),
        fontsize=8.4 if compact else 11.2,
    )
    ax.grid(True, which="both", alpha=0.18, linewidth=0.6)
    selected_finite = selected & finite
    zero_count = int(np.count_nonzero(selected_finite & (rate == 0.0)))
    ax.text(
        0.985,
        0.018,
        (
            f"positive-rate N = {int(np.count_nonzero(valid)):,}; "
            f"zero-rate N = {zero_count:,}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.7 if compact else 7.6,
        bbox={
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.76,
            "pad": 1.5,
        },
        zorder=10,
    )


def cumulative_constraint_legend_handles(
    masks: OrderedDict[str, np.ndarray],
    valid: np.ndarray,
) -> list[Line2D]:
    denominator = int(np.count_nonzero(valid))
    handles = []
    for key, style in CUMULATIVE_CONSTRAINT_STYLES.items():
        count = int(np.count_nonzero(valid & masks[key]))
        percent = 100.0 * count / denominator if denominator else 0.0
        handles.append(
            Line2D(
                [0],
                [0],
                linestyle="None",
                marker=style.marker,
                markersize=max(4.0, math.sqrt(style.size)),
                markerfacecolor=style.color,
                markeredgecolor=(
                    style.edgecolor
                    if style.edgecolor != "none"
                    else style.color
                ),
                markeredgewidth=max(style.linewidth, 0.35),
                color="none",
                label=f"{style.display}: {count:,} ({percent:.2f}%)",
            )
        )
    return handles


def render_cumulative_xy(
    ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    try:
        x = data.f(spec.x)
        y = data.f(spec.y)
    except KeyError as exc:
        missing = exc.args[0]
        raise PlotUnavailable(f"{missing} column is unavailable") from exc

    valid = finite_mask(x, y)
    if not np.any(valid):
        raise PlotUnavailable(f"{spec.x} and {spec.y} have no finite pairs")

    masks = cumulative_constraint_masks(data)
    for key, style in CUMULATIVE_CONSTRAINT_STYLES.items():
        mask = valid & masks[key]
        if not np.any(mask):
            continue
        ax.scatter(
            x[mask],
            y[mask],
            s=style.size * (0.78 if compact else 1.0),
            c=style.color,
            marker=style.marker,
            alpha=style.alpha,
            edgecolors=style.edgecolor,
            linewidths=style.linewidth,
            rasterized=True,
            zorder=style.zorder,
        )

    if spec.x == "M2" and spec.y == "M3":
        style_mass_axis(ax, data)
    else:
        y_finite = y[valid]
        if float(np.min(y_finite)) < 0.0 < float(np.max(y_finite)):
            ax.axhline(
                0.0,
                color="#666666",
                linewidth=0.7,
                alpha=0.55,
                zorder=0.1,
            )
        ax.set_xlabel(spec.xlabel)
        ax.set_ylabel(spec.ylabel)
        ax.grid(True, alpha=0.18, linewidth=0.6)

    ax.set_title(
        spec.title + "\n" + r"Nested sequence; EWPO not applied",
        fontsize=8.7 if compact else 11.5,
    )
    handles = cumulative_constraint_legend_handles(masks, valid)
    ax.legend(
        handles=handles,
        loc="best",
        frameon=True,
        framealpha=0.84,
        edgecolor="none",
        fontsize=5.7 if compact else 7.4,
        handletextpad=0.45,
        borderpad=0.42,
    )


def render_bsmpt_strength_xy(
    ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    x = data.f(spec.x)
    strength = data.f(spec.y)
    valid = finite_mask(x, strength) & (strength >= 0.0)
    if not np.any(valid):
        raise PlotUnavailable("no finite selected BSMPT FOPT strengths")

    categories, styles = category_styles(data, spec.scheme)
    for key, style in styles.items():
        mask = valid & (categories == key)
        if not np.any(mask):
            continue
        order = np.argsort(strength[mask], kind="stable")
        indices = np.flatnonzero(mask)[order]
        ax.scatter(
            x[indices],
            strength[indices],
            s=style.size * (0.82 if compact else 1.0),
            c=style.color,
            marker=style.marker,
            alpha=style.alpha,
            edgecolors=style.edgecolor,
            linewidths=style.linewidth,
            rasterized=True,
            zorder=style.zorder,
        )

    positive = strength[valid & (strength > 0.0)]
    if positive.size and np.count_nonzero(strength[valid] <= 0.0) == 0:
        ax.set_yscale("log")
    else:
        linthresh = 0.05
        if positive.size:
            linthresh = max(
                min(float(np.percentile(positive, 5.0)), 0.25),
                1.0e-4,
            )
        ax.set_yscale("symlog", linthresh=linthresh)
    ax.axhline(
        BSMPT_STRONG_EWPT_THRESHOLD,
        color="#222222",
        linestyle="--",
        linewidth=1.0,
        alpha=0.85,
        zorder=0.5,
    )
    ax.set_xlabel(spec.xlabel)
    ax.set_ylabel(spec.ylabel)
    ax.set_title(
        spec.title
        + "\n"
        + r"$T_*$ priority: nucleation, percolation, completion, critical",
        fontsize=8.8 if compact else 11.2,
    )
    ax.grid(True, which="both", alpha=0.18, linewidth=0.6)
    handles = category_legend_handles(
        categories[valid],
        styles,
        int(np.count_nonzero(valid)),
        include_empty=False,
    )
    ax.legend(
        handles=handles,
        loc="best",
        frameon=True,
        framealpha=0.82,
        edgecolor="none",
        fontsize=5.9 if compact else 7.2,
        handletextpad=0.4,
        borderpad=0.4,
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
        (
            "All non-DM constraints",
            int(np.count_nonzero(data.b("non_dm_viability"))),
            "#D55E00",
            "",
        ),
        (
            "DM evaluation available",
            int(np.count_nonzero(data.b("dm_result_available"))),
            "#56B4E9",
            "//",
        ),
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


def bsmpt_bar_metrics(data: ScanData):
    attempted = data.b("bsmpt_attempted")
    success = data.b("bsmpt_success")
    selected = data.b("bsmpt_selected_fopt")
    no_selected = success & ~selected
    return (
        (
            "BSMPT attempted",
            int(np.count_nonzero(attempted)),
            "#56B4E9",
            "//",
        ),
        (
            "BSMPT successful",
            int(np.count_nonzero(success)),
            "#0072B2",
            "",
        ),
        (
            "BSMPT failed",
            int(np.count_nonzero(data.b("bsmpt_failed"))),
            "#D55E00",
            "",
        ),
        (
            "Phase history available",
            int(np.count_nonzero(data.b("bsmpt_phase_available"))),
            "#CC79A7",
            "",
        ),
        (
            "Success / no selected FOPT",
            int(np.count_nonzero(no_selected)),
            "#0072B2",
            "..",
        ),
        (
            r"Selected FOPT: $v/T<1$",
            int(np.count_nonzero(data.b("bsmpt_weak_fopt"))),
            "#E69F00",
            "",
        ),
        (
            r"Selected FOPT: $v/T\geq1$",
            int(np.count_nonzero(data.b("bsmpt_strong_fopt"))),
            "#009E73",
            "",
        ),
        (
            r"Path includes $X$ breaking",
            int(np.count_nonzero(data.b("bsmpt_x_broken_path"))),
            "#CC79A7",
            "xx",
        ),
    )


def render_bsmpt_bars(
    ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    metrics = bsmpt_bar_metrics(data)
    labels = [metric[0] for metric in metrics]
    counts = [metric[1] for metric in metrics]
    colors = [metric[2] for metric in metrics]
    hatches = [metric[3] for metric in metrics]
    positions = np.arange(len(labels))
    bars = ax.barh(
        positions,
        counts,
        color=colors,
        alpha=0.82,
        edgecolor="#333333",
        linewidth=0.25,
    )
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ax.set_yticks(positions, labels)
    ax.invert_yaxis()
    maximum = max(max(counts, default=0), 1)
    ax.set_xlim(0.0, maximum * 1.22)
    ax.set_xlabel(f"Points (stored total N = {len(data):,})")
    ax.set_title(
        spec.title
        + "\n"
        + r"Strong-FOPT diagnostic uses $v_{\rm EW,true}(T_*)/T_*\geq1$",
        fontsize=8.8 if compact else 11.2,
    )
    ax.grid(True, axis="x", alpha=0.18, linewidth=0.6)
    for position, count in zip(positions, counts):
        ax.text(
            count + maximum * 0.018,
            position,
            f"{count:,} ({100.0 * count / len(data):.2f}%)",
            va="center",
            fontsize=6.3 if compact else 8.0,
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


def require_signal_mask(data: ScanData) -> np.ndarray:
    reason = signal_availability_reason(data)
    if reason is not None:
        raise PlotUnavailable(reason)
    return data.b("signal_viable_open")


def signal_width_legend_handles(
    data: ScanData,
    valid: np.ndarray,
) -> list[Line2D]:
    return category_legend_handles(
        data.derived["signal_width_category"][valid],
        SIGNAL_WIDTH_STYLES,
        int(np.count_nonzero(valid)),
        neutral_colors=True,
        include_empty=False,
    )


def render_signal_mass(
    fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    signal = require_signal_mask(data)
    m2 = data.f("M2")
    m3 = data.f("M3")
    rate = data.f("signal_dominant_rate_fb")
    valid = signal & finite_mask(m2, m3, rate) & (rate > 0.0)
    if not np.any(valid):
        raise PlotUnavailable("no finite positive full-viable signal rates")

    full = data.b("full_viability")
    context = finite_mask(m2, m3) & full & ~valid
    if np.any(context):
        ax.scatter(
            m2[context],
            m3[context],
            s=11.0 if compact else 14.0,
            facecolors="none",
            edgecolors="#9E9E9E",
            linewidths=0.45,
            marker="o",
            alpha=0.5,
            rasterized=True,
            zorder=0.8,
        )

    norm = robust_log_norm(rate[valid])
    cmap = plt.get_cmap(spec.cmap)
    categories = data.derived["signal_width_category"]
    for key, style in SIGNAL_WIDTH_STYLES.items():
        mask = valid & (categories == key)
        if not np.any(mask):
            continue
        order = np.argsort(rate[mask], kind="stable")
        indices = np.flatnonzero(mask)[order]
        ax.scatter(
            m2[indices],
            m3[indices],
            c=rate[indices],
            norm=norm,
            cmap=cmap,
            s=style.size * (0.78 if compact else 1.0),
            marker=style.marker,
            alpha=style.alpha,
            edgecolors=style.edgecolor,
            linewidths=style.linewidth,
            rasterized=True,
            zorder=style.zorder,
        )

    style_mass_axis(ax, data)
    ax.set_title(
        spec.title
        + "\n"
        + r"YR4 13.6 TeV NWA: $(\sigma_{\rm ggF}+\sigma_{\rm VBF})"
        + r"\,k_2^2\,\mathrm{BR}_{33}$",
        fontsize=8.4 if compact else 11.2,
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
    colorbar.set_label(spec.colorbar_label, fontsize=7.5 if compact else 9.2)
    colorbar.ax.tick_params(labelsize=6.7 if compact else 8.0)

    handles = []
    if np.any(context):
        handles.append(
            Line2D(
                [0],
                [0],
                linestyle="None",
                marker="o",
                markersize=5.0,
                markerfacecolor="none",
                markeredgecolor="#8F8F8F",
                markeredgewidth=0.6,
                label=(
                    "Full viable, closed/unsupported: "
                    f"{int(np.count_nonzero(context)):,}"
                ),
            )
        )
    handles.extend(signal_width_legend_handles(data, valid))
    ax.legend(
        handles=handles,
        loc="best",
        frameon=True,
        framealpha=0.82,
        edgecolor="none",
        fontsize=5.7 if compact else 6.9,
        handletextpad=0.4,
        borderpad=0.4,
    )


def render_signal_rates_m2(
    ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    signal = require_signal_mask(data)
    m2 = data.f("M2")
    categories = data.derived["signal_width_category"]
    modes = (
        ("ggf", "ggF", "#D55E00"),
        ("vbf", "VBF", "#0072B2"),
    )
    valid_any = np.zeros(len(data), dtype=bool)
    for mode, _label, color in modes:
        central = data.f(f"signal_{mode}_rate_fb")
        low = data.f(f"signal_{mode}_rate_low_fb")
        high = data.f(f"signal_{mode}_rate_high_fb")
        valid = (
            signal
            & finite_mask(m2, central, low, high)
            & (central > 0.0)
            & (low > 0.0)
            & (high > 0.0)
        )
        valid_any |= valid
        if np.any(valid):
            ax.errorbar(
                m2[valid],
                central[valid],
                yerr=np.vstack(
                    (
                        central[valid] - low[valid],
                        high[valid] - central[valid],
                    )
                ),
                fmt="none",
                ecolor=color,
                elinewidth=0.8 if compact else 1.0,
                capsize=1.5 if compact else 2.2,
                capthick=0.7 if compact else 0.9,
                alpha=0.48,
                rasterized=True,
                zorder=0.8,
            )
        for key, style in SIGNAL_WIDTH_STYLES.items():
            mask = valid & (categories == key)
            if not np.any(mask):
                continue
            ax.scatter(
                m2[mask],
                central[mask],
                s=style.size * (0.72 if compact else 0.9),
                c=color,
                marker=style.marker,
                alpha=style.alpha,
                edgecolors=style.edgecolor,
                linewidths=style.linewidth,
                rasterized=True,
                zorder=style.zorder,
            )
    if not np.any(valid_any):
        raise PlotUnavailable("no finite positive ggF/VBF signal rates")

    ax.set_yscale("log")
    ax.set_xlabel(r"$M_2$ [GeV]")
    ax.set_ylabel(r"$\sigma\,\mathrm{BR}(h_2\to h_3h_3)$ [fb]")
    ax.set_title(
        spec.title
        + "\n"
        + r"Vertical intervals: YR4 scale $\oplus$ PDF+$\alpha_s$",
        fontsize=8.4 if compact else 11.2,
    )
    ax.grid(True, which="both", alpha=0.18, linewidth=0.6)
    event_axis = ax.secondary_yaxis(
        "right",
        functions=(
            lambda rate_fb: rate_fb * HL_LHC_LUMINOSITY_FB,
            lambda events: events / HL_LHC_LUMINOSITY_FB,
        ),
    )
    event_axis.set_ylabel(
        r"Raw $h_3h_3$ events at $3\,\mathrm{ab}^{-1}$"
        + "\n(before acceptance)",
        fontsize=7.0 if compact else 8.6,
    )
    event_axis.tick_params(labelsize=6.5 if compact else 7.5)

    mode_handles = [
        Line2D(
            [0],
            [0],
            color=color,
            marker="o",
            linestyle="None",
            markersize=5.0,
            label=label,
        )
        for _mode, label, color in modes
    ]
    width_handles = signal_width_legend_handles(data, valid_any)
    ax.legend(
        handles=mode_handles + width_handles,
        loc="best",
        frameon=True,
        framealpha=0.82,
        edgecolor="none",
        fontsize=5.5 if compact else 6.8,
        handletextpad=0.4,
        borderpad=0.4,
    )


def render_signal_colored_xy(
    fig,
    ax,
    data: ScanData,
    *,
    x: np.ndarray,
    y: np.ndarray,
    color_values: np.ndarray,
    xlabel: str,
    ylabel: str,
    colorbar_label: str,
    title: str,
    cmap_name: str,
    color_log: bool,
    x_log: bool,
    y_log: bool,
    compact: bool,
) -> np.ndarray:
    signal = require_signal_mask(data)
    valid = signal & finite_mask(x, y, color_values)
    if x_log:
        valid &= x > 0.0
    if y_log:
        valid &= y > 0.0
    if color_log:
        valid &= color_values > 0.0
    if not np.any(valid):
        raise PlotUnavailable("signal observable has no finite plottable values")

    norm = (
        robust_log_norm_to_one(color_values[valid])
        if color_log
        else robust_linear_norm(color_values[valid])
    )
    cmap = plt.get_cmap(cmap_name)
    categories = data.derived["signal_width_category"]
    for key, style in SIGNAL_WIDTH_STYLES.items():
        mask = valid & (categories == key)
        if not np.any(mask):
            continue
        order = np.argsort(color_values[mask], kind="stable")
        indices = np.flatnonzero(mask)[order]
        ax.scatter(
            x[indices],
            y[indices],
            c=color_values[indices],
            norm=norm,
            cmap=cmap,
            s=style.size * (0.78 if compact else 1.0),
            marker=style.marker,
            alpha=style.alpha,
            edgecolors=style.edgecolor,
            linewidths=style.linewidth,
            rasterized=True,
            zorder=style.zorder,
        )
    if x_log:
        ax.set_xscale("log")
    if y_log:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=8.4 if compact else 11.2)
    ax.grid(True, which="both", alpha=0.18, linewidth=0.6)

    scalar_mappable = ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(
        scalar_mappable,
        ax=ax,
        extend="both",
        fraction=0.048 if compact else 0.046,
        pad=0.025,
    )
    colorbar.set_label(colorbar_label, fontsize=7.5 if compact else 9.2)
    colorbar.ax.tick_params(labelsize=6.7 if compact else 8.0)
    ax.legend(
        handles=signal_width_legend_handles(data, valid),
        loc="best",
        frameon=True,
        framealpha=0.82,
        edgecolor="none",
        fontsize=5.7 if compact else 6.9,
        handletextpad=0.4,
        borderpad=0.4,
    )
    return valid


def render_signal_rate_m3(
    fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    render_signal_colored_xy(
        fig,
        ax,
        data,
        x=data.f("M3"),
        y=data.f("signal_dominant_rate_fb"),
        color_values=data.f("M2"),
        xlabel=r"$M_3$ [GeV]",
        ylabel=(
            r"$[\sigma_{\rm ggF}+\sigma_{\rm VBF}]"
            r"\,\mathrm{BR}_{33}$ [fb]"
        ),
        colorbar_label=r"$M_2$ [GeV]",
        title=spec.title + "\nYR4 13.6 TeV NWA; full viability only",
        cmap_name="viridis",
        color_log=False,
        x_log=False,
        y_log=True,
        compact=compact,
    )


def render_signal_k2sq_br(
    fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    render_signal_colored_xy(
        fig,
        ax,
        data,
        x=data.f("k2_sq"),
        y=data.f("h2_h3h3_br"),
        color_values=data.f("signal_dominant_rate_fb"),
        xlabel=r"$k_2^2$",
        ylabel=r"$\mathrm{BR}(h_2\to h_3h_3)$",
        colorbar_label=(
            r"$[\sigma_{\rm ggF}+\sigma_{\rm VBF}]"
            r"\,\mathrm{BR}_{33}$ [fb]"
        ),
        title=spec.title + "\nRate is proportional to $k_2^2\\,\\mathrm{BR}_{33}$",
        cmap_name="viridis",
        color_log=True,
        x_log=True,
        y_log=True,
        compact=compact,
    )


def render_signal_width_rate(
    fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    render_signal_colored_xy(
        fig,
        ax,
        data,
        x=data.f("h2_width_fraction"),
        y=data.f("signal_dominant_rate_fb"),
        color_values=data.f("M2"),
        xlabel=r"$\Gamma_2/M_2$",
        ylabel=(
            r"$[\sigma_{\rm ggF}+\sigma_{\rm VBF}]"
            r"\,\mathrm{BR}_{33}$ [fb]"
        ),
        colorbar_label=r"$M_2$ [GeV]",
        title=spec.title + "\nYR4 production input assumes the NWA",
        cmap_name="viridis",
        color_log=False,
        x_log=True,
        y_log=True,
        compact=compact,
    )
    ax.axvline(
        0.01,
        color="#E69F00",
        linestyle="--",
        linewidth=0.9,
        alpha=0.8,
        zorder=0.5,
    )
    ax.axvline(
        0.10,
        color="#D55E00",
        linestyle="-.",
        linewidth=1.0,
        alpha=0.85,
        zorder=0.5,
    )


def render_signal_dm_complementarity(
    fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False
) -> None:
    render_signal_colored_xy(
        fig,
        ax,
        data,
        x=data.f("direct_ratio"),
        y=data.f("signal_dominant_rate_fb"),
        color_values=data.f("relic_ratio"),
        xlabel=r"$\sigma_{\rm SI}/\sigma_{\rm limit}$",
        ylabel=(
            r"$[\sigma_{\rm ggF}+\sigma_{\rm VBF}]"
            r"\,\mathrm{BR}_{33}$ [fb]"
        ),
        colorbar_label=r"$\Omega/\Omega_{\max}$",
        title=spec.title + "\nFull viability only; unit line is the DD limit",
        cmap_name="viridis",
        color_log=True,
        x_log=True,
        y_log=True,
        compact=compact,
    )
    ax.axvline(
        1.0,
        color="#222222",
        linestyle="--",
        linewidth=1.0,
        alpha=0.85,
        zorder=0.5,
    )


def render_spec(fig, ax, data: ScanData, spec: PlotSpec, compact: bool = False) -> None:
    if spec.requires_bsmpt and not has_bsmpt_results(data):
        raise PlotUnavailable("BSMPT was not run for any stored scan row")
    if spec.requires_signal and not has_signal_results(data):
        raise PlotUnavailable(signal_availability_reason(data) or "signal unavailable")
    if spec.kind == "categorical_mass":
        render_categorical_mass(ax, data, spec, compact=compact)
    elif spec.kind == "continuous_mass":
        render_continuous_mass(fig, ax, data, spec, compact=compact)
    elif spec.kind == "selected_continuous_mass":
        render_selected_continuous_mass(fig, ax, data, spec, compact=compact)
    elif spec.kind == "resonance_xy":
        render_resonance_xy(fig, ax, data, spec, compact=compact)
    elif spec.kind == "categorical_xy":
        render_categorical_xy(ax, data, spec, compact=compact)
    elif spec.kind == "rate_xy":
        render_rate_xy(ax, data, spec, compact=compact)
    elif spec.kind == "cumulative_xy":
        render_cumulative_xy(ax, data, spec, compact=compact)
    elif spec.kind == "bsmpt_strength_xy":
        render_bsmpt_strength_xy(ax, data, spec, compact=compact)
    elif spec.kind == "ratio_plane":
        render_ratio_plane(ax, data, spec, compact=compact)
    elif spec.kind == "bars":
        render_bars(ax, data, spec, compact=compact)
    elif spec.kind == "bsmpt_bars":
        render_bsmpt_bars(ax, data, spec, compact=compact)
    elif spec.kind == "signal_mass":
        render_signal_mass(fig, ax, data, spec, compact=compact)
    elif spec.kind == "signal_rates_m2":
        render_signal_rates_m2(ax, data, spec, compact=compact)
    elif spec.kind == "signal_rate_m3":
        render_signal_rate_m3(fig, ax, data, spec, compact=compact)
    elif spec.kind == "signal_k2sq_br":
        render_signal_k2sq_br(fig, ax, data, spec, compact=compact)
    elif spec.kind == "signal_width_rate":
        render_signal_width_rate(fig, ax, data, spec, compact=compact)
    elif spec.kind == "signal_dm_complementarity":
        render_signal_dm_complementarity(fig, ax, data, spec, compact=compact)
    else:
        raise ValueError(f"Unknown plot kind: {spec.kind}")


def extensions_for_format(plot_format: str) -> tuple[str, ...]:
    if plot_format == "both":
        return ("png", "pdf")
    return (plot_format,)


def all_figure_stems() -> tuple[str, ...]:
    return tuple(spec.stem for spec in PLOT_SPECS) + tuple(DASHBOARDS)


def has_observable(data: ScanData, name: str) -> bool:
    return (
        name in data.floats
        or name in data.bools
        or name in data.strings
        or name in data.derived
    )


def spec_unavailable_reason(data: ScanData, spec: PlotSpec) -> str | None:
    if spec.requires_bsmpt and not has_bsmpt_results(data):
        return "BSMPT was not run for any stored scan row"
    if spec.requires_signal and not has_signal_results(data):
        return signal_availability_reason(data) or "Signal results unavailable"
    missing = [
        column for column in spec.required_columns if not has_observable(data, column)
    ]
    if missing:
        return "missing scan column(s): " + ", ".join(missing)
    if spec.kind == "rate_xy":
        try:
            x = data.f(spec.x)
            rate = data.f(spec.y)
        except KeyError as exc:
            return f"missing scan column: {exc.args[0]}"
        valid = (
            selection_mask(data, spec.selection)
            & finite_mask(x, rate)
            & (rate > 0.0)
        )
        if not np.any(valid):
            return (
                f"{selection_label(spec.selection)} has no finite positive "
                f"{spec.y} values"
            )
    return None


def dashboard_available(data: ScanData, plot_stems: Sequence[str]) -> bool:
    return any(
        spec_unavailable_reason(data, PLOT_BY_STEM[plot_stem]) is None
        for plot_stem in plot_stems
    )


def figure_stems_for_data(data: ScanData) -> tuple[str, ...]:
    plot_stems = tuple(
        spec.stem
        for spec in PLOT_SPECS
        if spec_unavailable_reason(data, spec) is None
    )
    dashboard_stems = tuple(
        stem
        for stem, panels in DASHBOARDS.items()
        if dashboard_available(data, panels)
    )
    return plot_stems + dashboard_stems


def expected_figure_paths(
    output_dir: Path | str,
    plot_format: str,
    data: ScanData | None = None,
) -> list[Path]:
    output_dir = Path(output_dir)
    stems = all_figure_stems() if data is None else figure_stems_for_data(data)
    return [
        output_dir / f"{stem}.{extension}"
        for stem in stems
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
    if spec.kind in {"bars", "bsmpt_bars"}:
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
        ("non_dm_viability", "theory & experimental; DM not required"),
        ("dm", "Stored aggregate DM pass"),
        ("full_viability", "theory & experimental & dm"),
    ):
        rows.append(SummaryRow(name, int(np.count_nonzero(data.b(name))), n, note))

    cumulative_masks = cumulative_constraint_masks(data)
    cumulative_notes = {
        "all": "All stored rows in the collaborator diagnostic sequence",
        "hb": "HB",
        "hb_hs": "HB & HS",
        "hb_hs_wmass": "HB & HS & wmass",
        "hb_hs_wmass_dm": "HB & HS & wmass & dm",
    }
    for name, mask in cumulative_masks.items():
        rows.append(
            SummaryRow(
                f"cumulative_selection_{name}",
                int(np.count_nonzero(mask)),
                n,
                (
                    cumulative_notes[name]
                    + "; EWPO intentionally not applied in this diagnostic sequence"
                ),
            )
        )

    dm_result_available = data.b("dm_result_available")
    relic_available = data.b("relic_available")
    direct_available = data.b("direct_available")
    n_dm_available = int(np.count_nonzero(dm_result_available))
    n_relic_available = int(np.count_nonzero(relic_available))
    n_direct_available = int(np.count_nonzero(direct_available))
    rows.extend(
        [
            SummaryRow(
                "dm_result_available",
                n_dm_available,
                n,
                "Detailed microOMEGAs component result is available",
            ),
            SummaryRow(
                "dm_result_unavailable",
                n - n_dm_available,
                n,
                "Detailed microOMEGAs component result is missing; never treated as passing",
            ),
            SummaryRow(
                "relic_pass",
                int(np.count_nonzero(data.b("relic_pass"))),
                n_relic_available,
                "Not relic-density excluded among available relic results",
            ),
            SummaryRow(
                "relic_unavailable",
                n - n_relic_available,
                n,
                "Relic-density result unavailable",
            ),
            SummaryRow(
                "direct_pass",
                int(np.count_nonzero(data.b("direct_pass"))),
                n_direct_available,
                "Not direct-detection excluded among available direct results",
            ),
            SummaryRow(
                "direct_unavailable",
                n - n_direct_available,
                n,
                "Direct-detection result unavailable",
            ),
        ]
    )

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
    if np.any(data.derived["dm_failure"] == "DM unavailable"):
        dm_failure_names.append("DM unavailable")
    for category in dm_failure_names:
        rows.append(
            SummaryRow(
                f"dm_failure_{category.replace(' ', '_').replace('+', 'and')}",
                int(np.count_nonzero(data.derived["dm_failure"] == category)),
                n,
                "Stored DM result split by relic/direct flags",
            )
        )

    for category, metric_name in (
        ("DM unavailable", "dm_unavailable"),
        ("unavailable", "unavailable"),
        ("allowed", "allowed"),
        ("excluded", "excluded"),
    ):
        rows.append(
            SummaryRow(
                f"indirect_{metric_name}",
                int(np.count_nonzero(data.derived["indirect"] == category)),
                n,
                "Indirect-detection availability/status",
            )
        )

    bsmpt_attempted = int(np.count_nonzero(data.b("bsmpt_attempted")))
    bsmpt_success = int(np.count_nonzero(data.b("bsmpt_success")))
    bsmpt_failed = int(np.count_nonzero(data.b("bsmpt_failed")))
    bsmpt_selected = int(np.count_nonzero(data.b("bsmpt_selected_fopt")))
    bsmpt_strong = int(np.count_nonzero(data.b("bsmpt_strong_fopt")))
    bsmpt_phase_available = int(
        np.count_nonzero(data.b("bsmpt_phase_available"))
    )
    bsmpt_step_available = int(
        np.count_nonzero(
            np.isin(data.derived["bsmpt_step"], ["step 0", "step 1", "step 2+"])
        )
    )
    rows.extend(
        [
            SummaryRow(
                "bsmpt_attempted",
                bsmpt_attempted,
                n,
                (
                    "BSMPT run recorded for this stored row"
                    if bsmpt_attempted
                    else "No BSMPT runs recorded; BSMPT plots omitted"
                ),
            ),
            SummaryRow(
                "bsmpt_success",
                bsmpt_success,
                bsmpt_attempted,
                "Successful BSMPT evaluation among attempted rows",
            ),
            SummaryRow(
                "bsmpt_failed",
                bsmpt_failed,
                bsmpt_attempted,
                "Failed BSMPT evaluation among attempted rows",
            ),
            SummaryRow(
                "bsmpt_success_no_selected_fopt",
                bsmpt_success - bsmpt_selected,
                bsmpt_success,
                "Successful evaluation without a finite selected FOPT strength",
            ),
            SummaryRow(
                "bsmpt_selected_fopt",
                bsmpt_selected,
                bsmpt_success,
                (
                    "Finite selected FOPT result; temperature priority is "
                    "nucleation, percolation, completion, then critical"
                ),
            ),
            SummaryRow(
                "bsmpt_selected_weak_fopt",
                int(np.count_nonzero(data.b("bsmpt_weak_fopt"))),
                bsmpt_selected,
                (
                    "Selected FOPT with "
                    f"v_EW,true(T*)/T* < {BSMPT_STRONG_EWPT_THRESHOLD:g}"
                ),
            ),
            SummaryRow(
                "bsmpt_selected_strong_fopt",
                bsmpt_strong,
                bsmpt_selected,
                (
                    "Conventional strong-FOPT diagnostic: "
                    f"v_EW,true(T*)/T* >= {BSMPT_STRONG_EWPT_THRESHOLD:g}; "
                    "not an additional scan constraint"
                ),
            ),
            SummaryRow(
                "bsmpt_phase_history_available",
                bsmpt_phase_available,
                bsmpt_success,
                "MinimaTracer global-minimum phase route is stored",
            ),
            SummaryRow(
                "bsmpt_ew_entry_step_available",
                bsmpt_step_available,
                bsmpt_success,
                "Index of the transition entering an EW-broken phase is stored",
            ),
            SummaryRow(
                "bsmpt_direct_ew_entry",
                int(np.count_nonzero(data.b("bsmpt_direct_ew_entry"))),
                bsmpt_step_available,
                "EW breaking occurs at step 0",
            ),
            SummaryRow(
                "bsmpt_multistep_ew_entry",
                int(np.count_nonzero(data.b("bsmpt_multistep_ew_entry"))),
                bsmpt_step_available,
                "One or more phase transitions precede EW breaking",
            ),
            SummaryRow(
                "bsmpt_x_broken_path",
                int(np.count_nonzero(data.b("bsmpt_x_broken_path"))),
                bsmpt_phase_available,
                (
                    "Global cooling path includes X_BROKEN or EW_X_BROKEN; "
                    "inspect such points because the nominal dark Z2 is broken "
                    "at an intermediate or final phase"
                ),
            ),
            SummaryRow(
                "ewpt_finite",
                bsmpt_selected,
                n,
                "Backward-compatible alias for finite selected BSMPT FOPT strengths",
            ),
        ]
    )
    for category in BSMPT_PHASE_STYLES:
        rows.append(
            SummaryRow(
                f"bsmpt_phase_{category.replace(' ', '_').replace('/', 'or')}",
                int(np.count_nonzero(data.derived["bsmpt_phase"] == category)),
                n,
                "BSMPT/MinimaTracer global phase-route category",
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
    h2_to_h3h3_open = data.f("M2") > 2.0 * data.f("M3")
    m3_above_2m2 = data.f("M3") > 2.0 * data.f("M2")
    rows.extend(
        [
            SummaryRow(
                "h2_to_h3h3_kinematically_open",
                int(np.count_nonzero(h2_to_h3h3_open)),
                n,
                (
                    "Below the M2 = 2*M3 magenta dash-dot guide: "
                    "h2 -> h3 h3 has nonzero phase space"
                ),
            ),
            SummaryRow(
                "m3_above_2m2_reference",
                int(np.count_nonzero(m3_above_2m2)),
                n,
                (
                    "Above the M3 = 2*M2 blue dashed guide: "
                    "reciprocal mass-hierarchy reference, not an h3 decay "
                    "threshold because h3 is Z2-stable"
                ),
            ),
        ]
    )
    invisible_width_open = (
        (2.0 * data.f("M3") < SM_LIKE_HIGGS_MASS_GEV)
        | h2_to_h3h3_open
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

    full_viable = data.b("full_viability")
    n_full_viable = int(np.count_nonzero(full_viable))
    if "signal_viable_open" in data.derived:
        signal_viable = data.b("signal_viable_open")
        n_signal = int(np.count_nonzero(signal_viable))
        full_open = full_viable & data.b("h2_h3h3_kinematically_open")
        supported = full_open & data.b("signal_yr4_grid_available")
        rates = data.f("signal_dominant_rate_fb")[signal_viable]
        rate_note = (
            "YR4 13.6 TeV NWA ggF+VBF proxy; "
            "sigma = k2^2 * sigma_YR4(M2) * BR(h2 -> h3 h3)"
        )
        if rates.size:
            rate_note += (
                f"; range {float(np.min(rates)):.6g}--"
                f"{float(np.max(rates)):.6g} fb, "
                f"median {float(np.median(rates)):.6g} fb"
            )
        rows.extend(
            [
                SummaryRow(
                    "signal_full_viable_h2_to_h3h3_open",
                    int(np.count_nonzero(full_open)),
                    n_full_viable,
                    "Full-viable points above the h2 -> h3 h3 threshold",
                ),
                SummaryRow(
                    "signal_full_viable_yr4_supported",
                    int(np.count_nonzero(supported)),
                    int(np.count_nonzero(full_open)),
                    "Open full-viable points with M2 inside the 10--3000 GeV YR4 grid",
                ),
                SummaryRow(
                    "signal_full_viable_rate_available",
                    n_signal,
                    n_full_viable,
                    rate_note,
                ),
            ]
        )
        for category, label in (
            ("narrow", "Gamma2/M2 < 1%; NWA-compatible diagnostic region"),
            (
                "intermediate",
                "1% <= Gamma2/M2 < 10%; inspect finite-width effects",
            ),
            (
                "broad",
                "Gamma2/M2 >= 10%; YR4 NWA factorization is not reliable",
            ),
        ):
            rows.append(
                SummaryRow(
                    f"signal_width_{category}",
                    int(
                        np.count_nonzero(
                            signal_viable
                            & (
                                data.derived["signal_width_category"]
                                == category
                            )
                        )
                    ),
                    n_signal,
                    label,
                )
            )
        yr4_rows = len(load_yr4_cross_section_grid().mass_gev)
        rows.append(
            SummaryRow(
                "yr4_signal_table_rows",
                yr4_rows,
                yr4_rows,
                (
                    "Tracked LHCHXSWG YR4 BSM 13.6 TeV ggF/VBF grid; "
                    f"source {YR4_SIGNAL_SOURCE_URL}; commit "
                    f"{YR4_SIGNAL_REPOSITORY_COMMIT}; NWA, no EW corrections"
                ),
            )
        )
    else:
        rows.append(
            SummaryRow(
                "signal_full_viable_rate_available",
                0,
                n_full_viable,
                signal_availability_reason(data) or "Signal observables unavailable",
            )
        )

    component_dm = (
        ~data.bools["dm_relic_excluded"]
        & ~data.bools["dm_direct_detection_excluded"]
        & ~data.bools["dm_indirect_detection_excluded"]
    )
    component_mismatch = dm_result_available & (component_dm != data.bools["dm"])
    rows.append(
        SummaryRow(
            "dm_component_mismatch",
            int(np.count_nonzero(component_mismatch)),
            n_dm_available,
            "Among available results, stored dm differs from conjunction of stored component exclusions",
        )
    )
    rows.append(
        SummaryRow(
            "dm_pass_without_component_results",
            int(np.count_nonzero(data.bools["dm"] & ~dm_result_available)),
            n,
            "Aggregate DM pass is stored although detailed component results are unavailable",
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


def observed_parameter_range(
    data: ScanData, column: str
) -> tuple[float, float] | None:
    values = data.floats.get(column)
    if values is None:
        return None
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None
    return float(np.min(finite)), float(np.max(finite))


def format_index_number(value: object) -> str:
    if value is None:
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return str(value)
    return f"{number:.8g}"


def format_index_range(minimum: object, maximum: object) -> str:
    if minimum is None or maximum is None:
        return "—"
    lower = format_index_number(minimum)
    upper = format_index_number(maximum)
    if lower == upper:
        return lower
    return f"{lower} – {upper}"


def scan_information_html(data: ScanData) -> str:
    """Render configured scan provenance or a clearly labelled legacy fallback."""
    escaped = lambda value: html.escape(str(value), quote=True)
    metadata = data.metadata

    if metadata is not None:
        source = data.metadata_source or default_scan_metadata_path(data.source)
        notice = (
            '<div class="notice"><strong>Configured scan metadata</strong> loaded from '
            f'<code>{escaped(source)}</code>. Configured/effective bounds describe the '
            "generator; observed bounds describe only rows retained in this data file and "
            "can be narrower after selections and rounding.</div>"
        )
        range_entries = metadata.get("variable_ranges", [])
    else:
        expected = default_scan_metadata_path(data.source)
        warning = (
            f"{data.metadata_error} " if data.metadata_error else ""
        )
        notice = (
            '<div class="notice warning"><strong>Observed-range fallback.</strong> '
            f"{escaped(warning)}No usable scan-metadata sidecar was found at "
            f'<code>{escaped(expected)}</code>. The ranges below are extrema of stored rows, '
            "not the configured scan bounds; selections and finite sampling can narrow them.</div>"
        )
        range_entries = [
            {
                "variable": column,
                "column": column,
                "configured_min": None,
                "configured_max": None,
                "effective_min": None,
                "effective_max": None,
                "unit": unit,
                "sampling": "Unavailable for legacy input",
                "note": "Observed stored rows only",
            }
            for column, unit in OBSERVED_PARAMETER_COLUMNS
            if observed_parameter_range(data, column) is not None
        ]

    range_rows = []
    for entry in range_entries:
        if not isinstance(entry, dict):
            continue
        variable = str(entry.get("variable", ""))
        column = str(entry.get("column", variable))
        observed = observed_parameter_range(data, column)
        observed_text = (
            format_index_range(*observed) if observed is not None else "—"
        )
        range_rows.append(
            "<tr>"
            f"<td><code>{escaped(variable)}</code></td>"
            f"<td>{escaped(format_index_range(entry.get('configured_min'), entry.get('configured_max')))}</td>"
            f"<td>{escaped(format_index_range(entry.get('effective_min'), entry.get('effective_max')))}</td>"
            f"<td>{escaped(observed_text)}</td>"
            f"<td>{escaped(entry.get('unit', ''))}</td>"
            f"<td>{escaped(entry.get('sampling', ''))}</td>"
            f"<td>{escaped(entry.get('note', ''))}</td>"
            "</tr>"
        )

    facts = []
    if metadata is not None:
        mass_sampling = metadata.get("mass_sampling", {})
        portal_sampling = metadata.get("portal_sampling", {})
        scan_options = metadata.get("options", {})
        if not isinstance(mass_sampling, dict):
            mass_sampling = {}
        if not isinstance(portal_sampling, dict):
            portal_sampling = {}
        if not isinstance(scan_options, dict):
            scan_options = {}
        for label, value in (
            ("Seed", metadata.get("seed")),
            ("Requested points", metadata.get("requested_points")),
            ("Stopping rule", metadata.get("stopping_rule")),
            ("Rows written", metadata.get("output_selection")),
            ("Mass mode", mass_sampling.get("mode")),
            ("Portal mode", portal_sampling.get("mode")),
            ("Created (UTC)", metadata.get("created_utc")),
            ("Portal convention", metadata.get("portal_convention")),
            ("BSMPT requested", scan_options.get("run_ewpt")),
            ("BSMPT high temperature", scan_options.get("ewpt_thigh")),
        ):
            if value is not None:
                facts.append(
                    '<div class="fact"><span>'
                    f"{escaped(label)}</span><strong>{escaped(value)}</strong></div>"
                )
        descriptions = "".join(
            f"<p>{escaped(description)}</p>"
            for description in (
                mass_sampling.get("description"),
                portal_sampling.get("description"),
            )
            if description
        )
    else:
        facts = [
            '<div class="fact"><span>Input rows</span>'
            f"<strong>{len(data):,}</strong></div>"
        ]
        descriptions = ""

    fixed_html = ""
    fixed_parameters = metadata.get("fixed_parameters", []) if metadata else []
    if isinstance(fixed_parameters, list) and fixed_parameters:
        fixed_rows = []
        for entry in fixed_parameters:
            if not isinstance(entry, dict):
                continue
            fixed_rows.append(
                "<tr>"
                f"<td><code>{escaped(entry.get('variable', ''))}</code></td>"
                f"<td>{escaped(format_index_number(entry.get('value')))}</td>"
                f"<td>{escaped(entry.get('unit', ''))}</td>"
                f"<td>{escaped(entry.get('note', ''))}</td>"
                "</tr>"
            )
        fixed_html = (
            "<h3>Fixed parameters</h3>"
            '<div class="table-wrap"><table><thead><tr><th>Variable</th><th>Value</th>'
            "<th>Unit</th><th>Note</th></tr></thead>"
            f"<tbody>{''.join(fixed_rows)}</tbody></table></div>"
        )

    details_html = ""
    if metadata is not None:
        command_line = metadata.get("command_line")
        command_html = ""
        if isinstance(command_line, list):
            command = shlex.join(str(item) for item in command_line)
            command_html = f"<h4>Command line</h4><pre>{escaped(command)}</pre>"
        options = metadata.get("options")
        options_html = ""
        if isinstance(options, dict):
            options_json = json.dumps(options, indent=2, sort_keys=True)
            options_html = f"<h4>Parsed options</h4><pre>{escaped(options_json)}</pre>"
        if command_html or options_html:
            details_html = (
                "<details><summary>Complete generator invocation</summary>"
                f"{command_html}{options_html}</details>"
            )

    return (
        '<section id="scan"><h2>Scan configuration</h2>'
        f"{notice}<div class=\"fact-grid\">{''.join(facts)}</div>{descriptions}"
        "<h3>Variable ranges</h3>"
        '<div class="table-wrap"><table><thead><tr><th>Variable</th>'
        "<th>Configured</th><th>Effective/support</th><th>Observed stored</th>"
        "<th>Unit</th><th>Sampling</th><th>Note</th></tr></thead>"
        f"<tbody>{''.join(range_rows)}</tbody></table></div>"
        f"{fixed_html}{details_html}</section>"
    )


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
        plot_card(stem, DASHBOARD_TITLES[stem])
        for stem in DASHBOARDS
        if stem not in BSMPT_DASHBOARDS
        and stem not in SIGNAL_DASHBOARDS
    )
    standalone_cards = "\n".join(
        plot_card(spec.stem, spec.title)
        for spec in PLOT_SPECS
        if not spec.requires_bsmpt
        and not spec.requires_signal
    )
    bsmpt_attempted = int(np.count_nonzero(data.b("bsmpt_attempted")))
    if bsmpt_attempted:
        bsmpt_cards = [
            plot_card(stem, DASHBOARD_TITLES[stem])
            for stem in DASHBOARDS
            if stem in BSMPT_DASHBOARDS
        ]
        bsmpt_cards.extend(
            plot_card(spec.stem, spec.title)
            for spec in PLOT_SPECS
            if spec.requires_bsmpt
        )
        bsmpt_section = (
            '<section id="bsmpt"><h2>BSMPT electroweak phase-transition plots</h2>'
            "<p>The selected order parameter uses the first available temperature "
            "in the priority nucleation, percolation, completion, then critical. "
            f"The suite records {bsmpt_attempted:,} attempted BSMPT evaluations; "
            r"$v_{\rm EW,true}(T_*)/T_*\geq1$ is shown as a conventional strong-FOPT "
            "diagnostic, not as an additional scan constraint.</p>"
            f'<div class="grid">{"".join(bsmpt_cards)}</div></section>'
        )
    else:
        bsmpt_section = (
            '<section id="bsmpt"><h2>BSMPT electroweak phase-transition plots</h2>'
            '<div class="notice warning"><strong>No stored BSMPT evaluations.</strong> '
            "These plots are generated automatically when at least one row contains "
            "a BSMPT status, selected EWPT strength, phase history, or EW-entry "
            "step.</div></section>"
        )
    signal_cards = [
        plot_card(stem, DASHBOARD_TITLES[stem])
        for stem in DASHBOARDS
        if stem in SIGNAL_DASHBOARDS
    ]
    signal_cards.extend(
        plot_card(spec.stem, spec.title)
        for spec in PLOT_SPECS
        if spec.requires_signal
    )
    signal_reason = signal_availability_reason(data)
    signal_count = (
        int(np.count_nonzero(data.b("signal_viable_open")))
        if "signal_viable_open" in data.derived
        else 0
    )
    if signal_reason is None:
        signal_intro = (
            f"<p>The suite finds {signal_count:,} full-viable points with a "
            "positive <code>h2 -&gt; h3 h3</code> rate inside the 10--3000 GeV "
            "YR4 grid. It evaluates "
            "<code>sigma_P = k2^2 * sigma_P^YR4(M2) * BR(h2 -&gt; h3 h3)</code> "
            "for ggF and VBF at 13.6 TeV. Central cross sections are interpolated "
            "logarithmically in mass without extrapolation; scale and "
            "PDF+alpha_s components are combined in quadrature for the displayed "
            "intervals.</p>"
            "<p>The YR4 BSM inputs use the narrow-width approximation and omit "
            "electroweak corrections. Points with "
            "<code>Gamma2/M2 &gt;= 10%</code> are retained but marked as a region "
            "where this factorized rate is unreliable. Raw 3 ab<sup>-1</sup> "
            "event counts are before acceptance, triggering, reconstruction, "
            "and backgrounds. "
            f'<a href="{escaped(YR4_SIGNAL_SOURCE_URL)}" target="_blank">'
            "Official LHCHXSWG table</a> &middot; "
            f'<a href="{escaped(YR4_SIGNAL_CITATION_URL)}" target="_blank">'
            "YR4 citation</a>.</p>"
        )
    else:
        signal_intro = (
            '<div class="notice warning"><strong>Signal plots unavailable.</strong> '
            f"{escaped(signal_reason)}</div>"
        )
    signal_section = (
        '<section id="signal"><h2>Full-viability collider signal plots</h2>'
        f"{signal_intro}<div class=\"grid\">{''.join(signal_cards)}</div></section>"
    )
    scan_information = scan_information_html(data)
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
    non_dm_viable = int(np.count_nonzero(data.b("non_dm_viability")))
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
.notice {{ margin: 12px 0 18px; padding: 12px 14px; border-left: 4px solid #0072b2; background: #eef8fd; }}
.notice.warning {{ border-left-color: #e69f00; background: #fff8e8; }}
.fact-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 10px; margin: 16px 0; }}
.fact {{ display: grid; gap: 2px; padding: 10px 12px; border: 1px solid var(--line); border-radius: 7px; background: var(--panel); }}
.fact span {{ color: var(--muted); font-size: .82rem; }}
details {{ margin: 18px 0; padding: 12px 14px; border: 1px solid var(--line); border-radius: 7px; }}
summary {{ cursor: pointer; font-weight: 650; }}
pre {{ overflow-x: auto; padding: 12px; background: #f4f6f8; border-radius: 6px; }}
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
<p class="meta">Input: <code>{escaped(data.source)}</code> &middot; {len(data):,} rows &middot; experimental {experimental:,} &middot; non-DM viable {non_dm_viable:,} &middot; DM {dm_pass:,} &middot; full viability {full:,} &middot; signal points {signal_count:,} &middot; BSMPT attempted {bsmpt_attempted:,}</p>
<nav><a href="#scan">Scan configuration</a><a href="#dashboards">Dashboards</a><a href="#signal">Signals</a><a href="#bsmpt">BSMPT/EWPT</a><a href="#standalone">Individual plots</a><a href="#summary">Constraint summary</a><a href="constraint_summary.tsv">Download TSV</a></nav>
</header>
{scan_information}
<section id="dashboards"><h2>Dashboards</h2><div class="grid">{dashboard_cards}</div></section>
{signal_section}
{bsmpt_section}
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
        "--scan-metadata",
        type=Path,
        help=(
            "Optional scan-metadata JSON sidecar. By default the suite looks "
            "for <input-stem>.metadata.json beside the input scan."
        ),
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
    data = load_scan(args.input, args.scan_metadata)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loaded {len(data):,} rows from {data.source}")
    if data.metadata is not None:
        print(f"Loaded scan metadata from {data.metadata_source}")
    elif data.metadata_error is not None:
        print(f"Warning: {data.metadata_error}")
    else:
        print(
            "No scan metadata sidecar found; the index will label ranges as "
            "observed stored-row extrema."
        )
    print(
        "Selections: "
        f"experimental={np.count_nonzero(data.b('experimental')):,} "
        f"non_dm_viable={np.count_nonzero(data.b('non_dm_viability')):,} "
        f"dm={np.count_nonzero(data.b('dm')):,} "
        f"full={np.count_nonzero(data.b('full_viability')):,}"
    )

    paths: list[Path] = []
    skipped: list[tuple[str, str]] = []
    for spec in PLOT_SPECS:
        reason = spec_unavailable_reason(data, spec)
        if reason is not None:
            skipped.append((spec.stem, reason))
            print(f"Skipped {spec.stem}: {reason}")
            continue
        if spec.requires_signal and not has_signal_results(data):
            reason = signal_availability_reason(data) or "Signal results unavailable"
            skipped.append((spec.stem, reason))
            print(f"Skipped {spec.stem}: {reason}")
            continue
        try:
            paths.extend(
                render_standalone(data, spec, output_dir, args.format, args.dpi)
            )
        except PlotUnavailable as exc:
            reason = str(exc)
            skipped.append((spec.stem, reason))
            print(f"Skipped {spec.stem}: {reason}")

    for dashboard_stem, plot_stems in DASHBOARDS.items():
        if not dashboard_available(data, plot_stems):
            reason = "none of the dashboard observables are available"
            skipped.append((dashboard_stem, reason))
            print(f"Skipped {dashboard_stem}: {reason}")
            continue
        if dashboard_stem in SIGNAL_DASHBOARDS and not has_signal_results(data):
            reason = signal_availability_reason(data) or "Signal results unavailable"
            skipped.append((dashboard_stem, reason))
            print(f"Skipped {dashboard_stem}: {reason}")
            continue
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
