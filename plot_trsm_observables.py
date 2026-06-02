#!/usr/bin/env python3

import argparse
import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path


SM_GF = 1.1663787e-5
SM_VEV = math.sqrt(1.0 / math.sqrt(2.0) / SM_GF)
DEFAULT_M1 = 125.09
EWPT_STRENGTH_COLUMN = "ewpt_ew_true_over_T"

EQ418_PLOT_OBSERVABLES = (
    ("eq418_c_phiphi", r"$c_{\Phi\Phi}$"),
    ("eq418_c_xx", r"$c_{XX}$"),
    ("eq418_c_ss", r"$c_{SS}$"),
    ("eq418_C_phix", r"$C_{\Phi X}$"),
    ("eq418_C_phis", r"$C_{\Phi S}$"),
    ("eq418_C_xs", r"$C_{XS}$"),
    ("eq418_D", r"$D$"),
)

# Add future plots here.  In most cases you should only need to copy one of
# these PlotSpec entries and change x, y, labels, and optional color/size/marker
# observables.
@dataclass(frozen=True)
class PlotSpec:
    name: str
    x: str
    y: str
    xlabel: str | None = None
    ylabel: str | None = None
    title: str | None = None
    output_stem: str | None = None
    color_by: str | None = None
    size_by: str | None = None
    marker_by: str | None = None
    color: str = "tab:blue"
    marker: str = "o"
    size: float = 36.0
    alpha: float = 0.75
    cmap: str = "viridis"
    x_log: bool = False
    y_log: bool = False
    x_symlog: bool = False
    y_symlog: bool = False
    symlog_linthresh: float | None = None
    size_min: float = 24.0
    size_max: float = 120.0


PLOT_PRESETS = {
    "ewpt_vs_M2": PlotSpec(
        name="ewpt_vs_M2",
        x="M2",
        y=EWPT_STRENGTH_COLUMN,
        xlabel="M2 [GeV]",
        ylabel="EWPT strength ew_true/T",
        title="EWPT strength vs M2",
        output_stem="ewpt_ew_true_over_T_vs_M2",
    ),
    "ewpt_vs_lPhi": PlotSpec(
        name="ewpt_vs_lPhi",
        x="lPhi",
        y=EWPT_STRENGTH_COLUMN,
        xlabel=r"$\lambda_\Phi$",
        ylabel="EWPT strength ew_true/T",
        title="EWPT strength vs lPhi",
        output_stem="ewpt_ew_true_over_T_vs_lPhi",
    ),
    # Example to adapt:
    # "ewpt_vs_M2_colored_by_omega": PlotSpec(
    #     name="ewpt_vs_M2_colored_by_omega",
    #     x="M2",
    #     y="ewpt_ew_true_over_T",
    #     color_by="dm_omega",
    #     xlabel="M2 [GeV]",
    #     ylabel="EWPT strength ew_true/T",
    #     output_stem="ewpt_ew_true_over_T_vs_M2_colored_by_dm_omega",
    # ),
}

for observable, xlabel in EQ418_PLOT_OBSERVABLES:
    PLOT_PRESETS[f"ewpt_vs_{observable}"] = PlotSpec(
        name=f"ewpt_vs_{observable}",
        x=observable,
        y=EWPT_STRENGTH_COLUMN,
        xlabel=xlabel,
        ylabel="EWPT strength ew_true/T",
        title=f"EWPT strength vs {observable}",
        output_stem=f"{EWPT_STRENGTH_COLUMN}_vs_{observable}",
    )
    PLOT_PRESETS[f"ewpt_vs_{observable}_log"] = PlotSpec(
        name=f"ewpt_vs_{observable}_log",
        x=observable,
        y=EWPT_STRENGTH_COLUMN,
        xlabel=xlabel,
        ylabel="EWPT strength ew_true/T",
        title=f"EWPT strength vs {observable} (symlog x)",
        output_stem=f"{EWPT_STRENGTH_COLUMN}_vs_{observable}_symlogx",
        x_symlog=True,
    )

PLOT_GROUPS = {
    "eq418_all": [f"ewpt_vs_{observable}" for observable, _label in EQ418_PLOT_OBSERVABLES],
    "eq418_all_log": [
        f"ewpt_vs_{observable}_log"
        for observable, _label in EQ418_PLOT_OBSERVABLES
    ],
}

MARKER_CYCLE = ["o", "s", "^", "D", "P", "X", "v", "<", ">", "*"]


def finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isfinite(number):
        return number
    return None


def obs(row, column):
    value = finite_float(row.get(column))
    if value is None:
        raise ValueError(f"{column} is not finite")
    return value


def obs_any(row, *columns):
    for column in columns:
        value = finite_float(row.get(column))
        if value is not None:
            return value
    raise ValueError(f"none of {', '.join(columns)} is finite")


def obs_or_default(row, column, default):
    value = finite_float(row.get(column))
    if value is None:
        return default
    return value


def obs_any_or_default(row, columns, default):
    for column in columns:
        value = finite_float(row.get(column))
        if value is not None:
            return value
    return default


def safe_divide(numerator, denominator):
    if denominator == 0.0:
        raise ValueError("division by zero")
    return numerator / denominator


def safe_log10(value):
    if value <= 0.0:
        raise ValueError("log10 requires a positive value")
    return math.log10(value)


def lambda_phi(row):
    m1 = obs_any_or_default(row, ("M1", "m1"), DEFAULT_M1)
    m2 = obs_any(row, "M2", "m2")
    a12 = obs(row, "a12")
    return (
        m1**2 * math.cos(a12) ** 2
        + m2**2 * math.sin(a12) ** 2
    ) / (2.0 * SM_VEV**2)


def lambda_ss(row):
    m1 = obs_any_or_default(row, ("M1", "m1"), DEFAULT_M1)
    m2 = obs_any(row, "M2", "m2")
    vs = obs(row, "vs")
    a12 = obs(row, "a12")
    return (
        m2**2 * math.cos(a12) ** 2
        + m1**2 * math.sin(a12) ** 2
    ) / (2.0 * vs**2)


def lambda_phis(row):
    m1 = obs_any_or_default(row, ("M1", "m1"), DEFAULT_M1)
    m2 = obs_any(row, "M2", "m2")
    vs = obs(row, "vs")
    a12 = obs(row, "a12")
    return (
        (-m1**2 + m2**2) * math.cos(a12) * math.sin(a12)
    ) / (SM_VEV * vs)


def trsm_tree_lambdas(row):
    return {
        "lambda_phiphi": lambda_phi(row),
        "lambda_xx": obs_any(row, "lX", "lx"),
        "lambda_ss": lambda_ss(row),
        "lambda_phix": obs_any(row, "lPhiX", "lphix"),
        "lambda_phis": lambda_phis(row),
        "lambda_xs": obs_any(row, "lSX", "lsx"),
    }


def eq418_couplings(row):
    lambdas = trsm_tree_lambdas(row)
    c_phi = 0.5 * lambdas["lambda_phiphi"]
    c_x = 0.5 * lambdas["lambda_xx"]
    c_s = 0.5 * lambdas["lambda_ss"]
    c_phix = 0.5 * lambdas["lambda_phix"]
    c_phis = 0.5 * lambdas["lambda_phis"]
    c_xs = 0.5 * lambdas["lambda_xs"]
    c_phix_minor = c_phi * c_x - c_phix**2
    c_phis_minor = c_phi * c_s - c_phis**2
    c_xs_minor = c_x * c_s - c_xs**2
    determinant = (
        c_phi * c_x * c_s
        + 2.0 * c_phix * c_phis * c_xs
        - c_phis**2 * c_x
        - c_phi * c_xs**2
        - c_phix**2 * c_s
    )
    return {
        "eq418_c_phiphi": c_phi,
        "eq418_c_xx": c_x,
        "eq418_c_ss": c_s,
        "eq418_c_phix": c_phix,
        "eq418_c_phis": c_phis,
        "eq418_c_xs": c_xs,
        "eq418_C_phix": c_phix_minor,
        "eq418_C_phis": c_phis_minor,
        "eq418_C_xs": c_xs_minor,
        "eq418_D": determinant,
    }


def eq418_value(row, column):
    return eq418_couplings(row)[column]


# Derived observables are ordinary Python functions of one row.  Add your own
# here and then use the key name in --x, --y, --color-by, or --size-by.
DERIVED_OBSERVABLES = {
    "M2_over_M3": lambda row: safe_divide(obs(row, "M2"), obs(row, "M3")),
    "M3_minus_M2": lambda row: obs(row, "M3") - obs(row, "M2"),
    "abs_a12": lambda row: abs(obs(row, "a12")),
    "log10_dm_omega": lambda row: safe_log10(obs(row, "dm_omega")),
    "log10_dm_dir_det": lambda row: safe_log10(obs(row, "dm_dir_det")),
    "lPhi": lambda_phi,
    "lS": lambda_ss,
    "lPhiS": lambda_phis,
    "eq418_c_phiphi": lambda row: eq418_value(row, "eq418_c_phiphi"),
    "eq418_c_xx": lambda row: eq418_value(row, "eq418_c_xx"),
    "eq418_c_ss": lambda row: eq418_value(row, "eq418_c_ss"),
    "eq418_c_phix": lambda row: eq418_value(row, "eq418_c_phix"),
    "eq418_c_phis": lambda row: eq418_value(row, "eq418_c_phis"),
    "eq418_c_xs": lambda row: eq418_value(row, "eq418_c_xs"),
    "eq418_C_phix": lambda row: eq418_value(row, "eq418_C_phix"),
    "eq418_C_phis": lambda row: eq418_value(row, "eq418_C_phis"),
    "eq418_C_xs": lambda row: eq418_value(row, "eq418_C_xs"),
    "eq418_D": lambda row: eq418_value(row, "eq418_D"),
}

def format_number(value):
    if value is None:
        return "nan"
    return f"{value:.12g}"


def add_derived_observables(rows, derived_observables=None):
    if derived_observables is None:
        derived_observables = DERIVED_OBSERVABLES
    for row in rows:
        for name, function in derived_observables.items():
            try:
                row[name] = format_number(function(row))
            except (ArithmeticError, ValueError, TypeError):
                row[name] = "nan"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Make scatter plots from TRSM scan TSV files."
    )
    parser.add_argument("input", type=Path, help="TRSM scan TSV file.")
    parser.add_argument(
        "--plot",
        default="ewpt_vs_M2",
        choices=sorted(PLOT_PRESETS) + sorted(PLOT_GROUPS),
        help="Named preset from PLOT_PRESETS, or a plot group such as eq418_all.",
    )
    parser.add_argument("--x", help="Override x-axis observable.")
    parser.add_argument("--y", help="Override y-axis observable.")
    parser.add_argument("--color-by", help="Observable used for marker colors.")
    parser.add_argument("--size-by", help="Observable used for marker sizes.")
    parser.add_argument("--marker-by", help="Observable used for marker shapes.")
    parser.add_argument("--xlabel", help="Override x-axis label.")
    parser.add_argument("--ylabel", help="Override y-axis label.")
    parser.add_argument("--title", help="Override plot title.")
    parser.add_argument("--output-stem", help="Output basename without extension.")
    parser.add_argument("--output-dir", type=Path, default=Path("plots"))
    parser.add_argument("--format", choices=["png", "pdf", "both"], default="both")
    parser.add_argument("--marker", default=None, help="Matplotlib marker for one-shape plots.")
    parser.add_argument("--color", default=None, help="Matplotlib color for one-color plots.")
    parser.add_argument("--size", type=float, default=None, help="Marker size for one-size plots.")
    parser.add_argument("--alpha", type=float, default=None)
    parser.add_argument("--cmap", default=None, help="Matplotlib colormap for numeric --color-by.")
    parser.add_argument("--x-log", action="store_true")
    parser.add_argument("--y-log", action="store_true")
    parser.add_argument(
        "--x-symlog",
        action="store_true",
        help="Use a symmetric-log x-axis, preserving negative values.",
    )
    parser.add_argument(
        "--y-symlog",
        action="store_true",
        help="Use a symmetric-log y-axis, preserving negative values.",
    )
    parser.add_argument(
        "--symlog-linthresh",
        type=float,
        help="Linear threshold for symlog axes. Defaults to one tenth of the smallest nonzero absolute axis value.",
    )
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="Print available raw and derived columns and exit without plotting.",
    )
    parser.add_argument(
        "--list-plots",
        action="store_true",
        help="Print available preset names and exit without plotting.",
    )
    return parser.parse_args(argv)


def load_rows(path):
    path = Path(path)
    with path.open(encoding="ascii", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    add_derived_observables(rows)
    return rows


def rows_with_finite_columns(rows, columns):
    required = [column for column in columns if column is not None]
    return [
        row
        for row in rows
        if all(finite_float(row.get(column)) is not None for column in required)
    ]


def rows_for_plot(rows, spec):
    required_columns = [spec.x, spec.y, spec.color_by, spec.size_by]
    rows = rows_with_finite_columns(rows, required_columns)
    if spec.x_log and not spec.x_symlog:
        rows = [row for row in rows if finite_float(row.get(spec.x)) > 0.0]
    if spec.y_log and not spec.y_symlog:
        rows = [row for row in rows if finite_float(row.get(spec.y)) > 0.0]
    return rows


def resolve_one_plot_spec(plot_name, args):
    spec = PLOT_PRESETS[plot_name]
    overrides = {}
    for attr in [
        "x",
        "y",
        "color_by",
        "size_by",
        "marker_by",
        "xlabel",
        "ylabel",
        "title",
        "output_stem",
        "marker",
        "color",
        "size",
        "alpha",
        "cmap",
        "symlog_linthresh",
    ]:
        value = getattr(args, attr)
        if value is not None:
            overrides[attr] = value
    if args.x_log:
        overrides["x_log"] = True
        overrides["x_symlog"] = False
    if args.y_log:
        overrides["y_log"] = True
        overrides["y_symlog"] = False
    if args.x_symlog:
        overrides["x_symlog"] = True
        overrides["x_log"] = False
    if args.y_symlog:
        overrides["y_symlog"] = True
        overrides["y_log"] = False
    return replace(spec, **overrides)


def resolve_plot_spec(args):
    return resolve_one_plot_spec(args.plot, args)


def resolve_plot_specs(args):
    plot_names = PLOT_GROUPS.get(args.plot, [args.plot])
    return [resolve_one_plot_spec(plot_name, args) for plot_name in plot_names]


def output_paths(output_dir, output_stem, plot_format):
    output_dir = Path(output_dir)
    if plot_format == "both":
        extensions = ["png", "pdf"]
    else:
        extensions = [plot_format]
    return [output_dir / f"{output_stem}.{extension}" for extension in extensions]


def numeric_values(rows, column):
    return [finite_float(row[column]) for row in rows]


def infer_symlog_linthresh(values):
    absolute_nonzero_values = [abs(value) for value in values if value != 0.0]
    if not absolute_nonzero_values:
        return 1.0
    return min(absolute_nonzero_values) / 10.0


def scaled_sizes(rows, spec):
    if spec.size_by is None:
        return [spec.size for _row in rows]
    values = numeric_values(rows, spec.size_by)
    low = min(values)
    high = max(values)
    if math.isclose(low, high):
        return [(spec.size_min + spec.size_max) / 2.0 for _value in values]
    return [
        spec.size_min
        + (value - low) / (high - low) * (spec.size_max - spec.size_min)
        for value in values
    ]


def marker_groups(rows, spec):
    if spec.marker_by is None:
        return [(spec.marker, None, rows)]
    labels = sorted({row.get(spec.marker_by, "missing") for row in rows})
    marker_for_label = {
        label: MARKER_CYCLE[index % len(MARKER_CYCLE)]
        for index, label in enumerate(labels)
    }
    return [
        (
            marker_for_label[label],
            label,
            [row for row in rows if row.get(spec.marker_by, "missing") == label],
        )
        for label in labels
    ]


def color_values(rows, spec):
    if spec.color_by is None:
        return spec.color, None
    values = [finite_float(row.get(spec.color_by)) for row in rows]
    if all(value is not None for value in values):
        return values, "numeric"
    labels = sorted({row.get(spec.color_by, "missing") for row in rows})
    color_for_label = {label: f"C{index % 10}" for index, label in enumerate(labels)}
    return [color_for_label[row.get(spec.color_by, "missing")] for row in rows], "categorical"


def plot_scatter(rows, spec, output_dir, plot_format):
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required to make plots") from exc

    rows = rows_for_plot(rows, spec)
    if not rows:
        raise ValueError("No rows have finite values for the requested observables.")

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    last_numeric_scatter = None
    all_sizes = scaled_sizes(rows, spec)
    row_to_size = {id(row): size for row, size in zip(rows, all_sizes)}

    for marker, marker_label, group_rows in marker_groups(rows, spec):
        x_values = numeric_values(group_rows, spec.x)
        y_values = numeric_values(group_rows, spec.y)
        sizes = [row_to_size[id(row)] for row in group_rows]
        colors, color_kind = color_values(group_rows, spec)
        scatter = ax.scatter(
            x_values,
            y_values,
            s=sizes,
            c=colors,
            marker=marker,
            alpha=spec.alpha,
            cmap=spec.cmap if color_kind == "numeric" else None,
            edgecolors="black",
            linewidths=0.25,
            label=marker_label,
        )
        if color_kind == "numeric":
            last_numeric_scatter = scatter

    ax.set_xlabel(spec.xlabel or spec.x)
    ax.set_ylabel(spec.ylabel or spec.y)
    if spec.title:
        ax.set_title(spec.title)
    if spec.x_symlog:
        ax.set_xscale(
            "symlog",
            linthresh=spec.symlog_linthresh
            or infer_symlog_linthresh(numeric_values(rows, spec.x)),
        )
    elif spec.x_log:
        ax.set_xscale("log")
    if spec.y_symlog:
        ax.set_yscale(
            "symlog",
            linthresh=spec.symlog_linthresh
            or infer_symlog_linthresh(numeric_values(rows, spec.y)),
        )
    elif spec.y_log:
        ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    if spec.marker_by is not None:
        ax.legend(title=spec.marker_by, frameon=False)
    if last_numeric_scatter is not None:
        fig.colorbar(last_numeric_scatter, ax=ax, label=spec.color_by)
    fig.tight_layout()

    paths = output_paths(output_dir, spec.output_stem or f"{spec.y}_vs_{spec.x}", plot_format)
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return paths


def print_columns(rows):
    if not rows:
        return
    for column in rows[0].keys():
        print(column)


def print_presets():
    for name, spec in sorted(PLOT_PRESETS.items()):
        print(f"{name}: {spec.y} vs {spec.x}")
    for name, plot_names in sorted(PLOT_GROUPS.items()):
        print(f"{name}: {len(plot_names)} separate plots")


def run(argv=None):
    args = parse_args(argv)
    if args.list_plots:
        print_presets()
        return []
    rows = load_rows(args.input)
    if args.list_columns:
        print_columns(rows)
        return []
    paths = []
    for spec in resolve_plot_specs(args):
        spec_paths = plot_scatter(rows, spec, args.output_dir, args.format)
        paths.extend(spec_paths)
        for path in spec_paths:
            print(f"Saved {path}")
    return paths


if __name__ == "__main__":
    run()
