#!/usr/bin/env python3

import argparse
import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path


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
    size_min: float = 24.0
    size_max: float = 120.0


PLOT_PRESETS = {
    "ewpt_vs_M2": PlotSpec(
        name="ewpt_vs_M2",
        x="M2",
        y="ewpt_ew_true_over_T",
        xlabel="M2 [GeV]",
        ylabel="EWPT strength ew_true/T",
        title="EWPT strength vs M2",
        output_stem="ewpt_ew_true_over_T_vs_M2",
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


def safe_divide(numerator, denominator):
    if denominator == 0.0:
        raise ValueError("division by zero")
    return numerator / denominator


def safe_log10(value):
    if value <= 0.0:
        raise ValueError("log10 requires a positive value")
    return math.log10(value)


# Derived observables are ordinary Python functions of one row.  Add your own
# here and then use the key name in --x, --y, --color-by, or --size-by.
DERIVED_OBSERVABLES = {
    "M2_over_M3": lambda row: safe_divide(obs(row, "M2"), obs(row, "M3")),
    "M3_minus_M2": lambda row: obs(row, "M3") - obs(row, "M2"),
    "abs_a12": lambda row: abs(obs(row, "a12")),
    "log10_dm_omega": lambda row: safe_log10(obs(row, "dm_omega")),
    "log10_dm_dir_det": lambda row: safe_log10(obs(row, "dm_dir_det")),
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
        choices=sorted(PLOT_PRESETS),
        help="Named preset from PLOT_PRESETS.",
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


def resolve_plot_spec(args):
    spec = PLOT_PRESETS[args.plot]
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
    ]:
        value = getattr(args, attr)
        if value is not None:
            overrides[attr] = value
    if args.x_log:
        overrides["x_log"] = True
    if args.y_log:
        overrides["y_log"] = True
    return replace(spec, **overrides)


def output_paths(output_dir, output_stem, plot_format):
    output_dir = Path(output_dir)
    if plot_format == "both":
        extensions = ["png", "pdf"]
    else:
        extensions = [plot_format]
    return [output_dir / f"{output_stem}.{extension}" for extension in extensions]


def numeric_values(rows, column):
    return [finite_float(row[column]) for row in rows]


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

    required_columns = [spec.x, spec.y, spec.color_by, spec.size_by]
    rows = rows_with_finite_columns(rows, required_columns)
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
    if spec.x_log:
        ax.set_xscale("log")
    if spec.y_log:
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


def run(argv=None):
    args = parse_args(argv)
    if args.list_plots:
        print_presets()
        return []
    rows = load_rows(args.input)
    if args.list_columns:
        print_columns(rows)
        return []
    spec = resolve_plot_spec(args)
    paths = plot_scatter(rows, spec, args.output_dir, args.format)
    for path in paths:
        print(f"Saved {path}")
    return paths


if __name__ == "__main__":
    run()
