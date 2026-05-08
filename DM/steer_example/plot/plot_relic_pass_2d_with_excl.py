#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

PLOT_COLUMNS = [
    "index",
    "LX",
    "LHX",
    "LSX",
    "MX",
    "vevs",
    "SinT",
    "Mh2",
    "MDM",
    "Omega",
    "DirDet",
]

DM_COLUMNS = [
    "index",
    "LX",
    "LHX",
    "LSX",
    "MX",
    "vevs",
    "SinT",
    "Mh2",
    "MDM",
    "Omega",
    "DirDet",
    "DirDetLimit",
    "LUXBaseLimit",
]

RELIC_CENTRAL = 0.120
RELIC_TOLERANCE = 0.001
RELIC_MIN = RELIC_CENTRAL - RELIC_TOLERANCE
RELIC_MAX = RELIC_CENTRAL + RELIC_TOLERANCE

ACCEPTED_CMAP = plt.cm.rainbow
EXCLUDED_CMAP = plt.cm.Greys


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot accepted DM points and overlay experimentally excluded points "
            f"in greyscale for two chosen variables, keeping only Omega <= {RELIC_MAX}."
        )
    )
    parser.add_argument(
        "output_dir",
        help="Directory containing allall.dat and dmexcl.dat, typically dmonly/output.",
    )
    parser.add_argument(
        "x_variable",
        help="Variable for the x-axis. Choices: " + ", ".join(PLOT_COLUMNS),
    )
    parser.add_argument(
        "y_variable",
        help="Variable for the y-axis. Choices: " + ", ".join(PLOT_COLUMNS),
    )
    parser.add_argument(
        "--output",
        help=(
            "Optional output image path. Defaults to "
            "<plot_dir>/relic_pass_with_excl_<y_variable>_vs_<x_variable>.png"
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Also display the plot interactively.",
    )
    return parser.parse_args()


def load_rows(data_file: Path, columns):
    rows = []
    with data_file.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) != len(columns):
                raise ValueError(
                    f"{data_file}:{line_number}: expected {len(columns)} columns, "
                    f"found {len(parts)}"
                )
            try:
                rows.append(
                    {name: float(value) for name, value in zip(columns, parts)}
                )
            except ValueError as exc:
                raise ValueError(
                    f"{data_file}:{line_number}: could not parse numeric values "
                    f"from line: {stripped}"
                ) from exc
    return rows


def validate_variable(name: str):
    if name not in PLOT_COLUMNS:
        raise SystemExit(
            f"Unknown variable '{name}'. Valid choices: {', '.join(PLOT_COLUMNS)}"
        )


def build_normalization(rows, label):
    omega_values = [row["Omega"] for row in rows]
    if not omega_values:
        raise SystemExit(f"No points found for {label}")
    omega_min = min(omega_values)
    omega_max = max(omega_values)
    if omega_min == omega_max:
        omega_max = omega_min + 1e-12
    return colors.Normalize(vmin=omega_min, vmax=omega_max)


def split_coordinates(rows, x_variable, y_variable, log_y=False):
    x_values = [row[x_variable] for row in rows]
    y_source = [row[y_variable] for row in rows]
    if log_y:
        y_values = np.log10(y_source)
    else:
        y_values = y_source
    omega_values = [row["Omega"] for row in rows]
    return x_values, y_values, omega_values


def add_scatter_layers(
    axis,
    accepted_rows,
    excluded_rows,
    x_variable,
    y_variable,
    accepted_norm,
    excluded_norm,
    log_y=False,
):
    excluded_scatter = None
    accepted_scatter = None

    if excluded_rows:
        excl_x, excl_y, excl_omega = split_coordinates(
            excluded_rows, x_variable, y_variable, log_y=log_y
        )
        excluded_scatter = axis.scatter(
            excl_x,
            excl_y,
            c=excl_omega,
            cmap=EXCLUDED_CMAP,
            norm=excluded_norm,
            s=24,
            edgecolors="none",
            alpha=0.9,
            # label="Excluded points",
        )

    if accepted_rows:
        acc_x, acc_y, acc_omega = split_coordinates(
            accepted_rows, x_variable, y_variable, log_y=log_y
        )
        accepted_scatter = axis.scatter(
            acc_x,
            acc_y,
            c=acc_omega,
            cmap=ACCEPTED_CMAP,
            norm=accepted_norm,
            s=28,
            edgecolors="none",
            # label="Accepted points",
        )

        central_rows = [
            row for row in accepted_rows if RELIC_MIN <= row["Omega"] <= RELIC_MAX
        ]
        if central_rows:
            central_x, central_y, _ = split_coordinates(
                central_rows, x_variable, y_variable, log_y=log_y
            )
            axis.scatter(
                central_x,
                central_y,
                color="magenta",
                s=40,
                linewidths=0.5,
                label="Central relic range",
            )

    return accepted_scatter, excluded_scatter


def save_plot(
    accepted_rows,
    excluded_rows,
    x_variable,
    y_variable,
    accepted_norm,
    excluded_norm,
    output_path,
    log_y=False,
    show=False,
):
    figure, axis = plt.subplots(figsize=(7, 5))
    accepted_scatter, excluded_scatter = add_scatter_layers(
        axis,
        accepted_rows,
        excluded_rows,
        x_variable,
        y_variable,
        accepted_norm,
        excluded_norm,
        log_y=log_y,
    )

    axis.set_xlabel(x_variable)
    axis.set_ylabel(f"log10({y_variable})" if log_y else y_variable)
    axis.set_title(
        "Accepted points coloured by Omega, excluded points in greyscale"
    )
    axis.grid(True, alpha=0.3)

    if accepted_scatter is not None:
        figure.colorbar(accepted_scatter, ax=axis, label="Omega (accepted)")
    if excluded_scatter is not None:
        figure.colorbar(excluded_scatter, ax=axis, label="Omega (excluded)")

    handles, labels = axis.get_legend_handles_labels()
    if handles:
        axis.legend(loc="best")

    figure.tight_layout()
    figure.savefig(output_path, dpi=200)

    if show:
        plt.show()
    else:
        plt.close(figure)


def main():
    args = parse_args()
    validate_variable(args.x_variable)
    validate_variable(args.y_variable)

    output_dir = Path(args.output_dir)
    accepted_file = output_dir / "allall.dat"
    excluded_file = output_dir / "dmexcl.dat"

    if not accepted_file.is_file():
        raise SystemExit(f"Could not find {accepted_file}")

    script_dir = Path(__file__).resolve().parent
    outplots_dir = script_dir / "outplots"
    outplots_dir.mkdir(parents=True, exist_ok=True)

    accepted_rows = load_rows(accepted_file, DM_COLUMNS)
    excluded_rows = load_rows(excluded_file, DM_COLUMNS) if excluded_file.is_file() else []
    accepted_rows = [row for row in accepted_rows if row["Omega"] <= RELIC_MAX]
    excluded_rows = [row for row in excluded_rows if row["Omega"] <= RELIC_MAX]
    if not accepted_rows and not excluded_rows:
        raise SystemExit(
            f"No points in {output_dir} satisfy Omega <= {RELIC_MAX}"
        )
    accepted_norm = build_normalization(accepted_rows, accepted_file.name)
    excluded_norm = (
        build_normalization(excluded_rows, excluded_file.name)
        if excluded_rows
        else None
    )

    print(
        f"Loaded {len(accepted_rows)} accepted points from {accepted_file} "
        f"with Omega <= {RELIC_MAX}"
    )
    print(
        f"Loaded {len(excluded_rows)} excluded points from {excluded_file} "
        f"with Omega <= {RELIC_MAX}"
    )
    strict_count = sum(
        1 for row in accepted_rows if RELIC_MIN <= row["Omega"] <= RELIC_MAX
    )
    print(
        f"Found {strict_count} accepted points with "
        f"{RELIC_MIN} <= Omega <= {RELIC_MAX}"
    )

    default_output = (
        outplots_dir / f"relic_pass_with_excl_{args.y_variable}_vs_{args.x_variable}.png"
    )
    linear_output = Path(args.output) if args.output else default_output
    save_plot(
        accepted_rows,
        excluded_rows,
        args.x_variable,
        args.y_variable,
        accepted_norm,
        excluded_norm,
        linear_output,
        show=args.show,
    )

    log_output = (
        outplots_dir
        / f"relic_pass_with_excl_log10_{args.y_variable}_vs_{args.x_variable}.png"
    )
    save_plot(
        accepted_rows,
        excluded_rows,
        args.x_variable,
        args.y_variable,
        accepted_norm,
        excluded_norm,
        log_output,
        log_y=True,
        show=args.show,
    )

    print(f"Saved plots to {linear_output} and {log_output}")


if __name__ == "__main__":
    main()
