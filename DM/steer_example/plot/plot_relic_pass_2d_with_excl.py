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

INDIRECT_COLUMNS = [
    *DM_COLUMNS,
    "IndirAvailable",
    "IndirEnergy",
    "IndirFlux",
    "IndirLimit",
    "IndirRatio",
]

RELIC_CENTRAL = 0.120
RELIC_TOLERANCE = 0.001
RELIC_MIN = RELIC_CENTRAL - RELIC_TOLERANCE
RELIC_MAX = RELIC_CENTRAL + RELIC_TOLERANCE

ACCEPTED_CMAP = plt.cm.rainbow
DIR_EXCLUDED_CMAP = plt.cm.Greys
INDIR_EXCLUDED_CMAP = plt.cm.copper
# INDIRECT_FAIL_COLOR = "saddlebrown"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot relic- and direct-detection-passing DM points and overlay experimentally excluded points "
            f"in greyscale for two chosen variables, keeping only Omega <= {RELIC_MAX}."
        )
    )
    parser.add_argument(
        "output_dir",
        help="Directory containing all_dirpass.dat and dmexcl.dat, typically dmonly/output.",
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
    parser.add_argument(
        "--title-suffix",
        default="",
        help="Optional text appended to generated plot titles.",
    )
    parser.add_argument(
        "--mark-indirect-fail",
        action="store_true",
        help=(
            "Overlay relic- and direct-detection-passing points that fail indirect detection "
            "from indirexcl.dat in brown."
        ),
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


def row_key(row):
    return int(row["index"])


def add_scatter_layers(
    axis,
    accepted_rows,
    dir_excluded_rows,
    indir_excluded_rows,
    x_variable,
    y_variable,
    accepted_norm,
    dir_excluded_norm,
    indir_excluded_norm,
    log_y=False,
):
    dir_scatter = None
    accepted_scatter = None
    indir_scatter = None

    if dir_excluded_rows:
        excl_x, excl_y, excl_omega = split_coordinates(
            dir_excluded_rows, x_variable, y_variable, log_y=log_y
        )
        dir_scatter = axis.scatter(
            excl_x,
            excl_y,
            c=excl_omega,
            cmap=DIR_EXCLUDED_CMAP,
            norm=dir_excluded_norm,
            s=24,
            edgecolors="none",
            alpha=0.9,
            # label="Excluded points",
        )
    
    if indir_excluded_rows:
        excl_x, excl_y, excl_omega = split_coordinates(
            indir_excluded_rows, x_variable, y_variable, log_y=log_y
        )
        indir_scatter = axis.scatter(
            excl_x,
            excl_y,
            c=excl_omega,
            cmap=INDIR_EXCLUDED_CMAP,
            norm=indir_excluded_norm,
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

    # if indir_excluded_rows:
    #     indirect_x, indirect_y, _ = split_coordinates(
    #         indir_excluded_rows, x_variable, y_variable, log_y=log_y
    #     )
    #     indir_scatter = axis.scatter(
    #         indirect_x,
    #         indirect_y,
    #         color=INDIRECT_FAIL_COLOR,
    #         s=42,
    #         marker="s",
    #         edgecolors="black",
    #         linewidths=0.35,
    #         alpha=0.95,
    #         label="Fails indirect detection",
    #     )

    return accepted_scatter, dir_scatter, indir_scatter


def save_plot(
    accepted_rows,
    dir_excluded_rows,
    indir_excluded_rows,
    x_variable,
    y_variable,
    accepted_norm,
    dir_excluded_norm,
    indir_excluded_norm,
    output_path,
    title_suffix="",
    log_y=False,
    show=False,
):
    figure, axis = plt.subplots(figsize=(7, 5))
    accepted_scatter, dir_scatter, indir_scatter = add_scatter_layers(
        axis,
        accepted_rows,
        dir_excluded_rows,
        indir_excluded_rows,
        x_variable,
        y_variable,
        accepted_norm,
        dir_excluded_norm,
        indir_excluded_norm,
        log_y=log_y,
    )

    axis.set_xlabel(x_variable)
    axis.set_ylabel(f"log10({y_variable})" if log_y else y_variable)
    title = "Relic-passing points coloured by Omega"
    # if indir_excluded_rows:
    #     title += ", indirect failures in brown"
    axis.set_title(f"{title} {title_suffix}".strip())
    axis.grid(True, alpha=0.3)

    if accepted_scatter is not None:
        figure.colorbar(accepted_scatter, ax=axis, label="Omega (accepted)")
    if dir_scatter is not None:
        figure.colorbar(dir_scatter, ax=axis, label="Omega (excluded by direct detection)")
    if indir_scatter is not None:
        figure.colorbar(indir_scatter, ax=axis, label="Omega (excluded by indirect detection)")

    handles, _ = axis.get_legend_handles_labels()
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
    accepted_file = output_dir / "all_dirpass.dat"
    dir_excluded_file = output_dir / "dmexcl.dat"
    indirect_file = output_dir / "indirexcl.dat"

    if not (accepted_file.is_file() and dir_excluded_file.is_file() and indirect_file.is_file()):
        raise SystemExit(f"Could not find one or more required files")

    script_dir = Path(__file__).resolve().parent
    outplots_dir = script_dir / "outplots"
    outplots_dir.mkdir(parents=True, exist_ok=True)

    accepted_rows = load_rows(accepted_file, DM_COLUMNS)
    dir_excluded_rows = load_rows(dir_excluded_file, DM_COLUMNS) if dir_excluded_file.is_file() else []
    indirect_rows = (
        load_rows(indirect_file, INDIRECT_COLUMNS)
        if args.mark_indirect_fail and indirect_file.is_file()
        else []
    )
    accepted_rows = [row for row in accepted_rows if row["Omega"] <= RELIC_MAX]
    dir_excluded_rows = [row for row in dir_excluded_rows if row["Omega"] <= RELIC_MAX]
    indirect_rows = [row for row in indirect_rows if row["Omega"] <= RELIC_MAX]

    accepted_keys = {row_key(row) for row in accepted_rows}
    dir_excluded_rows = [row for row in dir_excluded_rows if row_key(row) not in accepted_keys]
    indir_excluded_rows = [
        row for row in indirect_rows if row_key(row) in accepted_keys
    ]

    if not accepted_rows and not dir_excluded_rows:
        raise SystemExit(
            f"No points in {output_dir} satisfy Omega <= {RELIC_MAX}"
        )
    accepted_norm = build_normalization(accepted_rows, accepted_file.name)
    dir_excluded_norm = (
        build_normalization(dir_excluded_rows, dir_excluded_file.name)
        if dir_excluded_rows
        else None
    )
    indir_excluded_norm = (
        build_normalization(indir_excluded_rows, indirect_file.name)
        if indir_excluded_rows
        else None
    )

    print(
        f"Loaded {len(accepted_rows)} relic- and direct-detection-passing points from {accepted_file} "
        f"with Omega <= {RELIC_MAX}"
    )
    print(
        f"Loaded {len(dir_excluded_rows)} excluded points from {dir_excluded_file} "
        f"with Omega <= {RELIC_MAX}"
    )
    if args.mark_indirect_fail:
        print(
            f"Loaded {len(indir_excluded_rows)} relic- and direct-detection-passing indirect-failed "
            f"points from {indirect_file}"
        )
    strict_count = sum(
        1 for row in accepted_rows if RELIC_MIN <= row["Omega"] <= RELIC_MAX
    )
    print(
        f"Found {strict_count} accepted points with "
        f"{RELIC_MIN} <= Omega <= {RELIC_MAX}"
    )

    output_tag = (
        "relic_pass_with_excl_indirect"
        if args.mark_indirect_fail
        else "relic_pass_with_excl"
    )
    default_output = outplots_dir / f"{output_tag}_{args.y_variable}_vs_{args.x_variable}.png"
    linear_output = Path(args.output) if args.output else default_output
    save_plot(
        accepted_rows,
        dir_excluded_rows,
        indir_excluded_rows,
        args.x_variable,
        args.y_variable,
        accepted_norm,
        dir_excluded_norm,
        indir_excluded_norm,
        linear_output,
        title_suffix=args.title_suffix,
        show=args.show,
    )

    log_output = (
        outplots_dir
        / f"{output_tag}_log10_{args.y_variable}_vs_{args.x_variable}.png"
    )
    save_plot(
        accepted_rows,
        dir_excluded_rows,
        indir_excluded_rows,
        args.x_variable,
        args.y_variable,
        accepted_norm,
        dir_excluded_norm,
        indir_excluded_norm,
        log_output,
        title_suffix=args.title_suffix,
        log_y=True,
        show=args.show,
    )

    print(f"Saved plots to {linear_output} and {log_output}")


if __name__ == "__main__":
    main()
