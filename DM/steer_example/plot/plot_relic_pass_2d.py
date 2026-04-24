#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

COLUMNS = [
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

RELIC_LIMIT = 0.1224
RELIC_CENTRAL = 0.1199
RELIC_TOLERANCE = 0.025
RELIC_MIN = RELIC_CENTRAL - RELIC_TOLERANCE
RELIC_MAX = RELIC_CENTRAL + RELIC_TOLERANCE

OMEGA_CMAP = LinearSegmentedColormap.from_list("omega", ["lightblue", "red"])


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot points that pass the relic-density cut "
            f"(Omega <= {RELIC_LIMIT}) for two chosen variables."
        )
    )
    parser.add_argument(
        "output_dir",
        help="Directory containing scan_results.dat, typically dmonly/output.",
    )
    parser.add_argument(
        "x_variable",
        help="Variable for the x-axis. Choices: " + ", ".join(COLUMNS),
    )
    parser.add_argument(
        "y_variable",
        help="Variable for the y-axis. Choices: " + ", ".join(COLUMNS),
    )
    parser.add_argument(
        "--output",
        help=(
            "Optional output image path. Defaults to "
            "<output_dir>/relic_pass_<x_variable>_vs_<y_variable>.png"
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Also display the plot interactively.",
    )
    return parser.parse_args()


def load_scan_results(scan_file: Path):
    rows = []
    with scan_file.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) != len(COLUMNS):
                raise ValueError(
                    f"{scan_file}:{line_number}: expected {len(COLUMNS)} columns, "
                    f"found {len(parts)}"
                )
            try:
                rows.append(
                    {name: float(value) for name, value in zip(COLUMNS, parts)}
                )
            except ValueError as exc:
                raise ValueError(
                    f"{scan_file}:{line_number}: could not parse numeric values "
                    f"from line: {stripped}"
                ) from exc
    if not rows:
        raise ValueError(f"{scan_file} is empty")
    return rows


def validate_variable(name: str):
    if name not in COLUMNS:
        raise SystemExit(
            f"Unknown variable '{name}'. Valid choices: {', '.join(COLUMNS)}"
        )


def main():
    args = parse_args()
    validate_variable(args.x_variable)
    validate_variable(args.y_variable)

    output_dir = Path(args.output_dir)
    scan_file = output_dir / "scan_results.dat"
    if not scan_file.is_file():
        raise SystemExit(f"Could not find {scan_file}")

    script_dir = Path(__file__).resolve().parent
    outplots_dir = script_dir / "outplots"
    outplots_dir.mkdir(parents=True, exist_ok=True)

    rows = load_scan_results(scan_file)
    passed_rows = [row for row in rows if row["Omega"] <= RELIC_LIMIT]
    if not passed_rows:
        raise SystemExit(
            f"No points in {scan_file} satisfy Omega <= {RELIC_LIMIT}"
        )
    
    x_values = [row[args.x_variable] for row in passed_rows]
    y_values = [row[args.y_variable] for row in passed_rows]
    omega_values = [row["Omega"] for row in passed_rows]
    central_rows = [
        row
        for row in rows
        if RELIC_MIN <= row["Omega"] <= RELIC_MAX
    ]
    central_x = [row[args.x_variable] for row in central_rows]
    central_y = [row[args.y_variable] for row in central_rows]
    print(f"Found {len(central_rows)} points in {scan_file} that satisfy {RELIC_MIN} <= Omega <= {RELIC_MAX}")

    figure, axis = plt.subplots(figsize=(7, 5))
    scatter = axis.scatter(
        x_values,
        y_values,
        c=omega_values,
        cmap=OMEGA_CMAP,
        s=28,
        edgecolors="none",
    )
    if central_rows:
        axis.scatter(
            central_x,
            central_y,
            color="purple",
            s=40,
            # edgecolors="purple",
            linewidths=0.5,
            label="Central relic range",
        )
    axis.set_xlabel(args.x_variable)
    axis.set_ylabel(args.y_variable)
    axis.set_title(
        f"Relic-density-passing points, purple points fit observed density"
    )
    axis.grid(True, alpha=0.3)
    figure.colorbar(scatter, ax=axis, label="Omega")
    figure.tight_layout()

    default_output = (
        outplots_dir
        / f"relic_pass_{args.y_variable}_vs_{args.x_variable}.png"
    )
    output_path = Path(args.output) if args.output else default_output
    figure.savefig(output_path, dpi=200)

    if args.show:
        plt.show()
    else:
        plt.close(figure)


    figure, axis = plt.subplots(figsize=(7, 5))
    scatter = axis.scatter(
        x_values,
        np.log10(y_values),
        c=omega_values,
        cmap=OMEGA_CMAP,
        s=28,
        edgecolors="none",
    )
    if central_rows:
        axis.scatter(
            central_x,
            np.log10(central_y),
            color="purple",
            s=40,
            # edgecolors="purple",
            linewidths=0.5,
            label="Central relic range",
        )
    axis.set_xlabel(args.x_variable)
    axis.set_ylabel(f"log10({args.y_variable})")
    axis.set_title(
        f"Relic-density-passing points, purple points fit observed density"
    )
    axis.grid(True, alpha=0.3)
    figure.colorbar(scatter, ax=axis, label="Omega")
    figure.tight_layout()

    default_output = (
        outplots_dir
        / f"relic_pass_log10_{args.y_variable}_vs_{args.x_variable}.png"
    )
    output_path = Path(args.output) if args.output else default_output
    figure.savefig(output_path, dpi=200)

    if args.show:
        plt.show()
    else:
        plt.close(figure)

    print(f"Saved plots to {output_path}")


if __name__ == "__main__":
    main()
