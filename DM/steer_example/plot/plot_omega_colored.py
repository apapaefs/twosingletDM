#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Omega against a variable, colored by another variable."
    )
    parser.add_argument(
        "output_dir",
        help="Directory containing scan_results.dat, typically dmonly/output.",
    )
    parser.add_argument(
        "x_variable",
        help="Variable to plot on the x-axis. Choices: "
        + ", ".join(name for name in COLUMNS if name != "Omega"),
    )
    parser.add_argument(
        "color_variable",
        help="Variable to use for color mapping. Choices: "
        + ", ".join(name for name in COLUMNS if name != "Omega"),
    )
    parser.add_argument(
        "--output",
        help="Optional output image path. Defaults to <output_dir>/omega_vs_<x_variable>_colored_by_<color_variable>.png",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Also display the plot interactively.",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap name. Defaults to 'viridis'.",
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


def main():
    args = parse_args()

    x_variable = args.x_variable
    color_variable = args.color_variable

    if x_variable not in COLUMNS:
        raise SystemExit(
            f"Unknown x_variable '{x_variable}'. Valid choices: {', '.join(COLUMNS)}"
        )
    if color_variable not in COLUMNS:
        raise SystemExit(
            f"Unknown color_variable '{color_variable}'. Valid choices: {', '.join(COLUMNS)}"
        )
    if x_variable == "Omega":
        raise SystemExit("Use a variable other than Omega for the x-axis.")
    if color_variable == "Omega":
        raise SystemExit("Use a variable other than Omega for the color mapping.")
    if x_variable == color_variable:
        raise SystemExit("x_variable and color_variable must be different.")

    output_dir = Path(args.output_dir)
    scan_file = output_dir / "scan_results.dat"
    if not scan_file.is_file():
        raise SystemExit(f"Could not find {scan_file}")

    script_dir = Path(__file__).resolve().parent
    outplots_dir = script_dir / "outplots"
    outplots_dir.mkdir(parents=True, exist_ok=True)

    rows = load_scan_results(scan_file)
    x_values = np.array([row[x_variable] for row in rows])
    y_values = np.array([row["Omega"] for row in rows])
    color_values = np.array([row[color_variable] for row in rows])

    # Linear scale plot
    figure, axis = plt.subplots(figsize=(8, 6))
    scatter = axis.scatter(x_values, y_values, c=color_values, s=28, cmap=args.cmap, alpha=0.7)
    axis.axhline(RELIC_LIMIT, color="red", linestyle="--", linewidth=1.5, label=f"Relic density upper limit = {RELIC_LIMIT}")
    axis.set_xlabel(x_variable)
    axis.set_ylabel(r'$\Omega h^2$')
    axis.set_title(f"Relic density vs {x_variable} (colored by {color_variable})")
    axis.grid(True, alpha=0.3)
    axis.legend()
    cbar = figure.colorbar(scatter, ax=axis)
    cbar.set_label(color_variable)
    figure.tight_layout()

    output_path = (
        Path(args.output)
        if args.output
        else outplots_dir / f"omega_vs_{x_variable}_colored_by_{color_variable}.png"
    )
    figure.savefig(output_path, dpi=200)
    print(f"Saved plot to {output_path}")

    # Log scale plot
    log_y_values = np.log10(y_values)
    figure, axis = plt.subplots(figsize=(8, 6))
    scatter = axis.scatter(x_values, log_y_values, c=color_values, s=28, cmap=args.cmap, alpha=0.7)
    axis.axhline(np.log10(RELIC_LIMIT), color="red", linestyle="--", linewidth=1.5, label=f"Relic density upper limit = {RELIC_LIMIT}")
    axis.set_xlabel(x_variable)
    axis.set_ylabel(r'$\log_{10}(\Omega h^2)$')
    axis.set_title(f"Relic density vs {x_variable} (colored by {color_variable})")
    axis.grid(True, alpha=0.3)
    axis.legend()
    cbar = figure.colorbar(scatter, ax=axis)
    cbar.set_label(color_variable)
    figure.tight_layout()

    log_output_path = (
        Path(args.output)
        if args.output
        else outplots_dir / f"log10_omega_vs_{x_variable}_colored_by_{color_variable}.png"
    )
    figure.savefig(log_output_path, dpi=200)
    print(f"Saved plot to {log_output_path}")

    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()