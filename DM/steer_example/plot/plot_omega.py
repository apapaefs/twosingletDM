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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Omega against a chosen variable from scan_results.dat."
    )
    parser.add_argument(
        "output_dir",
        help="Directory containing scan_results.dat, typically dmonly/output.",
    )
    parser.add_argument(
        "variable",
        help="Variable to plot on the x-axis. Choices: "
        + ", ".join(name for name in COLUMNS if name != "Omega"),
    )
    parser.add_argument(
        "--output",
        help="Optional output image path. Defaults to <output_dir>/omega_vs_<variable>.png",
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


def main():
    args = parse_args()

    variable = args.variable
    if variable not in COLUMNS:
        raise SystemExit(
            f"Unknown variable '{variable}'. Valid choices: {', '.join(COLUMNS)}"
        )
    if variable == "Omega":
        raise SystemExit("Use a variable other than Omega for the x-axis.")

    output_dir = Path(args.output_dir)
    scan_file = output_dir / "scan_results.dat"
    if not scan_file.is_file():
        raise SystemExit(f"Could not find {scan_file}")

    script_dir = Path(__file__).resolve().parent
    outplots_dir = script_dir / "outplots"
    outplots_dir.mkdir(parents=True, exist_ok=True)

    rows = load_scan_results(scan_file)
    x_values = [row[variable] for row in rows]
    y_values = [row["Omega"] for row in rows]
    log_y_values = [np.log10(row["Omega"]) for row in rows]

    figure, axis = plt.subplots(figsize=(7, 5))
    axis.scatter(x_values, y_values, s=28)
    axis.set_xlabel(variable)
    axis.set_ylabel("Omega")
    axis.set_title(f"Omega vs {variable}")
    axis.grid(True, alpha=0.3)
    figure.tight_layout()

    output_path = (
        Path(args.output)
        if args.output
        else outplots_dir / f"omega_vs_{variable}.png"
    )
    figure.savefig(output_path, dpi=200)

    if args.show:
        plt.show()
    else:
        plt.close(figure)


    figure, axis = plt.subplots(figsize=(7, 5))
    axis.scatter(x_values, log_y_values, s=28)
    axis.set_xlabel(variable)
    axis.set_ylabel("log10(Omega)")
    axis.set_title(f"Omega vs {variable}")
    axis.grid(True, alpha=0.3)
    figure.tight_layout()

    output_path = (
        Path(args.output)
        if args.output
        else outplots_dir / f"log10_omega_vs_{variable}.png"
    )
    figure.savefig(output_path, dpi=200)


    print(f"Saved plots to {output_path}")


if __name__ == "__main__":
    main()
