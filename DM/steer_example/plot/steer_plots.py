#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path


SCAN_COLUMNS = [
    "index",
    "LX",
    "LHX",
    "LSX",
    "MX",
    "vevs",
    "SinT",
    "Mh2",
]

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


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Infer plotting variables from run/oks.dat and run the plotting scripts."
        )
    )
    parser.add_argument(
        "--outdir",
        help=(
            "Directory containing scan_results.dat and related outputs. "
            "Defaults to ../output relative to this script."
        ),
    )
    parser.add_argument(
        "--xvar",
        help=(
            "Variable for the x-axis. If omitted, MX is preferred when it varies "
            "in run/oks.dat."
        ),
    )
    parser.add_argument(
        "--yvar",
        help=(
            "Secondary variable used as the y-axis for 2D plots and as the color "
            "variable for omega-colored plots."
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Pass --show to the underlying plotting scripts.",
    )
    return parser.parse_args()


def load_oks_rows(oks_file: Path):
    rows = []
    with oks_file.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) != len(SCAN_COLUMNS):
                raise ValueError(
                    f"{oks_file}:{line_number}: expected {len(SCAN_COLUMNS)} columns, "
                    f"found {len(parts)}"
                )
            row = {name: float(value) for name, value in zip(SCAN_COLUMNS, parts)}
            rows.append(row)
    if not rows:
        raise ValueError(f"{oks_file} is empty")
    return rows


def infer_varying_variables(rows):
    varying = []
    for name in SCAN_COLUMNS:
        if name == "index":
            continue
        values = {row[name] for row in rows}
        if len(values) > 1:
            varying.append(name)
    return varying


def validate_variable(name, option_name):
    if name not in PLOT_COLUMNS:
        raise SystemExit(
            f"Unknown {option_name} '{name}'. Valid choices: {', '.join(PLOT_COLUMNS)}"
        )


def choose_xvar(varying_variables, requested_xvar):
    if requested_xvar:
        validate_variable(requested_xvar, "x variable")
        if requested_xvar == "Omega":
            raise SystemExit("Omega cannot be used as --xvar.")
        return requested_xvar
    if "MX" in varying_variables:
        return "MX"
    if varying_variables:
        return varying_variables[0]
    raise SystemExit(
        "Could not infer an x variable from run/oks.dat because no scan variable varies. "
        "Use --xvar to choose one explicitly."
    )


def choose_yvar(varying_variables, xvar, requested_yvar):
    if requested_yvar:
        validate_variable(requested_yvar, "y variable")
        if requested_yvar == "Omega":
            raise SystemExit("Omega cannot be used as --yvar.")
        if requested_yvar == xvar:
            raise SystemExit("--xvar and --yvar must be different.")
        return requested_yvar
    for variable in varying_variables:
        if variable != xvar:
            return variable
    return None


def run_script(script_path: Path, outdir: Path, script_args, show=False):
    command = [sys.executable, str(script_path), str(outdir), *script_args]
    if show:
        command.append("--show")
    print("Running:", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def main():
    args = parse_args()

    plot_dir = Path(__file__).resolve().parent
    base_dir = plot_dir.parent
    outdir = Path(args.outdir).resolve() if args.outdir else (base_dir / "output")
    oks_file = base_dir / "run" / "oks.dat"

    if not oks_file.is_file():
        raise SystemExit(f"Could not find {oks_file}")
    if not outdir.is_dir():
        raise SystemExit(f"Could not find output directory {outdir}")
    if not (outdir / "scan_results.dat").is_file():
        raise SystemExit(f"Could not find {(outdir / 'scan_results.dat')}")

    oks_rows = load_oks_rows(oks_file)
    varying_variables = infer_varying_variables(oks_rows)
    xvar = choose_xvar(varying_variables, args.xvar)
    yvar = choose_yvar(varying_variables, xvar, args.yvar)

    print(f"Using output directory: {outdir}", flush=True)
    print(
        f"Detected varying scan variables: {', '.join(varying_variables) or 'none'}",
        flush=True,
    )
    print(f"Using x variable: {xvar}", flush=True)
    if yvar is not None:
        print(f"Using secondary variable: {yvar}", flush=True)
    else:
        print(
            "No secondary varying variable found; two-variable plots will be skipped.",
            flush=True,
        )

    run_script(plot_dir / "plot_omega.py", outdir, [xvar], show=args.show)

    if yvar is not None:
        run_script(
            plot_dir / "plot_omega_colored.py",
            outdir,
            [xvar, yvar],
            show=args.show,
        )
        run_script(
            plot_dir / "plot_relic_pass_2d.py",
            outdir,
            [xvar, yvar],
            show=args.show,
        )
        run_script(
            plot_dir / "plot_relic_strict_2d.py",
            outdir,
            [xvar, yvar],
            show=args.show,
        )

        accepted_file = outdir / "allall.dat"
        excluded_file = outdir / "dmexcl.dat"
        if accepted_file.is_file() and excluded_file.is_file():
            run_script(
                plot_dir / "plot_relic_pass_2d_with_excl.py",
                outdir,
                [xvar, yvar],
                show=args.show,
            )
        else:
            print(
                "Skipping plot_relic_pass_2d_with_excl.py because allall.dat or "
                "dmexcl.dat is missing.",
                flush=True,
            )


if __name__ == "__main__":
    main()
