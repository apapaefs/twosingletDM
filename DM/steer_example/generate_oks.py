#!/usr/bin/env python3

import argparse
import itertools
import math
from pathlib import Path

PARAMETERS = [
    ("lx", 0.2),
    ("lhx", 0.1),
    ("lsx", 0.1),
    ("mx", 800.0),
    ("vevs", 500.0),
    ("sint", 0.3),
    ("mh2", 500.0),
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate oks.dat for the updated write_mo.cpp layout."
    )
    parser.add_argument(
        "--output",
        default="oks.dat",
        help="Output file path. Defaults to oks.dat in the current directory.",
    )

    for name, default in PARAMETERS:
        label = name.upper() if name != "sint" else "SinT"
        parser.add_argument(
            f"--{name}-start",
            type=float,
            default=default,
            help=f"Minimum value for {label}. Defaults to {default}.",
        )
        parser.add_argument(
            f"--{name}-end",
            type=float,
            default=None,
            help=f"Maximum value for {label}. Defaults to the start value.",
        )
        parser.add_argument(
            f"--{name}-step",
            type=float,
            default=None,
            help=f"Step size for {label}. Required when start and end differ.",
        )

    for name, default in PARAMETERS:
        label = name.upper() if name != "sint" else "SinT"
        parser.add_argument(
            f"--{name}-log",
            action="store_true",
            help=f"Use logarithmic spacing for {label}."
        )
        parser.add_argument(
            f"--{name}-num-points",
            type=int,
            default=None,
            help=f"Number of points for {label} when using logarithmic spacing."
        )

    return parser.parse_args()


def build_scan_values(name, args):
    start = getattr(args, f"{name}_start")
    end = getattr(args, f"{name}_end")
    step = getattr(args, f"{name}_step")
    log = getattr(args, f"{name}_log")
    num_points = getattr(args, f"{name}_num_points")

    if log:
        if end is None:
            end = start

        if math.isclose(start, end, rel_tol=0.0, abs_tol=1e-12):
            return [start]

        if start <= 0 or end <= 0:
            raise ValueError(f"{name}: start and end must be positive for logarithmic spacing")

        if num_points is None or num_points < 1:
            raise ValueError(f"{name}: --{name}-num-points must be at least 1 for logarithmic spacing")

        log_start = math.log10(start)
        log_end = math.log10(end)
        if num_points == 1:
            return [start]
        step_log = (log_end - log_start) / (num_points - 1)
        values = [10 ** (log_start + i * step_log) for i in range(num_points)]
        return values

    else:
        if end is None:
            end = start

        if math.isclose(start, end, rel_tol=0.0, abs_tol=1e-12):
            return [start]

        if step is None or step == 0.0:
            raise ValueError(f"{name}: --{name}-step must be non-zero when start != end")

        if (end - start) * step < 0:
            raise ValueError(
                f"{name}: --{name}-step must move from start toward end"
            )

        values = []
        current = start
        tolerance = abs(step) * 1e-9 + 1e-12

        if step > 0:
            while current <= end + tolerance:
                values.append(current)
                current += step
        else:
            while current >= end - tolerance:
                values.append(current)
                current += step

        if not values:
            raise ValueError(f"{name}: no scan points generated")

        values[-1] = end
        return values


def main():
    args = parse_args()
    output_path = Path(args.output)
    scan_axes = [build_scan_values(name, args) for name, _ in PARAMETERS]

    with output_path.open("w", encoding="ascii") as stream:
        for index, values in enumerate(itertools.product(*scan_axes), start=1):
            stream.write(
                "{} {:.10g} {:.10g} {:.10g} {:.10g} {:.10g} {:.10g} {:.10g}\n".format(
                    index, *values
                )
            )


if __name__ == "__main__":
    main()
