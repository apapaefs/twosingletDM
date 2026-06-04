#!/usr/bin/env python3

import argparse
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


RELIC_UPPER_LIMIT = 0.1224
DEFAULT_LIMIT_MODEL = "legacy-output"

CARD_FIELDS = [
    ("LX", "lx"),
    ("LHX", "lhx"),
    ("LSX", "lsx"),
    ("MX", "mx"),
    ("vevs", "vevs"),
    ("SinT", "sint"),
    ("Mh2", "mh2"),
]

NUMBER_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


@dataclass(frozen=True)
class PointInput:
    index: int
    lx: float
    lhx: float
    lsx: float
    mx: float
    vevs: float
    sint: float
    mh2: float


@dataclass(frozen=True)
class MicromegasResult:
    mdm: float
    omega: float
    dir_det: float


@dataclass(frozen=True)
class PointSummary:
    point: PointInput
    result: MicromegasResult
    dir_det_limit: float
    lux_base_limit: float
    relic_excluded: bool
    direct_detection_excluded: bool

    @property
    def dm_excluded(self):
        return self.relic_excluded or self.direct_detection_excluded


def default_micromegas_main():
    script_dir = Path(__file__).resolve().parent
    return script_dir.parent / "micromegas_6.1.15" / "TRSM" / "main"


def default_output_dir():
    return Path(__file__).resolve().parent / "single_point_output"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run or post-process one two-singlet DM micrOMEGAs point. "
            "This is the single-point analogue of run/MOrun.sh plus "
            "source/mO_excluder.cpp."
        )
    )
    parser.add_argument(
        "--card",
        type=Path,
        help="Existing single-point micrOMEGAs card, for example run/cards/MO_inp100.dat.",
    )
    parser.add_argument(
        "--index",
        type=int,
        help="Point index. Defaults to the number in MO_inp<index>.dat, or 1 for explicit parameters.",
    )
    parser.add_argument(
        "--micromegas-main",
        type=Path,
        default=default_micromegas_main(),
        help="Path to the micrOMEGAs executable. Defaults to the TRSM/main executable in this tree.",
    )
    parser.add_argument(
        "--micromegas-output",
        type=Path,
        help=(
            "Read an existing micrOMEGAs text output instead of running the executable. "
            "Useful for checking an existing OUT_mO_<index> file."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir(),
        help="Directory for the single-point outputs. Defaults to steer_example/single_point_output.",
    )
    parser.add_argument(
        "--no-rescale",
        action="store_true",
        help="Do not rescale the direct-detection limit by the relic-density fraction.",
    )
    parser.add_argument(
        "--relic-upper-limit",
        type=float,
        default=RELIC_UPPER_LIMIT,
        help=f"Upper relic-density limit. Defaults to {RELIC_UPPER_LIMIT}.",
    )
    parser.add_argument(
        "--limit-model",
        choices=["legacy-output", "lz2025-source"],
        default=DEFAULT_LIMIT_MODEL,
        help=(
            "Direct-detection limit fit. 'legacy-output' reproduces the existing "
            "steer_example/output files; 'lz2025-source' follows the active "
            "DirDetexcl branch in source/mO_excluder.cpp."
        ),
    )

    for card_key, attr_name in CARD_FIELDS:
        parser.add_argument(
            f"--{attr_name}",
            type=float,
            help=f"{card_key} value for an explicit single point when --card is not used.",
        )

    return parser.parse_args()


def infer_index_from_card(card_path):
    match = re.search(r"MO_inp(\d+)\.dat$", card_path.name)
    return int(match.group(1)) if match else None


def read_card(card_path, index=None):
    values = {}
    with card_path.open("r", encoding="ascii") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                values[parts[0]] = float(parts[1])

    missing = [card_key for card_key, _ in CARD_FIELDS if card_key not in values]
    if missing:
        raise ValueError(f"{card_path} is missing required card fields: {', '.join(missing)}")

    resolved_index = index if index is not None else infer_index_from_card(card_path)
    if resolved_index is None:
        raise ValueError(f"Could not infer an index from {card_path}; pass --index.")

    return PointInput(
        index=resolved_index,
        lx=values["LX"],
        lhx=values["LHX"],
        lsx=values["LSX"],
        mx=values["MX"],
        vevs=values["vevs"],
        sint=values["SinT"],
        mh2=values["Mh2"],
    )


def point_from_explicit_args(args):
    missing = [
        attr_name
        for _, attr_name in CARD_FIELDS
        if getattr(args, attr_name) is None
    ]
    if missing:
        raise SystemExit(
            "Pass --card, or provide all explicit point values: "
            + ", ".join(f"--{name}" for name in missing)
        )

    return PointInput(
        index=args.index if args.index is not None else 1,
        lx=args.lx,
        lhx=args.lhx,
        lsx=args.lsx,
        mx=args.mx,
        vevs=args.vevs,
        sint=args.sint,
        mh2=args.mh2,
    )


def write_card(point, card_path):
    lines = [
        f"LX\t\t{format_value(point.lx)}",
        f"LHX\t\t{format_value(point.lhx)}",
        f"LSX\t\t{format_value(point.lsx)}",
        f"MX\t\t{format_value(point.mx)}",
        f"vevs\t{format_value(point.vevs)}",
        f"SinT\t{format_value(point.sint)}",
        f"Mh2\t\t{format_value(point.mh2)}",
    ]
    card_path.write_text("\n".join(lines) + "\n", encoding="ascii")


def find_number(pattern, text, label):
    match = re.search(pattern, text, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"Could not find {label} in micrOMEGAs output")
    return float(match.group(1))


def parse_micromegas_output(text):
    # These regexes mirror the awk snippets in run/MOrun.sh.
    mdm = find_number(
        rf"(?:^|\s)(?:MHX|MX)\s*=\s*({NUMBER_PATTERN})",
        text,
        "dark matter mass",
    )
    omega = find_number(rf"Omega=({NUMBER_PATTERN})", text, "Omega")

    neutron_match = re.search(
        rf"^[ \t]*neutron[ \t]+SI[ \t]+({NUMBER_PATTERN})\b",
        text,
        flags=re.MULTILINE,
    )
    if not neutron_match:
        raise ValueError("Could not find neutron SI direct-detection cross section")

    return MicromegasResult(
        mdm=mdm,
        omega=omega,
        dir_det=float(neutron_match.group(1)),
    )


def direct_detection_base_limit(mdm, model=DEFAULT_LIMIT_MODEL):
    if mdm <= 0.0:
        raise ValueError("Dark matter mass must be positive")

    if model == "legacy-output":
        return legacy_output_direct_detection_base_limit(mdm)
    if model == "lz2025-source":
        return lz2025_source_direct_detection_base_limit(mdm)
    raise ValueError(f"Unknown direct-detection limit model: {model}")


def legacy_output_direct_detection_base_limit(mdm):
    """Older DirDetexcl fit used by the checked-in steer_example/output files."""
    if mdm < 10.0:
        result = math.exp(-23.2949 + 2.2493 * (-3.60228 + math.log(mdm)) ** 2)
        return result * (1.0 - 0.465188 + 1.64349 * (-2.3605 + math.log(mdm)) ** 2)
    if mdm < 11.0:
        return 8.75459e-9 - 6.96272e-10 * mdm
    if mdm < 18.0:
        result = math.exp(-23.2938 + 0.0135204 * (-24.8742 + mdm) ** 2)
        return result * (1.0 - 0.0826054 + 0.0110134 * (-13.3867 + mdm) ** 2)
    if mdm < 20.0:
        return 5.0606e-10 - 1.91361e-11 * mdm
    if mdm < 38.0:
        return 7.64124e-11 + 1.34496e-13 * (-36.6817 + mdm) ** 2
    if mdm < 60.0:
        result = 4.84775e-11 + 7.28365e-13 * mdm
    elif mdm < 100.0:
        result = math.exp(0.671326 * math.log(mdm) - 25.8196)
    else:
        result = 1.49167e-12 * math.exp(0.977674 * math.log(mdm))

    if 40.0 < mdm < 100.0:
        result = (0.0774 * mdm + 1.382) * 1.0e-11
    if 50.0 < mdm < 70.0:
        result = (0.0767 * mdm + 1.423) * 1.0e-11
    elif 100.0 < mdm < 200.0:
        result = (0.0805 * mdm + 1.235) * 1.0e-11
    elif 200.0 < mdm < 400.0:
        result = 8.135e-13 * mdm + 8.8e-12
    elif mdm > 400.0:
        result = 8.168e-13 * mdm + 7.5e-12

    return result


def lz2025_source_direct_detection_base_limit(mdm):
    """Active DirDetexcl branch in the checked-in source/mO_excluder.cpp."""
    if mdm < 10.0:
        result = math.exp(-23.2949 + 2.2493 * (-3.60228 + math.log(mdm)) ** 2)
        return result * (1.0 - 0.465188 + 1.64349 * (-2.3605 + math.log(mdm)) ** 2)
    if mdm < 40.0:
        exponent = (
            1.65
            / (1.0 - math.log10(40.0)) ** 2
            * (math.log10(mdm) - math.log10(40.0)) ** 2
        )
        return 2.45e-12 * 10.0 ** exponent
    if mdm < 60.0:
        return 2.32775e-12 + 2.28365e-15 * mdm
    if mdm < 80.0:
        return 1.68775e-12 + 1.28365e-14 * mdm
    if mdm < 100.0:
        return math.exp(0.671326 * math.log(mdm) - 29.5696)
    return 3.49167e-14 * math.exp(0.977674 * math.log(mdm))


def build_summary(
    point,
    result,
    relic_upper_limit=RELIC_UPPER_LIMIT,
    rescale=True,
    limit_model=DEFAULT_LIMIT_MODEL,
):
    lux_base_limit = direct_detection_base_limit(result.mdm, model=limit_model)
    if rescale:
        if result.omega > 0.0:
            dir_det_limit = lux_base_limit * relic_upper_limit / result.omega
        else:
            dir_det_limit = math.inf
    else:
        dir_det_limit = lux_base_limit

    return PointSummary(
        point=point,
        result=result,
        dir_det_limit=dir_det_limit,
        lux_base_limit=lux_base_limit,
        relic_excluded=result.omega > relic_upper_limit,
        direct_detection_excluded=result.dir_det > dir_det_limit,
    )


def run_micromegas(micromegas_main, card_path):
    if not micromegas_main.is_file():
        raise FileNotFoundError(f"micrOMEGAs executable not found: {micromegas_main}")
    if not micromegas_main.exists():
        raise FileNotFoundError(f"micrOMEGAs executable not found: {micromegas_main}")

    completed = subprocess.run(
        [str(micromegas_main), str(card_path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return completed.stdout


def format_value(value):
    if math.isinf(value):
        return "inf"
    return f"{value:.6g}"


def summary_line(summary):
    point = summary.point
    result = summary.result
    values = [
        point.index,
        point.lx,
        point.lhx,
        point.lsx,
        point.mx,
        point.vevs,
        point.sint,
        point.mh2,
        result.mdm,
        result.omega,
        result.dir_det,
        summary.dir_det_limit,
        summary.lux_base_limit,
    ]
    return "\t".join(format_value(value) if isinstance(value, float) else str(value) for value in values)


def scan_result_line(point, result):
    values = [
        point.index,
        point.lx,
        point.lhx,
        point.lsx,
        point.mx,
        point.vevs,
        point.sint,
        point.mh2,
        result.mdm,
        result.omega,
        result.dir_det,
    ]
    return "\t".join(format_value(value) if isinstance(value, float) else str(value) for value in values)


def write_marker(marker_path):
    marker_path.write_text("", encoding="ascii")


def write_outputs(output_dir, raw_output, summary):
    output_dir.mkdir(parents=True, exist_ok=True)
    index = summary.point.index

    raw_path = output_dir / f"OUT_mO_{index}"
    dm_data_path = output_dir / f"DM_data_{index}"
    scan_result_path = output_dir / f"scan_result_{index}.dat"

    raw_path.write_text(raw_output, encoding="ascii")
    dm_data_path.write_text(summary_line(summary) + "\n", encoding="ascii")
    scan_result_path.write_text(scan_result_line(summary.point, summary.result) + "\n", encoding="ascii")

    # Refresh only this point's markers so re-running the same point cannot leave
    # stale exclusion flags from a previous single-point attempt.
    markers = [
        output_dir / f"DM_EXCLUDED_{index}",
        output_dir / f"RelDens_EXCLUDED_{index}",
        output_dir / f"DirDet_EXCLUDED_{index}",
    ]
    for marker in markers:
        marker.unlink(missing_ok=True)

    if summary.dm_excluded:
        write_marker(output_dir / f"DM_EXCLUDED_{index}")
    if summary.relic_excluded:
        write_marker(output_dir / f"RelDens_EXCLUDED_{index}")
    if summary.direct_detection_excluded:
        write_marker(output_dir / f"DirDet_EXCLUDED_{index}")

    return raw_path, dm_data_path, scan_result_path


def print_summary(summary, paths):
    raw_path, dm_data_path, scan_result_path = paths
    status = "excluded" if summary.dm_excluded else "accepted"
    print(f"Point {summary.point.index}: {status}")
    print(f"  Omega: {format_value(summary.result.omega)}")
    print(f"  DirDet: {format_value(summary.result.dir_det)}")
    print(f"  DirDetLimit: {format_value(summary.dir_det_limit)}")
    print(f"  LUXBaseLimit: {format_value(summary.lux_base_limit)}")
    print(f"  Raw output: {raw_path}")
    print(f"  DM summary: {dm_data_path}")
    print(f"  Scan-style summary: {scan_result_path}")


def main():
    args = parse_args()

    if args.card is not None:
        point = read_card(args.card, index=args.index)
        card_path = args.card
    else:
        point = point_from_explicit_args(args)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        card_path = args.output_dir / f"MO_inp{point.index}.dat"
        write_card(point, card_path)

    if args.micromegas_output is not None:
        raw_output = args.micromegas_output.read_text(encoding="ascii")
    else:
        raw_output = run_micromegas(args.micromegas_main, card_path)

    result = parse_micromegas_output(raw_output)
    summary = build_summary(
        point,
        result,
        relic_upper_limit=args.relic_upper_limit,
        rescale=not args.no_rescale,
        limit_model=args.limit_model,
    )
    paths = write_outputs(args.output_dir, raw_output, summary)
    print_summary(summary, paths)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"run_single_point.py: {exc}", file=sys.stderr)
        sys.exit(1)
