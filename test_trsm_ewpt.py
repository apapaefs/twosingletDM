#!/usr/bin/env python3

import argparse
import csv
import json
import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


DEFAULT_EXECUTABLE = Path(
    "/Users/apapaefs/Projects/TwoSingletDM/BSMPT/build/macos-armv8-release/bin/CalcTemps"
)
SM_GF = 1.1663787e-5
SM_VEV = math.sqrt(1 / math.sqrt(2) / SM_GF)
TRSM_COLUMNS = ["m1", "m2", "m3", "vs", "a12", "lx", "lphix", "lsx"]
SUMMARY_COLUMNS = [
    "status_nlo_stability",
    "status_ewsr",
    "status_tracing",
    "status_coex_pairs",
    "status_crit_0",
    "T_crit_0",
    "w1_crit_false_0",
    "wx_crit_false_0",
    "ws_crit_false_0",
    "w1_crit_true_0",
    "wx_crit_true_0",
    "ws_crit_true_0",
    "status_nucl_0",
    "T_nucl_0",
    "w1_nucl_false_0",
    "wx_nucl_false_0",
    "ws_nucl_false_0",
    "w1_nucl_true_0",
    "wx_nucl_true_0",
    "ws_nucl_true_0",
    "transition_history",
]


@dataclass(frozen=True)
class TRSMEWPTPoint:
    m2: float
    m3: float
    vs: float
    a12: float
    lx: float
    lphix: float
    lsx: float
    m1: float = 125.09
    index: int = 1

    def as_row(self):
        return [
            self.index,
            self.m1,
            self.m2,
            self.m3,
            self.vs,
            self.a12,
            self.lx,
            self.lphix,
            self.lsx,
        ]


@dataclass(frozen=True)
class EWPTConfig:
    executable: Path = DEFAULT_EXECUTABLE
    model: str = "trsm"
    firstline: int = 2
    lastline: int = 2
    multistepmode: str = "2"
    use_multithreading: bool = True


@dataclass(frozen=True)
class EWPTResult:
    rows: list
    input_path: Path
    output_path: Path
    command: list
    returncode: int
    stdout: str = ""
    stderr: str = ""


@dataclass(frozen=True)
class Eq418Check:
    satisfied: bool
    conditions: dict
    couplings: dict
    lambdas: dict


def format_value(value):
    return str(value)


def write_trsm_input(path, point):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow([""] + TRSM_COLUMNS)
        writer.writerow([format_value(value) for value in point.as_row()])


def build_calctemps_command(config, input_path, output_path):
    return [
        str(config.executable),
        config.model,
        str(input_path),
        str(output_path),
        str(config.firstline),
        str(config.lastline),
        f"--multistepmode={config.multistepmode}",
        f"--usemultithreading={str(config.use_multithreading).lower()}",
    ]


def coerce_output_value(value):
    if value == "":
        return value
    try:
        return float(value)
    except ValueError:
        return value


def normalize_header(header):
    return ["index" if column == "" else column for column in header]


def parse_calctemps_output(path):
    path = Path(path)
    with path.open("r", encoding="ascii", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        try:
            header = normalize_header(next(reader))
        except StopIteration:
            return []

        rows = []
        for raw_row in reader:
            if not raw_row:
                continue
            if len(raw_row) < len(header):
                raw_row = raw_row + [""] * (len(header) - len(raw_row))
            row = {
                column: coerce_output_value(raw_row[index])
                for index, column in enumerate(header)
            }
            rows.append(row)
        return rows


def trsm_tree_lambdas(row, sm_vev=SM_VEV):
    m1 = float(row["m1"])
    m2 = float(row["m2"])
    vs = float(row["vs"])
    a12 = float(row["a12"])
    lx = float(row["lx"])
    lphix = float(row["lphix"])
    lsx = float(row["lsx"])

    cos_a12 = math.cos(a12)
    sin_a12 = math.sin(a12)
    lphi = (m1**2 * cos_a12**2 + m2**2 * sin_a12**2) / (2.0 * sm_vev**2)
    ls = (m2**2 * cos_a12**2 + m1**2 * sin_a12**2) / (2.0 * vs**2)
    lphis = ((-m1**2 + m2**2) * cos_a12 * sin_a12) / (sm_vev * vs)

    return {
        "lambda_phiphi": lphi,
        "lambda_xx": lx,
        "lambda_ss": ls,
        "lambda_phix": lphix,
        "lambda_phis": lphis,
        "lambda_xs": lsx,
    }


def check_eq_4_18(row):
    lambdas = trsm_tree_lambdas(row)
    couplings = {
        "c_phiphi": 0.5 * lambdas["lambda_phiphi"],
        "c_xx": 0.5 * lambdas["lambda_xx"],
        "c_ss": 0.5 * lambdas["lambda_ss"],
        "c_phix": 0.5 * lambdas["lambda_phix"],
        "c_phis": 0.5 * lambdas["lambda_phis"],
        "c_xs": 0.5 * lambdas["lambda_xs"],
    }

    c_phi = couplings["c_phiphi"]
    c_x = couplings["c_xx"]
    c_s = couplings["c_ss"]
    c_phix = couplings["c_phix"]
    c_phis = couplings["c_phis"]
    c_xs = couplings["c_xs"]

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

    conditions = {
        "c_phiphi > 0": c_phi > 0.0,
        "c_xx > 0": c_x > 0.0,
        "c_ss > 0": c_s > 0.0,
        "C_phix > 0": c_phix_minor > 0.0,
        "C_phis > 0": c_phis_minor > 0.0,
        "C_xs > 0": c_xs_minor > 0.0,
        "D > 0": determinant > 0.0,
    }
    couplings.update(
        {
            "C_phix": c_phix_minor,
            "C_phis": c_phis_minor,
            "C_xs": c_xs_minor,
            "D": determinant,
        }
    )
    return Eq418Check(
        satisfied=all(conditions.values()),
        conditions=conditions,
        couplings=couplings,
        lambdas=lambdas,
    )


def run_trsm_ewpt(
    point,
    config=None,
    workdir=None,
    keep_files=False,
    runner=subprocess.run,
):
    if config is None:
        config = EWPTConfig()

    temporary_directory = None
    if workdir is None:
        if keep_files:
            run_dir = Path(tempfile.mkdtemp(prefix="trsm-ewpt-"))
        else:
            temporary_directory = tempfile.TemporaryDirectory(prefix="trsm-ewpt-")
            run_dir = Path(temporary_directory.name)
    else:
        run_dir = Path(workdir).expanduser().resolve()
        run_dir.mkdir(parents=True, exist_ok=True)

    input_path = run_dir / "TRSM_Input.tsv"
    output_path = run_dir / "test.output_new.csv"
    write_trsm_input(input_path, point)

    command = build_calctemps_command(config, input_path, output_path)
    completed = runner(
        command,
        check=True,
        capture_output=True,
        text=True,
        cwd=run_dir,
    )
    rows = parse_calctemps_output(output_path)

    if temporary_directory is not None:
        temporary_directory.cleanup()

    return EWPTResult(
        rows=rows,
        input_path=input_path,
        output_path=output_path,
        command=command,
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
    )


def run_vxzero_scan_point(
    m2,
    m3,
    vs,
    a12,
    lX,
    lPhiX,
    lSX,
    m1=125.09,
    **kwargs,
):
    point = TRSMEWPTPoint(
        m1=m1,
        m2=m2,
        m3=m3,
        vs=vs,
        a12=a12,
        lx=lX,
        lphix=lPhiX,
        lsx=lSX,
    )
    return run_trsm_ewpt(point, **kwargs)


def first_row(result):
    if not result.rows:
        return {}
    return result.rows[0]


def summarize_result(result):
    row = first_row(result)
    lines = [f"output: {result.output_path}"]
    for column in SUMMARY_COLUMNS:
        if column in row:
            lines.append(f"{column}: {row[column]}")
    try:
        eq_4_18 = check_eq_4_18(row)
    except KeyError:
        eq_4_18 = None
    if eq_4_18 is not None:
        condition_summary = "; ".join(
            f"{name}={value}" for name, value in eq_4_18.conditions.items()
        )
        lines.append(f"eq_4_18_satisfied: {eq_4_18.satisfied}")
        lines.append(f"eq_4_18_conditions: {condition_summary}")
    return "\n".join(lines)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run BSMPT CalcTemps for one TRSM EWPT point and parse the TSV output."
    )
    parser.add_argument("--m1", type=float, default=125.09)
    parser.add_argument("--m2", type=float, required=True)
    parser.add_argument("--m3", type=float, required=True)
    parser.add_argument("--vs", type=float, required=True)
    parser.add_argument("--a12", type=float, required=True)
    parser.add_argument("--lx", type=float, required=True)
    parser.add_argument("--lphix", type=float, required=True)
    parser.add_argument("--lsx", type=float, required=True)
    parser.add_argument("--index", type=int, default=1)
    parser.add_argument("--executable", type=Path, default=DEFAULT_EXECUTABLE)
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--keep-files", action="store_true")
    parser.add_argument("--json", action="store_true", help="Print the parsed row as JSON.")
    return parser.parse_args()


def main():
    args = parse_args()
    point = TRSMEWPTPoint(
        index=args.index,
        m1=args.m1,
        m2=args.m2,
        m3=args.m3,
        vs=args.vs,
        a12=args.a12,
        lx=args.lx,
        lphix=args.lphix,
        lsx=args.lsx,
    )
    config = EWPTConfig(executable=args.executable)
    result = run_trsm_ewpt(
        point,
        config=config,
        workdir=args.workdir,
        keep_files=args.keep_files or args.workdir is not None,
    )

    if args.json:
        print(json.dumps(first_row(result), indent=2, sort_keys=True))
    else:
        print(summarize_result(result))


if __name__ == "__main__":
    main()
