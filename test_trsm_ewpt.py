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
DEFAULT_MINIMA_EXECUTABLE = DEFAULT_EXECUTABLE.with_name("MinimaTracer")
SM_GF = 1.1663787e-5
SM_VEV = math.sqrt(1 / math.sqrt(2) / SM_GF)
TRSM_COLUMNS = ["m1", "m2", "m3", "vs", "a12", "lx", "lphix", "lsx"]
SUMMARY_BASE_COLUMNS = [
    "status_nlo_stability",
    "status_ewsr",
    "status_tracing",
    "status_coex_pairs",
    "transition_history",
]
FOPT_TEMPERATURE_KINDS = ("crit", "nucl_approx", "nucl", "perc", "compl")


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
    minima_executable: Path | None = None
    model: str = "trsm"
    firstline: int = 2
    lastline: int = 2
    multistepmode: str = "default"
    use_multithreading: bool = True
    thigh: float = 300.0
    plot_phases: bool = False
    plot_output: Path | None = None
    plot_format: str = "both"
    sym_threshold: float = 1.0
    w1_threshold: float = 5.0
    wx_threshold: float = 1.0
    ws_threshold: float = 1.0


@dataclass(frozen=True)
class EWPTResult:
    rows: list
    input_path: Path
    output_path: Path
    command: list
    returncode: int
    stdout: str = ""
    stderr: str = ""
    minima_output_path: Path | None = None
    minima_command: list | None = None
    minima_returncode: int | None = None
    minima_stdout: str = ""
    minima_stderr: str = ""
    minima_analysis: object | None = None
    plot_paths: tuple = ()


@dataclass(frozen=True)
class Eq418Check:
    satisfied: bool
    conditions: dict
    couplings: dict
    lambdas: dict


@dataclass(frozen=True)
class FOPTStrength:
    transition_index: int
    temperature_kind: str
    status: str | None
    temperature: float
    false_vev: dict
    true_vev: dict
    delta_vev: dict
    ew_jump: float
    field_jump: float
    ew_jump_over_T: float
    ew_true_over_T: float
    ew_false_over_T: float
    field_jump_over_T: float

    def to_dict(self):
        return {
            "transition_index": self.transition_index,
            "temperature_kind": self.temperature_kind,
            "status": self.status,
            "temperature": self.temperature,
            "false_vev": self.false_vev,
            "true_vev": self.true_vev,
            "delta_vev": self.delta_vev,
            "ew_jump": self.ew_jump,
            "field_jump": self.field_jump,
            "ew_jump_over_T": self.ew_jump_over_T,
            "ew_true_over_T": self.ew_true_over_T,
            "ew_false_over_T": self.ew_false_over_T,
            "field_jump_over_T": self.field_jump_over_T,
        }


@dataclass(frozen=True)
class PhaseClassificationThresholds:
    sym_threshold: float = 1.0
    w1_threshold: float = 5.0
    wx_threshold: float = 1.0
    ws_threshold: float = 1.0


@dataclass(frozen=True)
class PhaseSample:
    temp: float
    w1: float
    wx: float
    ws: float
    veff: float

    def to_dict(self):
        return {
            "temp": self.temp,
            "w1": self.w1,
            "wx": self.wx,
            "ws": self.ws,
            "veff": self.veff,
        }


@dataclass(frozen=True)
class PhaseTrace:
    index: int
    samples: list

    def to_dict(self):
        return {
            "index": self.index,
            "samples": [sample.to_dict() for sample in self.samples],
        }


@dataclass(frozen=True)
class GlobalBranchPoint:
    temp: float
    phase_index: int
    label: str
    w1: float
    wx: float
    ws: float
    veff: float

    def to_dict(self):
        return {
            "temp": self.temp,
            "phase_index": self.phase_index,
            "label": self.label,
            "w1": self.w1,
            "wx": self.wx,
            "ws": self.ws,
            "veff": self.veff,
        }


@dataclass(frozen=True)
class MinimaTracerAnalysis:
    phase_traces: list
    temperature_grid: list
    global_branch: list
    global_phase_path: list
    ew_step_index: int | None

    def to_dict(self):
        return {
            "phase_traces": [trace.to_dict() for trace in self.phase_traces],
            "temperature_grid": self.temperature_grid,
            "global_branch": [point.to_dict() for point in self.global_branch],
            "global_phase_path": self.global_phase_path,
            "ew_step_index": self.ew_step_index,
        }


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
        f"--model={config.model}",
        f"--input={input_path}",
        f"--output={output_path}",
        f"--firstline={config.firstline}",
        f"--lastline={config.lastline}",
        f"--multistepmode={config.multistepmode}",
        f"--usemultithreading={str(config.use_multithreading).lower()}",
        f"--thigh={config.thigh}",
    ]


def build_minimatracer_command(config, input_path, output_prefix):
    minima_executable = config.minima_executable
    if minima_executable is None:
        minima_executable = Path(config.executable).with_name("MinimaTracer")
    return [
        str(minima_executable),
        f"--model={config.model}",
        f"--input={input_path}",
        f"--output={output_prefix}",
        f"--firstline={config.firstline}",
        f"--lastline={config.lastline}",
        f"--multistepmode={config.multistepmode}",
        f"--usemultithreading={str(config.use_multithreading).lower()}",
        f"--thigh={config.thigh}",
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


def is_finite_number(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def calctemps_transition_indices(row):
    indices = set()
    for column in row:
        for kind in FOPT_TEMPERATURE_KINDS:
            prefix = f"T_{kind}_"
            if column.startswith(prefix):
                suffix = column[len(prefix) :]
                if suffix.isdigit():
                    indices.add(int(suffix))
    return sorted(indices)


def calculate_fopt_strengths(row):
    strengths = []
    for transition_index in calctemps_transition_indices(row):
        for kind in FOPT_TEMPERATURE_KINDS:
            temperature = finite_float(row.get(f"T_{kind}_{transition_index}"))
            if temperature is None or temperature <= 0.0:
                continue

            status = row.get(f"status_{kind}_{transition_index}")
            if status is not None and str(status).lower() != "success":
                continue

            false_vev = {}
            true_vev = {}
            for component in ("w1", "wx", "ws"):
                false_value = finite_float(
                    row.get(f"{component}_{kind}_false_{transition_index}")
                )
                true_value = finite_float(
                    row.get(f"{component}_{kind}_true_{transition_index}")
                )
                if false_value is None or true_value is None:
                    break
                false_vev[component] = false_value
                true_vev[component] = true_value
            else:
                delta_vev = {
                    component: true_vev[component] - false_vev[component]
                    for component in ("w1", "wx", "ws")
                }
                ew_jump = abs(delta_vev["w1"])
                field_jump = math.sqrt(
                    sum(delta_vev[component] ** 2 for component in ("w1", "wx", "ws"))
                )
                strengths.append(
                    FOPTStrength(
                        transition_index=transition_index,
                        temperature_kind=kind,
                        status=status,
                        temperature=temperature,
                        false_vev=false_vev,
                        true_vev=true_vev,
                        delta_vev=delta_vev,
                        ew_jump=ew_jump,
                        field_jump=field_jump,
                        ew_jump_over_T=ew_jump / temperature,
                        ew_true_over_T=abs(true_vev["w1"]) / temperature,
                        ew_false_over_T=abs(false_vev["w1"]) / temperature,
                        field_jump_over_T=field_jump / temperature,
                    )
                )
    return strengths


def calctemps_summary_columns(row):
    columns = [
        column
        for column in SUMMARY_BASE_COLUMNS
        if column != "transition_history"
    ]
    for transition_index in calctemps_transition_indices(row):
        for kind in FOPT_TEMPERATURE_KINDS:
            columns.append(f"status_{kind}_{transition_index}")
            columns.append(f"T_{kind}_{transition_index}")
            for phase in ("false", "true"):
                for component in ("w1", "wx", "ws"):
                    columns.append(
                        f"{component}_{kind}_{phase}_{transition_index}"
                    )
            if kind == "crit":
                columns.append(f"status_bounce_sol_{transition_index}")
    columns.append("transition_history")
    seen = set()
    ordered_columns = []
    for column in columns:
        if column not in seen:
            ordered_columns.append(column)
            seen.add(column)
    return ordered_columns


def parse_minimatracer_output(path):
    path = Path(path)
    with path.open("r", encoding="ascii", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        try:
            header = normalize_header(next(reader))
        except StopIteration:
            return []

        phase_indices = []
        for column in header:
            if column.startswith("Temp_"):
                try:
                    phase_indices.append(int(column.split("_", 1)[1]))
                except ValueError:
                    continue
        phase_indices = sorted(set(phase_indices))
        samples_by_phase = {index: {} for index in phase_indices}

        for raw_row in reader:
            if not raw_row:
                continue
            if len(raw_row) < len(header):
                raw_row = raw_row + [""] * (len(header) - len(raw_row))
            row = {
                column: coerce_output_value(raw_row[index])
                for index, column in enumerate(header)
            }

            for phase_index in phase_indices:
                temp = row.get(f"Temp_{phase_index}")
                w1 = row.get(f"w1(Temp_{phase_index})")
                wx = row.get(f"wx(Temp_{phase_index})")
                ws = row.get(f"ws(Temp_{phase_index})")
                veff = row.get(f"Veff(Temp_{phase_index})")
                values = [temp, w1, wx, ws, veff]
                if not all(is_finite_number(value) for value in values):
                    continue
                sample = PhaseSample(
                    temp=float(temp),
                    w1=float(w1),
                    wx=float(wx),
                    ws=float(ws),
                    veff=float(veff),
                )
                previous = samples_by_phase[phase_index].get(sample.temp)
                if previous is None or sample.veff < previous.veff:
                    samples_by_phase[phase_index][sample.temp] = sample

        traces = []
        for phase_index in phase_indices:
            samples = [
                samples_by_phase[phase_index][temp]
                for temp in sorted(samples_by_phase[phase_index])
            ]
            if samples:
                traces.append(PhaseTrace(index=phase_index, samples=samples))
        return traces


def interpolate_value(low, high, fraction, attribute):
    return getattr(low, attribute) + fraction * (getattr(high, attribute) - getattr(low, attribute))


def interpolate_phase_at(trace, temp):
    samples = trace.samples
    if not samples:
        return None
    if temp < samples[0].temp or temp > samples[-1].temp:
        return None
    for sample in samples:
        if math.isclose(sample.temp, temp, rel_tol=0.0, abs_tol=1e-9):
            return sample
    for low, high in zip(samples, samples[1:]):
        if low.temp <= temp <= high.temp:
            if math.isclose(low.temp, high.temp, rel_tol=0.0, abs_tol=1e-12):
                return low
            fraction = (temp - low.temp) / (high.temp - low.temp)
            return PhaseSample(
                temp=temp,
                w1=interpolate_value(low, high, fraction, "w1"),
                wx=interpolate_value(low, high, fraction, "wx"),
                ws=interpolate_value(low, high, fraction, "ws"),
                veff=interpolate_value(low, high, fraction, "veff"),
            )
    return None


def classify_phase(w1, wx, ws, thresholds):
    w1_abs = abs(w1)
    wx_abs = abs(wx)
    ws_abs = abs(ws)
    if (
        w1_abs < thresholds.sym_threshold
        and wx_abs < thresholds.sym_threshold
        and ws_abs < thresholds.sym_threshold
    ):
        return "SYM"

    ew_broken = w1_abs >= thresholds.w1_threshold
    x_broken = wx_abs >= thresholds.wx_threshold
    s_broken = ws_abs >= thresholds.ws_threshold

    if ew_broken and x_broken:
        return "EW_X_BROKEN"
    if ew_broken:
        return "EW"
    if x_broken:
        return "X_BROKEN"
    if s_broken:
        return "SINGLET_S"
    return "MIXED"


def compress_labels_for_cooling(global_branch):
    compressed = []
    for point in sorted(global_branch, key=lambda item: item.temp, reverse=True):
        if not compressed or compressed[-1] != point.label:
            compressed.append(point.label)
    return compressed


def find_ew_step_index(global_phase_path):
    for index in range(1, len(global_phase_path)):
        if global_phase_path[index].startswith("EW"):
            return index - 1
    return None


def analyze_minimatracer_phases(phase_traces, thresholds):
    temperature_grid = sorted(
        {
            sample.temp
            for trace in phase_traces
            for sample in trace.samples
            if is_finite_number(sample.temp)
        }
    )
    global_branch = []
    for temp in temperature_grid:
        candidates = []
        for trace in phase_traces:
            sample = interpolate_phase_at(trace, temp)
            if sample is not None:
                candidates.append((trace.index, sample))
        if not candidates:
            continue
        phase_index, sample = min(candidates, key=lambda item: item[1].veff)
        label = classify_phase(sample.w1, sample.wx, sample.ws, thresholds)
        global_branch.append(
            GlobalBranchPoint(
                temp=temp,
                phase_index=phase_index,
                label=label,
                w1=sample.w1,
                wx=sample.wx,
                ws=sample.ws,
                veff=sample.veff,
            )
        )

    global_phase_path = compress_labels_for_cooling(global_branch)
    return MinimaTracerAnalysis(
        phase_traces=phase_traces,
        temperature_grid=temperature_grid,
        global_branch=global_branch,
        global_phase_path=global_phase_path,
        ew_step_index=find_ew_step_index(global_phase_path),
    )


def phase_plot_paths(output_base, plot_format):
    output_base = Path(output_base)
    root = output_base.with_suffix("") if output_base.suffix in {".png", ".pdf"} else output_base
    if plot_format == "png":
        return (output_base if output_base.suffix == ".png" else root.with_suffix(".png"),)
    if plot_format == "pdf":
        return (output_base if output_base.suffix == ".pdf" else root.with_suffix(".pdf"),)
    if plot_format == "both":
        return (root.with_suffix(".png"), root.with_suffix(".pdf"))
    raise ValueError("plot_format must be one of: png, pdf, both")


def plot_minimatracer_phases(phase_traces, output_base, plot_format="both", thigh=None):
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required when --plot-phases is used") from exc

    paths = phase_plot_paths(output_base, plot_format)
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    line_styles = {"w1": "-", "wx": "--", "ws": ":"}
    colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["C0"])
    for trace in phase_traces:
        color = colors[trace.index % len(colors)]
        temps = [sample.temp for sample in trace.samples]
        for component, line_style in line_styles.items():
            values = [getattr(sample, component) for sample in trace.samples]
            ax.plot(
                temps,
                values,
                linestyle=line_style,
                color=color,
                label=f"phase {trace.index} {component}",
            )
    ax.set_xlabel("T [GeV]")
    ax.set_ylabel("VEV [GeV]")
    if thigh is not None:
        ax.set_xlim((0.0, float(thigh)))
    ax.legend(fontsize="small", ncol=2)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    for path in paths:
        fig.savefig(path)
    plt.close(fig)
    return paths


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
    minima_output_prefix = run_dir / "minima_trace"
    minima_output_path = Path(str(minima_output_prefix) + "_1.tsv")
    write_trsm_input(input_path, point)

    minima_command = build_minimatracer_command(config, input_path, minima_output_prefix)
    minima_completed = runner(
        minima_command,
        check=True,
        capture_output=True,
        text=True,
        cwd=run_dir,
    )
    phase_traces = parse_minimatracer_output(minima_output_path)
    thresholds = PhaseClassificationThresholds(
        sym_threshold=config.sym_threshold,
        w1_threshold=config.w1_threshold,
        wx_threshold=config.wx_threshold,
        ws_threshold=config.ws_threshold,
    )
    minima_analysis = analyze_minimatracer_phases(phase_traces, thresholds)
    plot_paths = ()
    if config.plot_phases:
        plot_base = config.plot_output
        if plot_base is None:
            plot_base = run_dir / "minima_phases"
        else:
            plot_base = Path(plot_base).expanduser()
            if not plot_base.is_absolute():
                plot_base = run_dir / plot_base
        plot_paths = plot_minimatracer_phases(
            phase_traces,
            plot_base,
            plot_format=config.plot_format,
            thigh=config.thigh,
        )

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
        minima_output_path=minima_output_path,
        minima_command=minima_command,
        minima_returncode=minima_completed.returncode,
        minima_stdout=minima_completed.stdout,
        minima_stderr=minima_completed.stderr,
        minima_analysis=minima_analysis,
        plot_paths=tuple(plot_paths),
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
    if result.minima_output_path is not None:
        lines.append(f"minima_output: {result.minima_output_path}")
    if result.plot_paths:
        for plot_path in result.plot_paths:
            lines.append(f"phase_plot: {plot_path}")
    if result.minima_analysis is not None:
        path = " -> ".join(result.minima_analysis.global_phase_path)
        lines.append(f"global_phase_path: {path if path else 'not_set'}")
        lines.append(f"ew_step_index: {result.minima_analysis.ew_step_index}")
    for column in calctemps_summary_columns(row):
        if column in row:
            lines.append(f"{column}: {row[column]}")
    for strength in calculate_fopt_strengths(row):
        lines.append(
            f"fopt_strength_{strength.temperature_kind}_{strength.transition_index}: "
            f"T={strength.temperature:.6g}, "
            f"ew_jump/T={strength.ew_jump_over_T:.6g}, "
            f"ew_true/T={strength.ew_true_over_T:.6g}, "
            f"field_jump/T={strength.field_jump_over_T:.6g}"
        )
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


def path_to_string(path):
    if path is None:
        return None
    return str(path)


def result_to_json(result):
    row = first_row(result)
    return {
        "calctemps": row,
        "transition_strengths": [
            strength.to_dict() for strength in calculate_fopt_strengths(row)
        ],
        "paths": {
            "input": path_to_string(result.input_path),
            "calctemps_output": path_to_string(result.output_path),
            "minima_output": path_to_string(result.minima_output_path),
            "phase_plots": [str(path) for path in result.plot_paths],
        },
        "commands": {
            "minimatracer": result.minima_command,
            "calctemps": result.command,
        },
        "minimatracer": (
            result.minima_analysis.to_dict()
            if result.minima_analysis is not None
            else None
        ),
    }


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
    parser.add_argument("--minima-executable", type=Path)
    parser.add_argument("--thigh", type=float, default=300.0)
    parser.add_argument("--multistepmode", default="default")
    parser.add_argument("--plot-phases", action="store_true")
    parser.add_argument("--plot-output", type=Path)
    parser.add_argument("--plot-format", choices=["png", "pdf", "both"], default="both")
    parser.add_argument("--sym-threshold", type=float, default=1.0)
    parser.add_argument("--w1-threshold", type=float, default=5.0)
    parser.add_argument("--wx-threshold", type=float, default=1.0)
    parser.add_argument("--ws-threshold", type=float, default=1.0)
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
    config = EWPTConfig(
        executable=args.executable,
        minima_executable=args.minima_executable,
        multistepmode=args.multistepmode,
        thigh=args.thigh,
        plot_phases=args.plot_phases,
        plot_output=args.plot_output,
        plot_format=args.plot_format,
        sym_threshold=args.sym_threshold,
        w1_threshold=args.w1_threshold,
        wx_threshold=args.wx_threshold,
        ws_threshold=args.ws_threshold,
    )
    result = run_trsm_ewpt(
        point,
        config=config,
        workdir=args.workdir,
        keep_files=args.keep_files or args.workdir is not None,
    )

    if args.json:
        print(json.dumps(result_to_json(result), indent=2, sort_keys=True))
    else:
        print(summarize_result(result))


if __name__ == "__main__":
    main()
