#!/usr/bin/env python3

import math
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


__test__ = False

RELIC_UPPER_LIMIT = 0.1224
DEFAULT_LIMIT_MODEL = "legacy-output"
NUMBER_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


@dataclass(frozen=True)
class DMPoint:
    lX: float
    lPhiX: float
    lSX: float
    M3: float
    vs: float
    a12: float
    M2: float

    @property
    def sinT(self):
        return math.sin(self.a12)


@dataclass(frozen=True)
class MicromegasResult:
    mdm: float
    omega: float
    dir_det: float


@dataclass(frozen=True)
class DMSummary:
    point: DMPoint
    result: MicromegasResult
    dir_det_limit: float
    lux_base_limit: float
    relic_excluded: bool
    direct_detection_excluded: bool

    @property
    def passed(self):
        return not (self.relic_excluded or self.direct_detection_excluded)


def default_micromegas_main():
    return Path(__file__).resolve().parents[1] / "micromegas_6.1.15" / "TRSM" / "main"


def format_value(value):
    if math.isinf(value):
        return "inf"
    return f"{value:.6g}"


def write_micromegas_card(point, card_path):
    # This card is the one-point equivalent of steer_example/run/cards/MO_inp*.dat.
    # In the vx=0 setup, h3 is the Z2-odd dark matter scalar, so M3 maps to MX.
    lines = [
        f"LX\t\t{format_value(point.lX)}",
        f"LHX\t\t{format_value(point.lPhiX)}",
        f"LSX\t\t{format_value(point.lSX)}",
        f"MX\t\t{format_value(point.M3)}",
        f"vevs\t{format_value(point.vs)}",
        f"SinT\t{format_value(point.sinT)}",
        f"Mh2\t\t{format_value(point.M2)}",
    ]
    card_path.write_text("\n".join(lines) + "\n", encoding="ascii")


def find_number(pattern, text, label):
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"Could not find {label} in micrOMEGAs output")
    return float(match.group(1))


def parse_micromegas_output(text):
    """Extract the same three quantities that MOrun.sh parsed with awk."""
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
    if neutron_match is None:
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
    # Older DirDetexcl fit used by the checked-in steer_example/output files.
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
    # Active DirDetexcl branch in DM/steer_example/source/mO_excluder.cpp.
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


def summarize_dm_result(
    point,
    result,
    relic_upper_limit=RELIC_UPPER_LIMIT,
    limit_model=DEFAULT_LIMIT_MODEL,
    rescale=True,
):
    lux_base_limit = direct_detection_base_limit(result.mdm, model=limit_model)
    if rescale:
        if result.omega > 0.0:
            dir_det_limit = lux_base_limit * relic_upper_limit / result.omega
        else:
            dir_det_limit = math.inf
    else:
        dir_det_limit = lux_base_limit

    return DMSummary(
        point=point,
        result=result,
        dir_det_limit=dir_det_limit,
        lux_base_limit=lux_base_limit,
        relic_excluded=result.omega > relic_upper_limit,
        direct_detection_excluded=result.dir_det > dir_det_limit,
    )


def run_micromegas(card_path, micromegas_main):
    completed = subprocess.run(
        [str(micromegas_main), str(card_path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return completed.stdout


def dm_info_string(summary):
    reasons = []
    if summary.relic_excluded:
        reasons.append("Omega above relic-density upper limit")
    if summary.direct_detection_excluded:
        reasons.append("DirDet above rescaled direct-detection limit")
    if not reasons:
        reasons.append("all DM checks passed")

    return (
        f"DM check: {'Pass' if summary.passed else 'Fail'}\n"
        f"  LX={format_value(summary.point.lX)} "
        f"LHX={format_value(summary.point.lPhiX)} "
        f"LSX={format_value(summary.point.lSX)} "
        f"MX={format_value(summary.point.M3)} "
        f"vevs={format_value(summary.point.vs)} "
        f"SinT={format_value(summary.point.sinT)} "
        f"Mh2={format_value(summary.point.M2)}\n"
        f"  MDM={format_value(summary.result.mdm)} "
        f"Omega={format_value(summary.result.omega)} "
        f"DirDet={format_value(summary.result.dir_det)} "
        f"DirDetLimit={format_value(summary.dir_det_limit)} "
        f"LUXBaseLimit={format_value(summary.lux_base_limit)}\n"
        f"  Reason: {', '.join(reasons)}"
    )


def dm_exclusion_info(summary, relic_upper_limit, limit_model, rescale):
    return {
        "dm_mdm": summary.result.mdm,
        "dm_omega": summary.result.omega,
        "dm_relic_upper_limit": relic_upper_limit,
        "dm_dir_det": summary.result.dir_det,
        "dm_dir_det_limit": summary.dir_det_limit,
        "dm_lux_base_limit": summary.lux_base_limit,
        "dm_relic_excluded": summary.relic_excluded,
        "dm_direct_detection_excluded": summary.direct_detection_excluded,
        "dm_limit_model": limit_model,
        "dm_rescale": rescale,
    }


def empty_dm_exclusion_info():
    return {
        "dm_mdm": None,
        "dm_omega": None,
        "dm_relic_upper_limit": None,
        "dm_dir_det": None,
        "dm_dir_det_limit": None,
        "dm_lux_base_limit": None,
        "dm_relic_excluded": None,
        "dm_direct_detection_excluded": None,
        "dm_limit_model": None,
        "dm_rescale": None,
    }


def print_dm_info(info):
    lines = str(info).strip().splitlines()
    if not lines:
        lines = ["No DM information available"]

    title = "Dark Matter Check"
    width = max([len(title)] + [len(line) for line in lines])
    border = "+" + "-" * (width + 2) + "+"

    print(border)
    print(f"| {title.center(width)} |")
    print(border)
    for line in lines:
        print(f"| {line.ljust(width)} |")
    print(border)


def test_dm(
    lX,
    lPhiX,
    lSX,
    M3,
    vs,
    a12,
    M2,
    *,
    point_index=1,
    micromegas_main=None,
    output_dir=None,
    raw_output=None,
    relic_upper_limit=RELIC_UPPER_LIMIT,
    limit_model=DEFAULT_LIMIT_MODEL,
    rescale=True,
):
    """Run one vx=0 TRSM dark-matter point through micrOMEGAs.

    Returns `(passed, info, exclusion_info)`, where `passed` is true only if the
    relic-density and direct-detection checks pass. `info` is suitable for
    printing in debug mode, and `exclusion_info` contains the numerical values
    used in the DM exclusion.
    """
    point = DMPoint(
        lX=float(lX),
        lPhiX=float(lPhiX),
        lSX=float(lSX),
        M3=float(M3),
        vs=float(vs),
        a12=float(a12),
        M2=float(M2),
    )

    try:
        if raw_output is None:
            main_path = Path(micromegas_main) if micromegas_main else default_micromegas_main()
            if not main_path.is_file():
                raise FileNotFoundError(f"micrOMEGAs executable not found: {main_path}")

            if output_dir is None:
                with tempfile.TemporaryDirectory(prefix="trsm_dm_") as tmpdir:
                    card_path = Path(tmpdir) / f"MO_inp{point_index}.dat"
                    write_micromegas_card(point, card_path)
                    raw_output = run_micromegas(card_path, main_path)
            else:
                out_path = Path(output_dir)
                out_path.mkdir(parents=True, exist_ok=True)
                card_path = out_path / f"MO_inp{point_index}.dat"
                raw_path = out_path / f"OUT_mO_{point_index}"
                write_micromegas_card(point, card_path)
                raw_output = run_micromegas(card_path, main_path)
                raw_path.write_text(raw_output, encoding="ascii")

        result = parse_micromegas_output(raw_output)
        summary = summarize_dm_result(
            point,
            result,
            relic_upper_limit=relic_upper_limit,
            limit_model=limit_model,
            rescale=rescale,
        )
        return (
            summary.passed,
            dm_info_string(summary),
            dm_exclusion_info(summary, relic_upper_limit, limit_model, rescale),
        )

    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        return False, f"DM check: Error\n  Reason: {exc}", empty_dm_exclusion_info()


if __name__ == "__main__":
    passed, info, _ = test_dm(0.2, 0.1, 0.1, 50.0, 500.0, math.asin(0.3), 100.0)
    print_dm_info(info)
    raise SystemExit(0 if passed else 1)
