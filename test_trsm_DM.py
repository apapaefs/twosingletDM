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


FERMI_LAT_R16_LINE_LIMITS = (
    (0.214, 45.5e-8),
    (0.234, 38.0e-8),
    (0.255, 34.0e-8),
    (0.278, 32.4e-8),
    (0.303, 30.9e-8),
    (0.329, 28.9e-8),
    (0.358, 26.3e-8),
    (0.388, 21.5e-8),
    (0.421, 18.8e-8),
    (0.456, 16.7e-8),
    (0.493, 14.7e-8),
    (0.533, 13.5e-8),
    (0.576, 11.8e-8),
    (0.620, 10.4e-8),
    (0.668, 9.84e-8),
    (0.718, 8.78e-8),
    (0.770, 7.80e-8),
    (0.826, 6.74e-8),
    (0.885, 6.36e-8),
    (0.947, 7.32e-8),
    (1.01, 8.68e-8),
    (1.08, 9.29e-8),
    (1.16, 8.68e-8),
    (1.24, 7.55e-8),
    (1.32, 6.39e-8),
    (1.41, 5.12e-8),
    (1.50, 4.03e-8),
    (1.60, 3.17e-8),
    (1.70, 2.92e-8),
    (1.81, 2.20e-8),
    (1.93, 1.84e-8),
    (2.06, 1.93e-8),
    (2.19, 1.87e-8),
    (2.33, 1.75e-8),
    (2.49, 1.35e-8),
    (2.65, 1.07e-8),
    (2.82, 0.844e-8),
    (3.00, 0.738e-8),
    (3.20, 0.780e-8),
    (3.40, 0.915e-8),
    (3.62, 1.02e-8),
    (3.85, 0.925e-8),
    (4.09, 0.764e-8),
    (4.35, 0.715e-8),
    (4.63, 0.561e-8),
    (4.91, 0.445e-8),
    (5.22, 0.263e-8),
    (5.54, 21.6e-10),
    (5.87, 20.6e-10),
    (6.23, 17.0e-10),
    (6.60, 22.5e-10),
    (6.99, 26.0e-10),
    (7.40, 33.4e-10),
    (7.83, 25.3e-10),
    (8.28, 18.6e-10),
    (8.76, 22.1e-10),
    (9.26, 12.7e-10),
    (9.79, 7.59e-10),
    (10.4, 8.14e-10),
    (10.9, 12.4e-10),
    (11.6, 16.0e-10),
    (12.2, 7.72e-10),
    (12.9, 7.83e-10),
    (13.6, 8.42e-10),
    (14.4, 5.91e-10),
    (15.2, 6.52e-10),
    (16.1, 6.66e-10),
    (17.0, 8.18e-10),
    (17.9, 10.2e-10),
    (18.9, 5.79e-10),
    (20.0, 3.65e-10),
    (21.1, 6.56e-10),
    (22.3, 3.66e-10),
    (23.6, 3.74e-10),
    (24.9, 2.97e-10),
    (26.4, 3.78e-10),
    (27.9, 4.56e-10),
    (29.5, 7.05e-10),
    (31.2, 4.37e-10),
    (33.0, 3.28e-10),
    (34.9, 4.17e-10),
    (36.9, 4.71e-10),
    (39.0, 3.18e-10),
    (41.3, 3.07e-10),
    (43.8, 4.71e-10),
    (46.4, 5.66e-10),
    (49.1, 6.40e-10),
    (52.1, 4.56e-10),
    (55.2, 3.96e-10),
    (58.6, 4.85e-10),
    (62.2, 3.32e-10),
    (66.0, 1.82e-10),
    (70.1, 1.90e-10),
    (74.5, 3.63e-10),
    (79.2, 1.48e-10),
    (84.2, 0.951e-10),
    (89.6, 0.947e-10),
    (95.4, 0.891e-10),
    (102.0, 2.29e-10),
    (108.0, 4.89e-10),
    (115.0, 4.92e-10),
    (123.0, 3.84e-10),
    (131.0, 3.11e-10),
    (140.0, 1.48e-10),
    (150.0, 0.765e-10),
    (160.0, 0.764e-10),
    (171.0, 1.11e-10),
    (183.0, 1.70e-10),
    (196.0, 2.22e-10),
    (210.0, 2.85e-10),
    (225.0, 1.59e-10),
    (241.0, 1.93e-10),
    (259.0, 0.867e-10),
    (276.0, 0.843e-10),
    (294.0, 1.32e-10),
    (321.0, 1.45e-10),
    (345.0, 1.17e-10),
    (367.0, 0.646e-10),
    (396.0, 0.613e-10),
    (427.0, 0.560e-10),
    (462.0, 0.487e-10),
)


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
class IndirectLineChannel:
    energy_gev: float
    flux_cm2_s: float


@dataclass(frozen=True)
class MicromegasResult:
    mdm: float
    omega: float
    dir_det: float
    indirect_line_channels: tuple[IndirectLineChannel, ...] = ()


@dataclass(frozen=True)
class IndirectLimitResult:
    available: bool
    excluded: bool
    channels_seen: int
    channels_used: int
    max_ratio: float
    energy_gev: float
    flux_cm2_s: float
    limit_cm2_s: float


@dataclass(frozen=True)
class DMSummary:
    point: DMPoint
    result: MicromegasResult
    dir_det_limit: float
    lux_base_limit: float
    indirect_limit: IndirectLimitResult
    relic_excluded: bool
    direct_detection_excluded: bool

    @property
    def passed(self):
        return not (
            self.relic_excluded
            or self.direct_detection_excluded
            or self.indirect_limit.excluded
        )


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


def parse_indirect_line_channels(text):
    channels = []
    for line in text.splitlines():
        if "FermiLAT_line_channel" not in line:
            continue
        energy_match = re.search(rf"E_gamma=({NUMBER_PATTERN})", line)
        flux_match = re.search(rf"Phi_R16=({NUMBER_PATTERN})", line)
        if energy_match is None or flux_match is None:
            continue
        channels.append(
            IndirectLineChannel(
                energy_gev=float(energy_match.group(1)),
                flux_cm2_s=float(flux_match.group(1)),
            )
        )
    return tuple(channels)


def fermi_lat_r16_line_limit(energy_gev):
    limits = FERMI_LAT_R16_LINE_LIMITS
    if (
        not limits
        or energy_gev < limits[0][0]
        or energy_gev > limits[-1][0]
    ):
        return math.inf

    for (low_energy, low_limit), (high_energy, high_limit) in zip(limits, limits[1:]):
        if energy_gev == low_energy:
            return low_limit
        if energy_gev <= high_energy:
            log_energy = (
                (math.log(energy_gev) - math.log(low_energy))
                / (math.log(high_energy) - math.log(low_energy))
            )
            return math.exp(
                math.log(low_limit)
                + log_energy * (math.log(high_limit) - math.log(low_limit))
            )
    return limits[-1][1]


def assess_indirect_limit(channels):
    best = IndirectLimitResult(
        available=False,
        excluded=False,
        channels_seen=0,
        channels_used=0,
        max_ratio=0.0,
        energy_gev=math.nan,
        flux_cm2_s=math.nan,
        limit_cm2_s=math.nan,
    )
    channels_used = 0

    for channel in channels:
        limit = fermi_lat_r16_line_limit(channel.energy_gev)
        if (
            not math.isfinite(limit)
            or limit <= 0.0
            or not math.isfinite(channel.flux_cm2_s)
            or channel.flux_cm2_s < 0.0
        ):
            continue

        channels_used += 1
        ratio = channel.flux_cm2_s / limit
        if not best.available or ratio > best.max_ratio:
            best = IndirectLimitResult(
                available=True,
                excluded=ratio > 1.0,
                channels_seen=len(channels),
                channels_used=channels_used,
                max_ratio=ratio,
                energy_gev=channel.energy_gev,
                flux_cm2_s=channel.flux_cm2_s,
                limit_cm2_s=limit,
            )

    if best.available:
        return IndirectLimitResult(
            available=True,
            excluded=best.max_ratio > 1.0,
            channels_seen=len(channels),
            channels_used=channels_used,
            max_ratio=best.max_ratio,
            energy_gev=best.energy_gev,
            flux_cm2_s=best.flux_cm2_s,
            limit_cm2_s=best.limit_cm2_s,
        )
    return IndirectLimitResult(
        available=False,
        excluded=False,
        channels_seen=len(channels),
        channels_used=0,
        max_ratio=0.0,
        energy_gev=math.nan,
        flux_cm2_s=math.nan,
        limit_cm2_s=math.nan,
    )


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
        indirect_line_channels=parse_indirect_line_channels(text),
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

    indirect_limit = assess_indirect_limit(result.indirect_line_channels)
    return DMSummary(
        point=point,
        result=result,
        dir_det_limit=dir_det_limit,
        lux_base_limit=lux_base_limit,
        indirect_limit=indirect_limit,
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
    if summary.indirect_limit.excluded:
        reasons.append("Fermi-LAT gamma-line flux above limit")
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
        f"  IndirAvailable={summary.indirect_limit.available} "
        f"IndirEnergy={format_value(summary.indirect_limit.energy_gev)} "
        f"IndirFlux={format_value(summary.indirect_limit.flux_cm2_s)} "
        f"IndirLimit={format_value(summary.indirect_limit.limit_cm2_s)} "
        f"IndirRatio={format_value(summary.indirect_limit.max_ratio)}\n"
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
        "dm_indirect_available": summary.indirect_limit.available,
        "dm_indirect_energy": summary.indirect_limit.energy_gev,
        "dm_indirect_flux": summary.indirect_limit.flux_cm2_s,
        "dm_indirect_limit": summary.indirect_limit.limit_cm2_s,
        "dm_indirect_ratio": summary.indirect_limit.max_ratio,
        "dm_indirect_detection_excluded": summary.indirect_limit.excluded,
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
        "dm_indirect_available": None,
        "dm_indirect_energy": None,
        "dm_indirect_flux": None,
        "dm_indirect_limit": None,
        "dm_indirect_ratio": None,
        "dm_indirect_detection_excluded": None,
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
    relic-density, direct-detection, and indirect-detection checks pass. `info`
    is suitable for printing in debug mode, and `exclusion_info` contains the
    numerical values used in the DM exclusion.
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
