"""Screen fixed-vacuum relic inputs against the sampled global thermal branch.

These observables describe equilibrium minima and zero-temperature mass gaps.
They do not establish the cosmological phase history or validate a relic density.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Mapping


M1_GEV = 125.09
VEV_SHIFT_SCREEN = 0.10
WINDOW_LOW_FACTOR = 0.5
WINDOW_HIGH_FACTOR = 2.0

RESONANCE_COLUMNS = (
    "dm_resonance_h1_mass_gap_GeV",
    "dm_resonance_h2_mass_gap_GeV",
    "dm_resonance_h1_abs_gap_over_Tf",
    "dm_resonance_h2_abs_gap_over_Tf",
    "dm_resonance_nearest_mediator",
)

THERMAL_VEV_COLUMNS = (
    "dm_relic_thermal_vev_window_low_T_GeV",
    "dm_relic_thermal_vev_window_high_T_GeV",
    "dm_relic_thermal_ew_vev_Tf_over_T0",
    "dm_relic_thermal_s_vev_Tf_over_T0",
    "dm_relic_thermal_ew_vev_max_fractional_shift",
    "dm_relic_thermal_s_vev_max_fractional_shift",
    "dm_relic_thermal_phase_boundary_bracket_overlaps_window",
    "dm_relic_thermal_vev_shift_ge_10pct",
)


def finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def resonance_proximity_updates(m2, m3, freezeout_temperature=None, *, m1=M1_GEV):
    """Record signed zero-T mediator gaps from the two-DM threshold.

    A positive gap puts the pole above threshold. Gap/Tf is a useful screening
    scale, not a resonance significance: widths and thermal masses are absent.
    """
    updates = dict.fromkeys(RESONANCE_COLUMNS)
    m2 = finite_number(m2)
    m3 = finite_number(m3)
    m1 = finite_number(m1)
    tf = finite_number(freezeout_temperature)
    if m1 is None or m2 is None or m3 is None or min(m1, m2, m3) <= 0:
        return updates
    gaps = {"h1": m1 - 2.0 * m3, "h2": m2 - 2.0 * m3}
    updates["dm_resonance_h1_mass_gap_GeV"] = gaps["h1"]
    updates["dm_resonance_h2_mass_gap_GeV"] = gaps["h2"]
    updates["dm_resonance_nearest_mediator"] = min(gaps, key=lambda key: abs(gaps[key]))
    if tf is not None and tf > 0:
        updates["dm_resonance_h1_abs_gap_over_Tf"] = abs(gaps["h1"]) / tf
        updates["dm_resonance_h2_abs_gap_over_Tf"] = abs(gaps["h2"]) / tf
    return updates


def _global_samples(payload):
    minimatracer = payload.get("minimatracer") or {}
    if not isinstance(minimatracer, Mapping):
        return None
    branch = minimatracer.get("global_branch")
    if not isinstance(branch, list) or not branch:
        return None
    samples = []
    for point in branch:
        if not isinstance(point, Mapping):
            return None
        temperature = finite_number(point.get("temp"))
        ew_vev = finite_number(point.get("w1"))
        s_vev = finite_number(point.get("ws"))
        if temperature is None or temperature < 0 or ew_vev is None or s_vev is None:
            return None
        phase = point.get("phase_index")
        if phase is None:
            phase = point.get("label")
        samples.append((temperature, abs(ew_vev), abs(s_vev), phase))
    samples.sort(key=lambda sample: sample[0])
    if any(left[0] == right[0] for left, right in zip(samples, samples[1:])):
        return None
    return samples


def _vevs_at(samples, temperature):
    temperatures = [sample[0] for sample in samples]
    position = bisect.bisect_left(temperatures, temperature)
    if position < len(samples) and math.isclose(
        temperatures[position], temperature, rel_tol=0.0, abs_tol=1e-9
    ):
        return samples[position][1:3]
    if position == 0 or position == len(samples):
        return None
    lower, upper = samples[position - 1], samples[position]
    if lower[3] is None or upper[3] is None or lower[3] != upper[3]:
        return None  # Never interpolate across different minima.
    fraction = (temperature - lower[0]) / (upper[0] - lower[0])
    return tuple(lower[index] + fraction * (upper[index] - lower[index]) for index in (1, 2))


def thermal_vev_updates(payload, freezeout_temperature):
    """Compare VEV magnitudes around Tf with the traced T=0 vacuum.

    The window [Tf/2, 2 Tf] is a screening convention. Maxima use samples in
    that window and same-phase interpolation at its endpoints. A boundary
    bracket overlapping the window is reported separately; its location is
    unresolved between samples. Missing coverage never produces a false flag.
    """
    updates = dict.fromkeys(THERMAL_VEV_COLUMNS)
    tf = finite_number(freezeout_temperature)
    if tf is None or tf <= 0:
        return updates
    low, high = WINDOW_LOW_FACTOR * tf, WINDOW_HIGH_FACTOR * tf
    updates["dm_relic_thermal_vev_window_low_T_GeV"] = low
    updates["dm_relic_thermal_vev_window_high_T_GeV"] = high
    samples = _global_samples(payload)
    if not samples or samples[0][0] > 1e-9:
        return updates
    zero = samples[0][1:3]
    if zero[0] <= 1e-9 or zero[1] <= 1e-9:
        return updates

    at_tf = _vevs_at(samples, tf)
    if at_tf is not None:
        updates["dm_relic_thermal_ew_vev_Tf_over_T0"] = at_tf[0] / zero[0]
        updates["dm_relic_thermal_s_vev_Tf_over_T0"] = at_tf[1] / zero[1]

    boundary = any(
        left[3] != right[3] and left[0] <= high and right[0] >= low
        for left, right in zip(samples, samples[1:])
        if left[3] is not None and right[3] is not None
    )
    full_coverage = low >= samples[0][0] and high <= samples[-1][0]
    if full_coverage and all(sample[3] is not None for sample in samples):
        updates["dm_relic_thermal_phase_boundary_bracket_overlaps_window"] = boundary

    endpoint_vevs = [_vevs_at(samples, low), _vevs_at(samples, high)]
    vevs = [sample[1:3] for sample in samples if low <= sample[0] <= high]
    vevs.extend(value for value in endpoint_vevs if value is not None)
    if vevs:
        ew_shift = max(abs(value[0] / zero[0] - 1.0) for value in vevs)
        s_shift = max(abs(value[1] / zero[1] - 1.0) for value in vevs)
        updates["dm_relic_thermal_ew_vev_max_fractional_shift"] = ew_shift
        updates["dm_relic_thermal_s_vev_max_fractional_shift"] = s_shift
        if max(ew_shift, s_shift) >= VEV_SHIFT_SCREEN:
            updates["dm_relic_thermal_vev_shift_ge_10pct"] = True
        elif (
            updates["dm_relic_thermal_phase_boundary_bracket_overlaps_window"] is False
            and all(value is not None for value in endpoint_vevs)
        ):
            updates["dm_relic_thermal_vev_shift_ge_10pct"] = False
    return updates
