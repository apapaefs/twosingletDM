"""Sampled global X-breaking history and its relation to DM freeze-out.

MinimaTracer can retain local minima that are not the cosmological vacuum.
Only the reconstructed lowest-potential branch is used here. The result is a
diagnostic on the available temperature grid, not a relic-density cut.
"""

from __future__ import annotations

import bisect
import json
import math
from collections.abc import Mapping


X_HISTORY_COLUMNS = (
    "ewpt_x_broken_min_T_GeV",
    "ewpt_x_broken_max_T_GeV",
    "ewpt_x_broken_intervals_GeV",
    "ewpt_x_final_restoration_low_T_GeV",
    "ewpt_x_final_restoration_high_T_GeV",
    "ewpt_x_phase_at_freezeout",
    "ewpt_x_broken_at_or_after_freezeout",
    "dm_relic_z2_freezeout_compatible",
)


def finite_number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def x_history_updates(payload, freezeout_temperature=None):
    """Summarize sampled global X-breaking windows and freeze-out overlap.

    Intervals store the lowest and highest *sampled* temperatures in each
    contiguous X-broken run, not an interpolated transition temperature.
    A freeze-out point between samples of different phases is unresolved.
    """
    updates = dict.fromkeys(X_HISTORY_COLUMNS)
    minimatracer = payload.get("minimatracer") or {}
    if not isinstance(minimatracer, Mapping):
        return updates
    raw_branch = minimatracer.get("global_branch") or []
    if not isinstance(raw_branch, list):
        return updates

    samples = []
    for point in raw_branch:
        if not isinstance(point, Mapping):
            continue
        temp = finite_number(point.get("temp"))
        label = point.get("label")
        if temp is None or temp < 0 or not isinstance(label, str):
            continue
        samples.append((temp, "X_BROKEN" in label))
    if not samples:
        return updates
    samples.sort()
    if any(a[0] == b[0] for a, b in zip(samples, samples[1:])):
        return updates

    intervals = []
    start = end = None
    for temp, broken in samples:
        if broken:
            if start is None:
                start = temp
            end = temp
        elif start is not None:
            intervals.append([start, end])
            start = end = None
    if start is not None:
        intervals.append([start, end])

    updates["ewpt_x_broken_intervals_GeV"] = json.dumps(intervals, separators=(",", ":"))
    if intervals:
        updates["ewpt_x_broken_min_T_GeV"] = intervals[0][0]
        updates["ewpt_x_broken_max_T_GeV"] = intervals[-1][1]

    # The first unbroken-to-broken boundary in ascending temperature is the
    # final restoration on cooling, provided the trace reaches an unbroken
    # T=0 vacuum. Its two sampled temperatures bracket, rather than locate,
    # the actual transition.
    if samples[0][0] <= 1e-9 and not samples[0][1]:
        for lower_sample, upper_sample in zip(samples, samples[1:]):
            if not lower_sample[1] and upper_sample[1]:
                updates["ewpt_x_final_restoration_low_T_GeV"] = lower_sample[0]
                updates["ewpt_x_final_restoration_high_T_GeV"] = upper_sample[0]
                break

    tf = finite_number(freezeout_temperature)
    if tf is None or tf <= 0:
        return updates
    temperatures = [temp for temp, _ in samples]
    position = bisect.bisect_left(temperatures, tf)
    if position < len(samples) and math.isclose(temperatures[position], tf, rel_tol=0, abs_tol=1e-9):
        phase = "broken" if samples[position][1] else "unbroken"
    elif position == 0 or position == len(samples):
        phase = "outside_traced_range"
    elif samples[position - 1][1] == samples[position][1]:
        phase = "broken" if samples[position][1] else "unbroken"
    else:
        phase = "boundary_unresolved"
    updates["ewpt_x_phase_at_freezeout"] = phase

    # A broken sample at or below Tf is evidence that Z2 is broken during or
    # after nominal freeze-out. To assert the opposite, demand coverage down
    # to T=0 and an unbroken, resolved phase at Tf.
    if any(broken and temp <= tf for temp, broken in samples):
        updates["ewpt_x_broken_at_or_after_freezeout"] = True
    elif temperatures[0] <= 1e-9 and phase == "unbroken":
        updates["ewpt_x_broken_at_or_after_freezeout"] = False
    overlap = updates["ewpt_x_broken_at_or_after_freezeout"]
    if overlap is not None:
        updates["dm_relic_z2_freezeout_compatible"] = not overlap
    return updates
