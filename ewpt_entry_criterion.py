"""Classify qualitative baryogenesis and gravitational-wave FOPT candidates.

CalcTemps transition indices identify FOPTs independently of MinimaTracer's
global cooling-path steps. These flags do not test that a transition follows
the global path, calculate a sphaleron rate, or predict a GW signal.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping


EW_ENTRY_COLUMNS = (
    "ewpt_ew_entry_true_over_T",
    "ewpt_ew_entry_false_over_T",
    "ewpt_ew_entry_jump_over_T",
    "ewpt_ew_entry_temperature_GeV",
    "ewpt_ew_entry_temperature_kind",
    "ewpt_ew_entry_transition_index",
    "ewpt_ew_entry_nucl_jump_over_T",
    "ewpt_ew_entry_nucl_temperature_GeV",
    "ewpt_ew_entry_perc_jump_over_T",
    "ewpt_ew_entry_perc_temperature_GeV",
    "ewpt_ew_entry_percolated",
    "ewpt_ew_entry_completed",
    "ewpt_baryo_candidate",
    "ewpt_gw_crit_field_jump_over_T",
    "ewpt_gw_crit_temperature_GeV",
    "ewpt_gw_nucl_field_jump_over_T",
    "ewpt_gw_nucl_temperature_GeV",
    "ewpt_gw_perc_field_jump_over_T",
    "ewpt_gw_perc_temperature_GeV",
    "ewpt_gw_max_field_jump_over_T",
    "ewpt_gw_max_temperature_kind",
    "ewpt_gw_max_transition_index",
    "ewpt_gw_candidate",
)

_KINDS = ("crit", "nucl", "perc", "nucl_approx", "compl")
_GW_SELECTION_KINDS = ("crit", "nucl", "perc")


def _finite(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _transition_index(value):
    number = _finite(value)
    if number is None or number < 0 or not number.is_integer():
        return None
    return int(number)


def _field_jump_over_T(strength, false_vev, true_vev, temperature):
    value = _finite(strength.get("field_jump_over_T"))
    if value is not None and value >= 0:
        return value
    coordinates = []
    for name in ("w1", "wx", "ws"):
        before = _finite(false_vev.get(name))
        after = _finite(true_vev.get(name))
        if before is None or after is None:
            return None
        coordinates.append(after - before)
    return math.sqrt(sum(change * change for change in coordinates)) / temperature


def ew_entry_updates(payload, *, w1_threshold=5.0, strength_threshold=1.0):
    """Return independent critical-EW and any-field FOPT candidate flags.

    The baryogenesis flag requires an EW-symmetric-to-broken *critical*
    transition with |delta w1|/Tc > threshold. Nucleation, percolation, and
    completion never gate it. The GW flag requires the Euclidean jump across
    all three scalar VEVs divided by T to exceed threshold for any FOPT at
    critical, nucleation, or percolation temperature. Both are qualitative
    candidate flags rather than physical baryogenesis or GW predictions.
    """
    if not math.isfinite(w1_threshold) or w1_threshold <= 0:
        raise ValueError("w1_threshold must be positive and finite")
    if not math.isfinite(strength_threshold) or strength_threshold <= 0:
        raise ValueError("strength_threshold must be positive and finite")

    updates = dict.fromkeys(EW_ENTRY_COLUMNS)
    # Called only after successful BSMPT runs. No transition is a negative
    # result, while a run that was not attempted retains null scan fields.
    updates["ewpt_baryo_candidate"] = False
    updates["ewpt_gw_candidate"] = False
    groups = defaultdict(dict)
    by_kind = defaultdict(list)
    for strength in payload.get("transition_strengths") or []:
        if not isinstance(strength, Mapping):
            continue
        kind = strength.get("temperature_kind")
        index = _transition_index(strength.get("transition_index"))
        temperature = _finite(strength.get("temperature"))
        false_vev = strength.get("false_vev")
        true_vev = strength.get("true_vev")
        status = strength.get("status")
        if (
            kind not in _KINDS
            or index is None
            or temperature is None
            or temperature <= 0
            or not isinstance(false_vev, Mapping)
            or not isinstance(true_vev, Mapping)
            or (status is not None and str(status).lower() != "success")
        ):
            continue
        false_w1 = _finite(false_vev.get("w1"))
        true_w1 = _finite(true_vev.get("w1"))
        ew_jump = (
            abs(true_w1 - false_w1) / temperature
            if false_w1 is not None and true_w1 is not None
            else None
        )
        field_jump = _field_jump_over_T(strength, false_vev, true_vev, temperature)
        entry = {
            "index": index,
            "kind": kind,
            "temperature": temperature,
            "false_w1": false_w1,
            "true_w1": true_w1,
            "ew_jump": ew_jump,
            "field_jump": field_jump,
        }
        groups[index].setdefault(kind, []).append(entry)
        if kind in _GW_SELECTION_KINDS and field_jump is not None:
            by_kind[kind].append(entry)

    critical_entries = [
        entry
        for kinds in groups.values()
        for entry in kinds.get("crit", [])
        if entry["ew_jump"] is not None
        and abs(entry["false_w1"]) < w1_threshold <= abs(entry["true_w1"])
    ]
    if critical_entries:
        selected = max(
            critical_entries,
            key=lambda item: (item["ew_jump"], item["temperature"], -item["index"]),
        )
        temperature = selected["temperature"]
        index = selected["index"]
        kinds = groups[index]
        updates.update(
            {
                "ewpt_ew_entry_true_over_T": abs(selected["true_w1"]) / temperature,
                "ewpt_ew_entry_false_over_T": abs(selected["false_w1"]) / temperature,
                "ewpt_ew_entry_jump_over_T": selected["ew_jump"],
                "ewpt_ew_entry_temperature_GeV": temperature,
                "ewpt_ew_entry_temperature_kind": "crit",
                "ewpt_ew_entry_transition_index": index,
                "ewpt_ew_entry_percolated": bool(kinds.get("perc") or kinds.get("compl")),
                "ewpt_ew_entry_completed": bool(kinds.get("compl")),
                "ewpt_baryo_candidate": selected["ew_jump"] > strength_threshold,
            }
        )
        for kind in ("nucl", "perc"):
            entries = [entry for entry in kinds.get(kind, []) if entry["ew_jump"] is not None]
            if entries:
                diagnostic = max(entries, key=lambda item: item["ew_jump"])
                updates[f"ewpt_ew_entry_{kind}_jump_over_T"] = diagnostic["ew_jump"]
                updates[f"ewpt_ew_entry_{kind}_temperature_GeV"] = diagnostic["temperature"]

    gw_entries = []
    for kind in _GW_SELECTION_KINDS:
        if not by_kind[kind]:
            continue
        strongest = max(by_kind[kind], key=lambda item: item["field_jump"])
        updates[f"ewpt_gw_{kind}_field_jump_over_T"] = strongest["field_jump"]
        updates[f"ewpt_gw_{kind}_temperature_GeV"] = strongest["temperature"]
        gw_entries.append(strongest)
    if gw_entries:
        strongest = max(
            gw_entries,
            key=lambda item: (item["field_jump"], -_GW_SELECTION_KINDS.index(item["kind"])),
        )
        updates["ewpt_gw_max_field_jump_over_T"] = strongest["field_jump"]
        updates["ewpt_gw_max_temperature_kind"] = strongest["kind"]
        updates["ewpt_gw_max_transition_index"] = strongest["index"]
        updates["ewpt_gw_candidate"] = strongest["field_jump"] > strength_threshold

    return updates
