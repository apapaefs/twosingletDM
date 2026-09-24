#!/usr/bin/env python3

"""Re-evaluate HiggsTools and dark-matter results in an existing vx=0 scan.

The input is never modified.  Progress is stored transactionally in a SQLite
checkpoint named ``<output>.partial``; the requested TSV appears only after
all rows have completed and passed result validation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sqlite3
from pathlib import Path
from typing import Callable, Mapping, Sequence

from dm_thermal_relic_diagnostic import (
    RESONANCE_COLUMNS,
    THERMAL_VEV_COLUMNS,
    resonance_proximity_updates,
)
from trsm_cmb import CMB_COLUMNS, add_cmb_arguments, cmb_configuration, cmb_diagnostics, require_cmb_capability
from trsm_direct_detection import DEFAULT_LIMIT_MODEL, DEFAULT_LIMIT_TABLE, load_si_limit_table
from trsm_micromegas import require_v2_capability, DEFAULT_MICROMEGAS_VERSION, MICROMEGAS_VERSIONS, micromegas_configuration, normalize_micromegas_version
from trsm_scan_campaign import atomic_write_json


EXPECTED_CONVENTION_ID = "trsm_vxzero_canonical_v1"
from trsm_inputs import M1 as M1_GEV, PHYSICS_VERSION
from trsm_constraint_profile import profile_updates, physics_manifest, precision_updates, validate_v2_result
from trsm_theory_diagnostics import theory_diagnostics
from scan_output import V2_COLUMNS
from trsm_flavour import FLAVOUR_COLUMNS, generated_flavour_updates
from ewpt_entry_criterion import EW_ENTRY_COLUMNS
CHECKPOINT_SCHEMA_VERSION = 1

REQUIRED_INPUT_COLUMNS = (
    "M2",
    "M3",
    "vs",
    "vx",
    "a12",
    "lX",
    "lPhiX",
    "lSX",
)

DM_COLUMNS = (
    "dm_mdm",
    "dm_omega",
    "dm_relic_upper_limit",
    "dm_dir_det",
    "dm_dir_det_limit",
    "dm_lux_base_limit",
    "dm_relic_excluded",
    "dm_direct_detection_excluded",
    "dm_indirect_available",
    "dm_indirect_channels_seen",
    "dm_indirect_channels_used",
    "dm_indirect_energy",
    "dm_indirect_flux",
    "dm_indirect_limit",
    "dm_indirect_ratio",
    "dm_indirect_detection_excluded",
    "dm_limit_model",
    "dm_rescale",
    "dm_xf",
    "dm_freezeout_temperature_GeV",
)

HIGGSTOOLS_COLUMNS = (
    "higgstools_hb_selected_limits",
    "higgstools_hb_top_obs",
    "higgstools_hs_chi2",
    "higgstools_hs_delta_chi2",
    "higgstools_hs_top_chi2",
)

PROVENANCE_COLUMNS = (
    "h1_h3h3_width",
    "h1_h3h3_br",
    "h1_h2h2_width",
    "h1_h2h2_br",
    "h2_h3h3_width",
    "h2_h3h3_br",
    "higgs_invisible_widths_included",
    "portal_convention",
    "micromegas_model_convention",
)

UPDATED_COLUMNS = (
    "K133",
    "K233",
    "w1",
    "w2",
    "w3",
    "hb",
    "hs",
    "dm",
) + DM_COLUMNS + HIGGSTOOLS_COLUMNS + PROVENANCE_COLUMNS

FREEZEOUT_DEPENDENT_COLUMNS = (
    "ewpt_x_phase_at_freezeout",
    "ewpt_x_broken_at_or_after_freezeout",
    "dm_relic_z2_freezeout_compatible",
    *THERMAL_VEV_COLUMNS,
)

BOOLEAN_RESULT_COLUMNS = (
    "hb",
    "hs",
    "dm",
    "dm_relic_excluded",
    "dm_direct_detection_excluded",
    "dm_indirect_available",
    "dm_indirect_detection_excluded",
    "dm_rescale",
    "higgs_invisible_widths_included",
)

FINITE_RESULT_COLUMNS = (
    "K133",
    "K233",
    "w1",
    "w2",
    "w3",
    "h1_h3h3_width",
    "h1_h3h3_br",
    "h1_h2h2_width",
    "h1_h2h2_br",
    "h2_h3h3_width",
    "h2_h3h3_br",
    "dm_mdm",
    "dm_omega",
    "dm_relic_upper_limit",
    "dm_dir_det",
    "dm_dir_det_limit",
    "dm_lux_base_limit",
    "dm_indirect_ratio",
    "higgstools_hs_chi2",
    "higgstools_hs_delta_chi2",
)


class ReEvaluationError(RuntimeError):
    """Raised when a row cannot be safely re-evaluated."""


def compact_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def require_float(row: Mapping[str, str], column: str, row_number: int) -> float:
    try:
        value = float(row[column])
    except (KeyError, TypeError, ValueError) as exc:
        raise ReEvaluationError(
            f"invalid {column!r} value on input row {row_number}"
        ) from exc
    if not math.isfinite(value):
        raise ReEvaluationError(
            f"non-finite {column!r} value on input row {row_number}"
        )
    return value


def validate_updates(updates: Mapping[str, object], row_number: int) -> None:
    if updates.get("constraint_version") == PHYSICS_VERSION:
        validate_v2_result(updates)
        return
    missing = [name for name in UPDATED_COLUMNS if name not in updates]
    if missing:
        raise ReEvaluationError(
            f"row {row_number} result is missing columns: {', '.join(missing)}"
        )

    cmb_pass = True
    if any(name in updates for name in CMB_COLUMNS):
        if any(name not in updates for name in CMB_COLUMNS):
            raise ReEvaluationError(f"row {row_number} has incomplete CMB diagnostics")
        enabled, available = updates["dm_cmb_enabled"], updates["dm_cmb_available"]
        if type(enabled) is not bool or type(available) is not bool:
            raise ReEvaluationError(f"row {row_number} has invalid CMB status flags")
        if available:
            if not enabled or updates["dm_cmb_status"] != "ok":
                raise ReEvaluationError(f"row {row_number} has inconsistent CMB status")
            for name in ("dm_cmb_ratio_raw", "dm_cmb_abundance_fraction", "dm_cmb_ratio"):
                value = updates[name]
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                    raise ReEvaluationError(f"row {row_number} has invalid {name}")
            fraction = updates["dm_cmb_abundance_fraction"]
            ratio = updates["dm_cmb_ratio"]
            if fraction > 1 or not math.isclose(ratio, updates["dm_cmb_ratio_raw"] * fraction**2, rel_tol=1e-12):
                raise ReEvaluationError(f"row {row_number} has inconsistent CMB rescaling")
            if type(updates["dm_cmb_excluded"]) is not bool or updates["dm_cmb_excluded"] != (ratio > 1):
                raise ReEvaluationError(f"row {row_number} has inconsistent CMB exclusion")
        else:
            if any(updates[name] is not None for name in ("dm_cmb_excluded", "dm_cmb_ratio_raw", "dm_cmb_ratio")):
                raise ReEvaluationError(f"row {row_number} retains unavailable CMB results")
            if not enabled and any(updates[name] != value for name, value in cmb_diagnostics().items()):
                raise ReEvaluationError(f"row {row_number} retains disabled CMB results")
            if enabled and (updates["dm_cmb_status"] in ("ok", "disabled", "") or not updates["dm_cmb_reason"]):
                raise ReEvaluationError(f"row {row_number} has no CMB failure reason")
        cmb_pass = not enabled or (available and not updates["dm_cmb_excluded"])

    # micrOMEGAs can return a non-finite relic density for an otherwise valid
    # scan point.  The original scan records that as DM-failed with all DM
    # details unavailable; keep the same semantics while updating HiggsTools.
    dm_unavailable = all(updates[column] is None for column in DM_COLUMNS)
    if updates.get("dm_cmb_available"):
        if dm_unavailable:
            raise ReEvaluationError(f"row {row_number} has CMB results without core DM data")
        expected_fraction = min(1.0, float(updates["dm_omega"]) / 0.12) if updates["dm_rescale"] else 1.0
        if not math.isclose(updates["dm_cmb_abundance_fraction"], expected_fraction, rel_tol=1e-12):
            raise ReEvaluationError(f"row {row_number} has inconsistent CMB abundance fraction")
    if dm_unavailable and updates["dm"] is not False:
        raise ReEvaluationError(
            f"row {row_number} has unavailable DM details but is marked DM-passing"
        )

    for column in BOOLEAN_RESULT_COLUMNS:
        if dm_unavailable and column in DM_COLUMNS:
            continue
        if type(updates[column]) is not bool:
            raise ReEvaluationError(
                f"row {row_number} result {column!r} is not boolean"
            )

    for column in FINITE_RESULT_COLUMNS:
        if dm_unavailable and column in DM_COLUMNS:
            continue
        if column == "dm_dir_det_limit" and updates["dm_rescale"] and updates["dm_omega"] == 0 and updates[column] == math.inf:
            continue  # At zero abundance, the relic-rescaled cross-section bound is unbounded.
        try:
            finite = math.isfinite(float(updates[column]))
        except (TypeError, ValueError):
            finite = False
        if not finite:
            raise ReEvaluationError(
                f"row {row_number} result {column!r} is not finite"
            )

    for column in ("w1", "w2"):
        if float(updates[column]) <= 0.0:
            raise ReEvaluationError(f"row {row_number} result {column!r} is not positive")
    if float(updates["w3"]) < 0.0:
        raise ReEvaluationError(f"row {row_number} result 'w3' is negative")
    for column in ("h1_h3h3_width", "h1_h2h2_width", "h2_h3h3_width"):
        if float(updates[column]) < 0.0:
            raise ReEvaluationError(f"row {row_number} result {column!r} is negative")
    for column in ("h1_h3h3_br", "h1_h2h2_br", "h2_h3h3_br"):
        if not 0.0 <= float(updates[column]) <= 1.0:
            raise ReEvaluationError(
                f"row {row_number} result {column!r} is outside [0, 1]"
            )
    for width_column, br_column, total_column in (
        ("h1_h3h3_width", "h1_h3h3_br", "w1"),
        ("h1_h2h2_width", "h1_h2h2_br", "w1"),
        ("h2_h3h3_width", "h2_h3h3_br", "w2"),
    ):
        expected_br = float(updates[width_column]) / float(updates[total_column])
        if not math.isclose(
            float(updates[br_column]), expected_br, rel_tol=1.0e-10, abs_tol=1.0e-12
        ):
            raise ReEvaluationError(
                f"row {row_number} result {br_column!r} is inconsistent with its widths"
            )
    if not dm_unavailable:
        if float(updates["dm_mdm"]) <= 0.0:
            raise ReEvaluationError(f"row {row_number} dark-matter mass is not positive")
        xf = updates["dm_xf"]
        freezeout_temperature = updates["dm_freezeout_temperature_GeV"]
        if (xf is None) != (freezeout_temperature is None):
            raise ReEvaluationError(f"row {row_number} has incomplete freeze-out diagnostics")
        if xf is not None:
            try:
                valid_freezeout = (
                    math.isfinite(float(xf)) and float(xf) > 0.0
                    and math.isfinite(float(freezeout_temperature))
                    and float(freezeout_temperature) > 0.0
                    and math.isclose(
                        float(freezeout_temperature), float(updates["dm_mdm"]) / float(xf),
                        rel_tol=1.0e-10,
                    )
                )
            except (TypeError, ValueError, ZeroDivisionError):
                valid_freezeout = False
            if not valid_freezeout:
                raise ReEvaluationError(f"row {row_number} has inconsistent freeze-out diagnostics")
        if float(updates["dm_omega"]) < 0.0 or float(updates["dm_dir_det"]) < 0.0:
            raise ReEvaluationError(f"row {row_number} has a negative core DM result")
        if float(updates["dm_relic_upper_limit"]) <= 0.0:
            raise ReEvaluationError(f"row {row_number} relic upper limit is not positive")
        if float(updates["dm_dir_det_limit"]) <= 0.0:
            raise ReEvaluationError(
                f"row {row_number} direct-detection limit is not positive"
            )

        for column in ("dm_indirect_channels_seen", "dm_indirect_channels_used"):
            value = updates[column]
            try:
                valid_integer = (
                    not isinstance(value, bool)
                    and int(value) == value
                    and int(value) >= 0
                )
            except (TypeError, ValueError, OverflowError):
                valid_integer = False
            if not valid_integer:
                raise ReEvaluationError(
                    f"row {row_number} result {column!r} is not a non-negative integer"
                )
        if int(updates["dm_indirect_channels_used"]) > int(
            updates["dm_indirect_channels_seen"]
        ):
            raise ReEvaluationError(
                f"row {row_number} uses more indirect channels than were seen"
            )

        if updates["dm_indirect_available"]:
            for column in ("dm_indirect_energy", "dm_indirect_flux", "dm_indirect_limit"):
                try:
                    finite = math.isfinite(float(updates[column]))
                except (TypeError, ValueError):
                    finite = False
                if not finite:
                    raise ReEvaluationError(
                        f"row {row_number} available indirect result {column!r} is not finite"
                    )
        elif updates["dm_indirect_detection_excluded"]:
            raise ReEvaluationError(
                f"row {row_number} marks unavailable indirect data as excluded"
            )

        component_pass = not (
            updates["dm_relic_excluded"]
            or updates["dm_direct_detection_excluded"]
            or updates["dm_indirect_detection_excluded"]
        )
        if updates["dm"] != (component_pass and cmb_pass):
            raise ReEvaluationError(
                f"row {row_number} aggregate DM result disagrees with its components"
            )
        if not isinstance(updates["dm_limit_model"], str) or not updates["dm_limit_model"]:
            raise ReEvaluationError(f"row {row_number} has no DM limit-model provenance")

    for column in ("portal_convention", "micromegas_model_convention"):
        if updates[column] != EXPECTED_CONVENTION_ID:
            raise ReEvaluationError(
                f"row {row_number} has unexpected {column}: {updates[column]!r}"
            )


class CoreEvaluator:
    """Adapter around the production TRSM, HiggsTools, and microOMEGAs APIs."""

    def __init__(self, micromegas_main: Path | None = None, *, planck_cmb=False, limit_table=None, limit_model=DEFAULT_LIMIT_MODEL):
        import generate_trsm_info as info
        import test_trsm_DM as dm_provider
        import test_trsm_higgstools as higgs_provider

        convention = getattr(info, "PORTAL_CONVENTION_ID", None)
        if convention != EXPECTED_CONVENTION_ID:
            raise ReEvaluationError(
                "generate_trsm_info does not expose the expected canonical convention"
            )
        for name in (
            "generate_lams",
            "vxzero_portal_couplings",
            "scalar_to_identical_scalar_width",
            "vxzero_invisible_decay_info",
        ):
            if not hasattr(info, name):
                raise ReEvaluationError(f"generate_trsm_info is missing {name}()")

        if micromegas_main is None:
            micromegas_main = dm_provider.default_micromegas_main()
        self.micromegas_main = Path(micromegas_main).expanduser().resolve()
        if not self.micromegas_main.is_file():
            raise FileNotFoundError(
                f"micrOMEGAs executable not found: {self.micromegas_main}"
            )

        self.info = info
        self.dm_provider = dm_provider
        self.higgs = higgs_provider
        self.planck_cmb = planck_cmb
        self.limit_table = limit_table
        self.limit_model = limit_model
        if planck_cmb:
            require_cmb_capability(self.micromegas_main)

    def __call__(self, row: Mapping[str, str], point_index: int) -> dict[str, object]:
        row_number = point_index + 1
        m2 = require_float(row, "M2", row_number)
        m3 = require_float(row, "M3", row_number)
        vs_input = require_float(row, "vs", row_number)
        vx = require_float(row, "vx", row_number)
        a12_input = require_float(row, "a12", row_number)
        lx = require_float(row, "lX", row_number)
        lphix = require_float(row, "lPhiX", row_number)
        lsx = require_float(row, "lSX", row_number)
        if not math.isclose(vx, 0.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise ReEvaluationError(
                f"input row {row_number} has vx={vx}; only vx=0 scans are supported"
            )

        generated = self.info.generate_lams(
            point_index,
            m2,
            m3,
            vs_input,
            0.0,
            a12_input,
            0.0,
            0.0,
            False,
            lX=lx,
            lPhiX=lphix,
            lSX=lsx,
        )
        if len(generated) != 28:
            raise ReEvaluationError(
                f"generate_lams returned {len(generated)} fields; expected 28"
            )
        (
            vs,
            _vx,
            mass2,
            mass3,
            a12,
            _a13,
            _a23,
            w1,
            w2,
            w3,
            _k111,
            _k112,
            _k113,
            _k123,
            k122,
            _k1111,
            _k1112,
            _k1113,
            generated_k133,
            k1,
            k2,
            k3,
            h1_brs,
            h2_brs,
            h3_brs,
            _xs1,
            _xs2,
            _xs3,
        ) = generated

        k133, k233 = self.info.vxzero_portal_couplings(lphix, lsx, vs, a12)
        if not math.isclose(
            float(generated_k133), float(k133), rel_tol=1.0e-10, abs_tol=1.0e-10
        ):
            raise ReEvaluationError(
                f"input row {row_number} has inconsistent generated K133"
            )
        invisible = self.info.vxzero_invisible_decay_info(
            M1_GEV,
            mass2,
            mass3,
            k133,
            k233,
            h1_brs[-1],
            h2_brs[-1],
            K122=k122,
        )
        gamma1 = invisible["h1_h3h3_width"]
        gamma2 = invisible["h2_h3h3_width"]
        gamma1_h2h2 = invisible["h1_h2h2_width"]
        expected_w1 = invisible["w1"]
        expected_w2 = invisible["w2"]
        if not math.isclose(float(w1), expected_w1, rel_tol=1.0e-10, abs_tol=1.0e-12):
            raise ReEvaluationError(
                f"input row {row_number} has inconsistent physical h1 width"
            )
        if not math.isclose(float(w2), expected_w2, rel_tol=1.0e-10, abs_tol=1.0e-12):
            raise ReEvaluationError(
                f"input row {row_number} has inconsistent physical h2 width"
            )

        higgs_error = None
        try:
            higgs_result = self.higgs.analyze_parampoint(
                self.higgs.pred,
                self.higgs.H1,
                self.higgs.H2,
                self.higgs.H3,
                M1_GEV,
                mass2,
                mass3,
                k1,
                k2,
                k3,
                h1_brs,
                h2_brs,
                h3_brs,
                h1_direct_invisible_width=gamma1,
                h2_direct_invisible_width=gamma2,
                h1_h2h2_width=gamma1_h2h2,
                return_details=True,
            )
        except (ValueError, RuntimeError, ArithmeticError) as error:
            higgs_error = str(error)
            higgs_result = (None, None, {"higgsbounds": {}, "higgssignals": {}})
        if len(higgs_result) < 3:
            raise ReEvaluationError(
                f"input row {row_number} produced no HiggsTools details"
            )
        hb, hs, details = higgs_result[:3]
        if (hb is not None and type(hb) is not bool) or (hs is not None and type(hs) is not bool) or not isinstance(details, dict):
            raise ReEvaluationError(
                f"input row {row_number} produced unusable HB/HS results"
            )
        hb_details = details.get("higgsbounds")
        hs_details = details.get("higgssignals")
        if not isinstance(hb_details, dict) or not isinstance(hs_details, dict):
            raise ReEvaluationError(
                f"input row {row_number} produced malformed HiggsTools details"
            )

        dm_passed, dm_info, dm_values = self.dm_provider.test_dm(
            lx,
            lphix,
            lsx,
            mass3,
            vs,
            a12,
            mass2,
            point_index=point_index,
            micromegas_main=self.micromegas_main,
            planck_cmb=self.planck_cmb,
            limit_table=self.limit_table,
            limit_model=self.limit_model,
        )
        if dm_passed is not None and type(dm_passed) is not bool or not isinstance(dm_values, dict):
            raise ReEvaluationError(
                f"input row {row_number} produced unusable DM results: {dm_info}"
            )
        if dm_values and all(value is None for value in dm_values.values()):
            print(f"Input row {row_number}: {dm_info}", flush=True)

        updates: dict[str, object] = {
            "K133": k133,
            "K233": k233,
            "w1": w1,
            "w2": w2,
            "w3": w3,
            "hb": hb,
            "hs": hs,
            "dm": dm_passed,
            "higgstools_hb_selected_limits": compact_json(
                hb_details.get("selected_limits", {})
            ),
            "higgstools_hb_top_obs": compact_json(
                hb_details.get("top_observed_ratios", [])
            ),
            "higgstools_hs_chi2": hs_details.get("chi2"),
            "higgstools_hs_delta_chi2": hs_details.get("delta_chi2"),
            "higgstools_hs_top_chi2": compact_json(
                hs_details.get("top_chi2_contributors", [])
            ),
            "h1_h3h3_width": gamma1,
            "h1_h3h3_br": invisible["h1_h3h3_br"],
            "h1_h2h2_width": gamma1_h2h2,
            "h1_h2h2_br": invisible["h1_h2h2_br"],
            "h2_h3h3_width": gamma2,
            "h2_h3h3_br": invisible["h2_h3h3_br"],
            "higgs_invisible_widths_included": True,
            "portal_convention": EXPECTED_CONVENTION_ID,
            "micromegas_model_convention": EXPECTED_CONVENTION_ID,
        }
        updates.update(dm_values)
        updates.update(generated_flavour_updates(mass2, mass3, (k1, k2, k3),
            (h1_brs, h2_brs, h3_brs), (w1, w2, w3), vx=0))
        updates.update(resonance_proximity_updates(
            mass2, mass3, updates.get("dm_freezeout_temperature_GeV")
        ))
        if not self.planck_cmb:
            updates.update(cmb_diagnostics())
        from test_trsm_theory_constraints import theory_constraints_vxzero
        updates['thc']=bool(theory_constraints_vxzero(vs,mass2,mass3,a12,lx,lphix,lsx))
        updates.update(theory_diagnostics(vs,mass2,mass3,a12,lx,lphix,lsx))
        updates.update(precision_updates(mass2,a12))
        updates['evo']=updates['rg_integration_success']  # legacy numeric field
        point={**dict(row),**updates,'M2':mass2,'M3':mass3}
        updates.update(profile_updates(point))
        if higgs_error:
            updates["point_assessment_reason"] += "; HiggsTools: " + higgs_error
        updates['point_index']=point_index
        updates.update(dict.fromkeys(EW_ENTRY_COLUMNS))
        updates['ewpt_constraint_version']=None
        updates['ewpt_status']='requires_v2_reevaluation'
        validate_updates(updates, row_number)
        return updates


def read_input(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ReEvaluationError(f"input is empty: {path}") from exc
        if not header or any(not name for name in header):
            raise ReEvaluationError("input header contains an empty column name")
        duplicates = sorted({name for name in header if header.count(name) > 1})
        if duplicates:
            raise ReEvaluationError(
                f"input has duplicate columns: {', '.join(duplicates)}"
            )
        missing = [name for name in REQUIRED_INPUT_COLUMNS if name not in header]
        if missing:
            raise ReEvaluationError(
                f"input is missing columns: {', '.join(missing)}"
            )

        rows = []
        for row_number, values in enumerate(reader, start=2):
            if len(values) != len(header):
                raise ReEvaluationError(
                    f"input row {row_number} has {len(values)} fields; expected {len(header)}"
                )
            rows.append(dict(zip(header, values)))
    if not rows:
        raise ReEvaluationError(f"input has a header but no data rows: {path}")
    return header, rows


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def input_identity(path: Path, header: Sequence[str], row_count: int) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "size": path.stat().st_size,
        "header": list(header),
        "row_count": row_count,
    }


def output_header(input_header: Sequence[str], *, planck_cmb=False) -> list[str]:
    result = list(input_header)
    result.extend(name for name in UPDATED_COLUMNS if name not in result)
    result.extend(name for name in FLAVOUR_COLUMNS if name not in result)
    result.extend(name for name in RESONANCE_COLUMNS if name not in result)
    result.extend(name for name in (*V2_COLUMNS,*EW_ENTRY_COLUMNS,"thc","evo","ewpo","wmass","ewpt_status","ewpt_constraint_version","legacy_ewpt_baryo_candidate","legacy_ewpt_gw_candidate") if name not in result)
    if planck_cmb or any(name in result for name in CMB_COLUMNS):
        result.extend(name for name in CMB_COLUMNS if name not in result)
    return result


def format_value(value: object) -> str:
    if value is None:
        return "nan"
    if type(value) is bool:
        return "True" if value else "False"
    return str(value)


def merged_payload(
    input_row: Mapping[str, str], updates: Mapping[str, object], header: Sequence[str]
) -> str:
    fields = []
    for column in header:
        if column in updates:
            fields.append(format_value(updates[column]))
        else:
            fields.append(input_row.get(column,"nan"))
    return "\t".join(fields)


def checkpoint_path(output: Path) -> Path:
    return output.with_name(output.name + ".partial")


def checkpoint_metadata(
    identity: Mapping[str, object], header: Sequence[str], output: Path,
    evaluation_configuration=None,
) -> dict[str, object]:
    metadata = {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "convention": EXPECTED_CONVENTION_ID,
        "input_identity": dict(identity),
        "output_header": list(header),
        "output_path": str(output.resolve()),
    }
    if evaluation_configuration is not None:
        metadata["schema_version"] = 2
        metadata["evaluation_configuration"] = evaluation_configuration
    return metadata


def read_checkpoint_metadata(path):
    # Read-only: a misspelled checkpoint must never create an empty database.
    with sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True) as connection:
        return {key: json.loads(value) for key, value in connection.execute("SELECT key, value FROM metadata")}


def create_checkpoint(path: Path, metadata: Mapping[str, object]) -> sqlite3.Connection:
    connection = sqlite3.connect(path)
    try:
        with connection:
            connection.execute(
                "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
            )
            connection.execute(
                "CREATE TABLE rows (row_index INTEGER PRIMARY KEY, payload TEXT NOT NULL)"
            )
            connection.executemany(
                "INSERT INTO metadata(key, value) VALUES (?, ?)",
                [
                    (key, json.dumps(value, sort_keys=True, separators=(",", ":")))
                    for key, value in metadata.items()
                ],
            )
    except Exception:
        connection.close()
        raise
    return connection


def open_checkpoint(
    path: Path, expected_metadata: Mapping[str, object]
) -> tuple[sqlite3.Connection, int]:
    try:
        connection = sqlite3.connect(path)
        stored = {
            key: json.loads(value)
            for key, value in connection.execute("SELECT key, value FROM metadata")
        }
        if stored != dict(expected_metadata):
            raise ReEvaluationError(
                "checkpoint does not match the current input, output, convention, or evaluation configuration"
            )
        count, minimum, maximum = connection.execute(
            "SELECT COUNT(*), MIN(row_index), MAX(row_index) FROM rows"
        ).fetchone()
        if count and (minimum != 0 or maximum != count - 1):
            raise ReEvaluationError("checkpoint rows are not contiguous")
        return connection, int(count)
    except (sqlite3.DatabaseError, json.JSONDecodeError, ReEvaluationError) as exc:
        try:
            connection.close()
        except UnboundLocalError:
            pass
        if isinstance(exc, ReEvaluationError):
            raise
        raise ReEvaluationError(f"invalid checkpoint {path}: {exc}") from exc


def export_completed_output(
    connection: sqlite3.Connection,
    output: Path,
    header: Sequence[str],
    expected_rows: int,
) -> None:
    temporary = output.with_name("." + output.name + ".complete.tmp")
    count = connection.execute("SELECT COUNT(*) FROM rows").fetchone()[0]
    if count != expected_rows:
        raise ReEvaluationError(
            f"checkpoint contains {count} rows; expected {expected_rows}"
        )
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        stream.write("\t".join(header) + "\n")
        for expected_index, (row_index, payload) in enumerate(
            connection.execute("SELECT row_index, payload FROM rows ORDER BY row_index")
        ):
            if row_index != expected_index:
                raise ReEvaluationError("checkpoint rows are not contiguous")
            stream.write(payload + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output)


def reevaluate(
    input_path: Path,
    output: Path,
    evaluator: Callable[[Mapping[str, str], int], Mapping[str, object]],
    *,
    resume: bool = False,
    checkpoint_every: int = 25,
    evaluation_configuration=None,
) -> Path:
    input_path = input_path.expanduser().resolve()
    output = output.expanduser().resolve()
    if input_path == output:
        raise ReEvaluationError("input and output paths must differ")
    if not input_path.is_file():
        raise FileNotFoundError(f"input scan not found: {input_path}")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if checkpoint_every <= 0:
        raise ValueError("checkpoint_every must be positive")

    header, rows = read_input(input_path)
    identity = input_identity(input_path, header, len(rows))
    planck_cmb = bool((evaluation_configuration or {}).get("planck_cmb", {}).get("enabled", False))
    final_header = output_header(header, planck_cmb=planck_cmb)
    metadata = checkpoint_metadata(identity, final_header, output, evaluation_configuration)
    provenance_path = output.with_suffix(".metadata.json")
    if evaluation_configuration is not None and provenance_path.exists():
        raise FileExistsError(f"refusing to overwrite existing metadata: {provenance_path}")
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = checkpoint_path(output)

    if resume:
        if not partial.is_file():
            raise FileNotFoundError(f"checkpoint not found: {partial}")
        stored = read_checkpoint_metadata(partial)
        if evaluation_configuration is not None and "evaluation_configuration" not in stored:
            raise ReEvaluationError("Legacy checkpoint has no physics configuration fingerprint; restart reevaluation")
        connection, completed = open_checkpoint(partial, metadata)
    else:
        if partial.exists():
            raise FileExistsError(
                f"checkpoint already exists; use --resume or move it aside: {partial}"
            )
        connection = create_checkpoint(partial, metadata)
        completed = 0

    if completed > len(rows):
        connection.close()
        raise ReEvaluationError("checkpoint has more rows than the input")
    if completed:
        print(f"Resuming at row {completed + 1:,}/{len(rows):,}")

    try:
        for batch_start in range(completed, len(rows), checkpoint_every):
            batch_stop = min(batch_start + checkpoint_every, len(rows))
            current_index = batch_start
            try:
                connection.execute("BEGIN")
                for index in range(batch_start, batch_stop):
                    current_index = index
                    updates = dict(evaluator(rows[index], index + 1))
                    for flag in ("ewpt_baryo_candidate","ewpt_gw_candidate"):
                        updates["legacy_"+flag]=rows[index].get("legacy_"+flag,rows[index].get(flag))
                    updates.update(resonance_proximity_updates(
                        rows[index].get("M2"), rows[index].get("M3"),
                        updates.get("dm_freezeout_temperature_GeV"),
                        widths=(updates.get("dm_h1_width_GeV"),updates.get("dm_h2_width_GeV")),
                    ))
                    # A changed micrOMEGAs Xf invalidates comparisons made with
                    # an earlier BSMPT result. Retain the X history itself.
                    for column in FREEZEOUT_DEPENDENT_COLUMNS:
                        if column in rows[index]:
                            updates[column] = None
                    if any(name in final_header for name in CMB_COLUMNS):
                        if planck_cmb and (any(name not in updates for name in CMB_COLUMNS) or updates["dm_cmb_enabled"] is not True):
                            raise ReEvaluationError("CMB-enabled evaluator returned no complete, enabled CMB diagnostics")
                        if not planck_cmb:
                            updates.update(cmb_diagnostics())
                    validate_updates(updates, index + 2)
                    payload = merged_payload(rows[index], updates, final_header)
                    connection.execute(
                        "INSERT INTO rows(row_index, payload) VALUES (?, ?)",
                        (index, payload),
                    )
                connection.commit()
            except Exception as exc:
                connection.rollback()
                if isinstance(exc, ReEvaluationError):
                    raise
                raise ReEvaluationError(
                    f"re-evaluation failed on input row {current_index + 2}: {exc}"
                ) from exc
            print(f"Checkpointed {batch_stop:,}/{len(rows):,} rows")

        if input_identity(input_path, header, len(rows)) != identity:
            raise ReEvaluationError("input changed while re-evaluation was running")
        export_completed_output(connection, output, final_header, len(rows))
        if evaluation_configuration is not None:
            atomic_write_json(provenance_path, {
                "schema": "trsm_reevaluation_metadata_v1",
                "input_identity": identity, "scan_file": output.name,
                "evaluation_configuration": evaluation_configuration,
                **evaluation_configuration,
            })
    finally:
        connection.close()

    partial.unlink()
    print(f"Wrote {len(rows):,} re-evaluated rows to {output}")
    return output


def parse_args(argv: Sequence[str] | None = None):
    parser = argparse.ArgumentParser(
        description=(
            "Recompute physical widths, HiggsTools constraints, and all DM "
            "results in an existing vx=0 TRSM scan."
        )
    )
    add_cmb_arguments(parser)
    parser.add_argument("--micromegas-version", type=normalize_micromegas_version,
                        choices=MICROMEGAS_VERSIONS, default=DEFAULT_MICROMEGAS_VERSION)
    limits = parser.add_mutually_exclusive_group()
    limits.add_argument("--dm-limit-table", type=Path, help="Override the default observed LZ WS2024 SI table")
    limits.add_argument("--dm-limit-model", choices=(DEFAULT_LIMIT_MODEL,"lz2025-source", "legacy-output"),
                        help="Explicitly select an analytic SI fit instead of input provenance")
    parser.add_argument("input", type=Path, help="Existing TRSM scan TSV")
    parser.add_argument("--output", type=Path, required=True, help="New versioned TSV")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume the transactionally checkpointed <output>.partial file",
    )
    parser.add_argument(
        "--micromegas-main",
        type=Path,
        help="Path to the rebuilt canonical microOMEGAs TRSM/main executable",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=25,
        metavar="N",
        help="Commit a restartable checkpoint after N successful points (default: 25)",
    )
    args = parser.parse_args(argv)
    if args.checkpoint_every <= 0:
        parser.error("--checkpoint-every must be positive")
    return args


def resolve_evaluation_configuration(args):
    """Freeze backend, CMB settings and the input's SI treatment before any output."""
    if args.planck_cmb is None:
        partial = checkpoint_path(args.output.expanduser().resolve())
        if args.resume and partial.is_file():
            saved = read_checkpoint_metadata(partial)
            args.planck_cmb = saved.get("evaluation_configuration", {}).get("planck_cmb", {}).get("enabled", False)
        else:
            args.planck_cmb = args.micromegas_version == "7.1.4"
    backend = micromegas_configuration(args)
    require_v2_capability(backend["executable"])
    if args.planck_cmb:
        args._cmb_driver = require_cmb_capability(backend["executable"])
    # A v2 reevaluation uses the agreed v2 default, regardless of the old fit.
    model=args.dm_limit_model or DEFAULT_LIMIT_MODEL
    table_path=args.dm_limit_table or (DEFAULT_LIMIT_TABLE if model==DEFAULT_LIMIT_MODEL else None)
    table=load_si_limit_table(table_path) if table_path else None
    dd_configuration=table.metadata() if table else {"model":model}
    configuration = {"physics_version":PHYSICS_VERSION,
                     "physics_manifest":physics_manifest(backend["executable"]),
                     "micromegas":backend,"direct_detection":dd_configuration,
                     "planck_cmb":cmb_configuration(args)}
    return configuration, table, model


def run(
    argv: Sequence[str] | None = None,
    *,
    evaluator_factory: Callable[[Path | None], CoreEvaluator] = CoreEvaluator,
) -> Path:
    args = parse_args(argv)
    configuration, table, model = resolve_evaluation_configuration(args)
    evaluator = evaluator_factory(
        Path(configuration["micromegas"]["executable"]), planck_cmb=args.planck_cmb,
        limit_table=table, limit_model=model,
    )
    return reevaluate(
        args.input,
        args.output,
        evaluator,
        resume=args.resume,
        checkpoint_every=args.checkpoint_every,
        evaluation_configuration=configuration,
    )


def main() -> None:
    run()


if __name__ == "__main__":
    main()
