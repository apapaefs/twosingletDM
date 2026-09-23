"""Versioned scientific provenance and independent assessment subsets."""

import hashlib
import json
import math
import subprocess
import os
import importlib.metadata
from functools import lru_cache
from pathlib import Path
from trsm_inputs import PHYSICS_VERSION, SCHEMA_VERSION, M1, sm_inputs, RG_BOUNDARY, nullable_and

PROFILE_COLUMNS = (
    "constraint_version", "constraint_schema_version", "point_index",
    "point_assessment_status", "point_assessment_reason",
    "experimental_subset", "dm_subset", "ewpt_eligible",
    "wmass_covered", "wmass_status",
    "h1_width_over_mass", "h2_width_over_mass", "h1_ctau_mm", "h2_ctau_mm",
)


def profile_updates(point):
    m2 = float(point["M2"])
    covered = 133 <= m2 <= 999
    experimental = nullable_and(point.get(k) for k in ("hb", "hs", "ewpo"))
    if covered:
        experimental = nullable_and((experimental, point.get("wmass")))
    result = {"constraint_version": PHYSICS_VERSION, "constraint_schema_version": SCHEMA_VERSION,
              "experimental_subset": experimental, "dm_subset": point.get("dm"),
              "wmass_covered": covered, "wmass_status": ("assessed" if point.get("wmass") is not None else "nonfinite_or_failed") if covered else "outside_coverage",
              "ewpt_eligible": point.get("thc") is True and experimental is not False}
    for i, mass in ((1, M1), (2, m2)):
        width = point.get(f"w{i}")
        valid = isinstance(width, (int, float)) and math.isfinite(width) and width >= 0
        result[f"h{i}_width_over_mass"] = width/mass if valid else None
        result[f"h{i}_ctau_mm"] = 1.973269804e-13/width if valid and width > 0 else None
    unassessed=[key for key in ("thc","hb","hs","ewpo","dm","vacuum_tree_global","rg_bfb","rg_unitarity") if point.get(key) is None]
    if covered and point.get("wmass") is None:unassessed.append("wmass")
    result["point_assessment_status"]="partially_assessed" if unassessed else "assessed"
    result["point_assessment_reason"]="unassessed: "+", ".join(unassessed) if unassessed else "independent_zero_temperature_subsets_assessed"
    return result


def sha256(path):
    path = Path(path)
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_revision(path):
    try:
        return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


@lru_cache(maxsize=16)
def physics_manifest(micromegas_executable, calctemps_executable=None, minima_executable=None):
    root = Path(__file__).resolve().parent
    source_names = ("generate_trsm_points", "generate_trsm_info", "test_trsm_DM", "test_trsm_ewpt",
        "test_trsm_evolution", "test_trsm_theory_constraints", "test_trsm_higgstools", "singlet_EWPO",
        "trsm_inputs", "trsm_constraint_profile", "trsm_theory_diagnostics", "trsm_direct_detection",
        "trsm_micromegas", "trsm_cmb", "trsm_scan_campaign", "trsm_kstoalphas", "scan_output",
        "ewpt_assessment", "ewpt_entry_criterion", "ewpt_equilibrium", "ewpt_x_history",
        "dm_thermal_relic_diagnostic", "reevaluate_trsm_dm_higgs", "reprocess_trsm_ewpt",
        "run_trsm_seed_campaign", "mg5_process_runner", "generate_mg5_trsm_xsecs")
    paths = [root/(name+".py") for name in source_names] + list((root/"DM/models/h4GOn").glob("*.mdl"))
    paths += [root/"DM/main.c", root/"DM/trsm_loop.c", root/"DM/models/lanhep_mdl/TRSM_mixed.mdl"]
    paths += list((root/"config").glob("*.json")) + list((root/"DM/data").rglob("*.json"))
    paths += list((root/"BSMPT").rglob("*.cpp")) + list((root/"BSMPT").rglob("*.h"))
    paths += list((root/"datafiles").glob("*"))
    for directory in ("YR", "couplings_vxzero"):
        paths += [p for p in (root/directory).rglob('*') if p.is_file() and p.suffix in ('.dat','.txt')]
    hashes = {str(p.relative_to(root)): sha256(p) for p in sorted(paths) if p.is_file()}
    executables = {"micromegas": str(micromegas_executable)}
    if calctemps_executable:
        executables["CalcTemps"] = str(calctemps_executable)
    if minima_executable:
        executables["MinimaTracer"] = str(minima_executable)
    if calctemps_executable:
        executables["PhaseProbe"] = str(Path(calctemps_executable).with_name("PhaseProbe"))
    manifest = {"physics_version": PHYSICS_VERSION, "schema_version": SCHEMA_VERSION,
                "repository_commit": git_revision(root), "sources_sha256": hashes,
                "sm_inputs": sm_inputs(), "rg_boundary": RG_BOUNDARY,
                "executables": {k: {"path": v, "sha256": sha256(v)} for k, v in executables.items()},
                "runtime_build_manifest_sha256": sha256(Path(micromegas_executable).resolve().parents[2]/"runtime-manifest.json") if len(Path(micromegas_executable).resolve().parents)>2 else None,
                "python_packages": {name:package_version(name) for name in ("HiggsTools","numpy","scipy")},
                "virtual_WZ_decays": "TRSM_LEGACY_VIRTUAL_OFF" not in os.environ,
                "datasets": {name: {"commit":git_revision(root.parent/name),"tree_sha256":dataset_digest(root.parent/name)} for name in ("hbdataset", "hsdataset")},
                "experimental_prescriptions": {
                    "HiggsSignals":"pipeline SM reference, delta chi2 < 4; heuristic fixed threshold",
                    "STU":"Jens Erler private communication 2025-05-16; 3 dof chi2 <= 7.82",
                    "W_mass":"Tania Robens Snowmass table, cubic interpolation, no extrapolation (133--999 GeV)"},
                "thermal_settings":{"crossing_bracket_GeV":.001,"EW_field_threshold_GeV":5,"X_field_threshold_GeV":1,
                                    "S_field_threshold_GeV":1,"freezeout_window":[.5,2],"VEV_fractional_screen":.1}}
    manifest["sha256"] = hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest()
    return manifest


def package_version(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def precision_updates(m2, angle):
    import singlet_EWPO as e
    from singlet_EWPO import check_wmass_tania
    try:
        ewpo=e.check_EWPO_wU(M1,m2,math.sin(angle),e.Mz,e.Mw,
            e.Delta_S_central_wU,e.Delta_T_central_wU,e.Delta_U_central_wU,
            e.errS_wU,e.errT_wU,e.errU_wU,e.covST_wU,e.covSU_wU,e.covTU_wU)
    except (ValueError,ArithmeticError):ewpo=None
    return {'ewpo':ewpo,'wmass':check_wmass_tania(m2,math.sin(angle))}


def validate_v2_result(point):
    """Check status/value consistency without coercing missing results to failure."""
    for key in ('dm','hb','hs','ewpo','thc','experimental_subset','vacuum_tree_global',
                'rg_bfb','rg_unitarity','ewpt_baryo_candidate','ewpt_gw_candidate'):
        if point.get(key) is not None and type(point[key]) is not bool:
            raise ValueError(f'{key} must be bool or None')
    if point.get('dm_calculation_status')!='success' and point.get('dm') is not None:
        raise ValueError('Unassessed solver result has a DM verdict')
    if point.get('dm_direct_detection_available') is False and point.get('dm_direct_detection_excluded') is not None:
        raise ValueError('Uncovered DD result has an exclusion verdict')
    if point.get('dm_calculation_status')=='success':
        for key in ('dm_omega','dm_dir_det','dm_mdm'):
            value=point.get(key)
            if not isinstance(value,(float,int)) or not math.isfinite(value) or value<0:
                raise ValueError(f'Invalid successful DM value: {key}')
    xf,tf=point.get('dm_xf'),point.get('dm_freezeout_temperature_GeV')
    if (xf is None)!=(tf is None) or (xf is not None and (xf<=0 or not math.isclose(point['dm_mdm']/xf,tf,rel_tol=1e-12))):
        raise ValueError('Inconsistent freeze-out diagnostics')
    for width,br,total in [('h1_h3h3_width','h1_h3h3_br','w1'),('h2_h3h3_width','h2_h3h3_br','w2')]:
        if point.get(total,0)>0 and not math.isclose(point[width]/point[total],point[br],rel_tol=1e-10,abs_tol=1e-12):
            raise ValueError('Inconsistent width and branching fraction')


@lru_cache(maxsize=8)
def dataset_digest(path):
    files=sorted(p for p in Path(path).rglob('*') if p.is_file() and '.git' not in p.parts)
    if not files:return None
    h=hashlib.sha256()
    for p in files:h.update((str(p.relative_to(path))+':'+sha256(p)+'\n').encode())
    return h.hexdigest()
