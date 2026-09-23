#!/usr/bin/env python3
"""Create reviewable historical and cross-host comparisons from saved run evidence."""
import argparse
import csv
import hashlib
import json
import math
import shutil
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from trsm_inputs import M1, VEV, PHYSICS_VERSION
from test_trsm_theory_constraints import _copositive_quartic_matrix, _unitarity_eigenvalues


def read_tsv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_json(path):
    return json.loads(Path(path).read_text())


def verdict(value):
    return value if value in ("True", "False") else "unassessed"


def number(value):
    try:
        value = float(value)
    except (ValueError, TypeError):
        return None
    return value if math.isfinite(value) else None


def numerical_comparison(left, right, keys, rtol=2e-6, atol=1e-30):
    result = {}
    for key in keys:
        maximum = 0.0
        available = 0
        for a, b in zip(left, right, strict=True):
            x, y = number(a.get(key)), number(b.get(key))
            assert (x is None) == (y is None), (key, "coverage mismatch")
            if x is None:
                continue
            assert math.isclose(x, y, rel_tol=rtol, abs_tol=atol), (key, x, y)
            available += 1
            maximum = max(maximum, abs(x-y)/max(abs(x), abs(y), 1e-300))
        result[key] = {"compared": available, "max_relative_difference": maximum}
    return result


def tree_failure(row):
    m, vs, a, lx, hx, sx = (float(row[k]) for k in ("M2", "vs", "a12", "lX", "lPhiX", "lSX"))
    c, s = math.cos(a), math.sin(a)
    quartics = ((M1*M1*c*c + m*m*s*s)/(2*VEV*VEV),
                (M1*M1*s*s + m*m*c*c)/(2*vs*vs), lx,
                (m*m-M1*M1)*s*c/(VEV*vs), hx, sx)
    if not _copositive_quartic_matrix(*quartics):
        return "stored couplings fail current tree boundedness"
    return "tree quartic/unitarity bound; max eigenvalue=" + str(max(abs(_unitarity_eigenvalues(*quartics))))


def run(base, destination):
    base, destination = Path(base), Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    old = read_tsv(ROOT / "benchmarks/v2/stored100.tsv")
    local = read_tsv(base / "stored100-local-canonical.tsv")
    remote = read_tsv(base / "manto/stored100-manto-canonical.tsv")
    new = read_tsv(base / "stored100-ewpt/assessed100.tsv")
    assert len(old) == len(local) == len(remote) == len(new) == 100
    identity = ("M2", "M3", "vs", "vx", "a12", "lX", "lPhiX", "lSX", "selection_key")
    for rows in zip(old, local, remote, new, strict=True):
        for key in identity:
            assert all(r[key] == rows[0][key] for r in rows), (key, rows[0]["selection_key"])
    flags = ("dm", "hb", "hs", "ewpo", "wmass", "thc", "experimental_subset", "ewpt_eligible",
             "vacuum_tree_global", "rg_integration_success", "rg_bfb", "rg_unitarity", "theory_strict_subset")
    for a, b in zip(local, remote, strict=True):
        for key in flags:
            assert verdict(a.get(key)) == verdict(b.get(key)), (key, a["selection_key"])
        for key in ("point_assessment_status", "dm_calculation_status", "dm_assessment_reason"):
            assert a[key] == b[key]
    numerical = numerical_comparison(local, remote, (
        "dm_omega", "dm_xf", "dm_freezeout_temperature_GeV", "dm_dir_det", "dm_h1_width_GeV", "dm_h2_width_GeV"))
    numerical.update(numerical_comparison(local, remote, ("higgstools_hs_chi2", "higgstools_hs_delta_chi2"), 1e-8, 1e-10))
    old_flags = ("dm", "hb", "hs", "ewpo", "wmass", "thc", "evo")
    history = []
    for index, (a, b) in enumerate(zip(old, new, strict=True), 1):
        record = {"validation_index": index, **{k: a.get(k) for k in
                  ("selection_key", "selection_stratum", "source_file", "source_line", "original_point_index", *identity[:-1])}}
        notes = []
        for key in old_flags:
            record["old_" + key] = verdict(a.get(key))
            record["new_" + key] = verdict(b.get(key))
            if record["old_" + key] != record["new_" + key]:
                if record["old_" + key] == "unassessed":
                    notes.append(key + ": no stored assessment; newly evaluated")
                elif key == "thc":
                    notes.append("thc: " + tree_failure(b) + "; old source provenance does not establish why its saved flag passed")
                elif key == "dm":
                    changed = [k for k in ("dm_relic_excluded", "dm_direct_detection_excluded", "dm_indirect_detection_excluded")
                               if verdict(a.get(k)) != verdict(b.get(k))]
                    notes.append("dm: changed " + ",".join(changed) + "; shared inputs, corrected backend and official LZ table reevaluated together")
                elif key in ("hb", "hs"):
                    notes.append(key + ": recomputed Higgs widths/BR inputs and common SM reference; no one-change causal attribution")
                elif key == "wmass":
                    notes.append("wmass: recomputed from stored M2/mixing; " + b["wmass_status"])
                else:
                    notes.append(key + ": recomputed with v2 inputs")
        for key in ("dm_omega", "dm_dir_det", "dm_dir_det_limit", "dm_limit_model", "higgstools_hs_chi2", "higgstools_hs_delta_chi2"):
            record["old_" + key], record["new_" + key] = a.get(key), b.get(key)
        for key in ("experimental_subset", "ewpt_eligible", "vacuum_tree_global", "vacuum_tree_status",
                    "rg_integration_success", "rg_bfb", "rg_unitarity", "rg_first_bfb_failure_GeV", "rg_first_unitarity_failure_GeV",
                    "theory_strict_subset", "ewpt_baryo_candidate", "ewpt_gw_candidate", "ewpt_status", "ewpt_candidate_assessment_reason",
                    "dm_assessment_reason", "point_assessment_status", "point_assessment_reason"):
            record["new_" + key] = b.get(key)
        record["legacy_gw_critical_ratio"] = a.get("legacy_gw_critical_ratio")
        record["explanation"] = "; ".join(notes) or "Existing comparable verdicts agree; vacuum/RG/EWPT subset definitions are new assessments"
        history.append(record)
    with (destination / "historical-verdicts.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(history[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(history)
    runtime_local = read_json(base / "runtime-local/runtime-comparison.json")
    runtime_remote = read_json(base / "manto/runtime-comparison.json")
    assert len(runtime_local) == len(runtime_remote) == 28
    for a, b in zip(runtime_local, runtime_remote, strict=True):
        assert (a['case'], a['backend'], a['virtual_WZ']) == (b['case'], b['backend'], b['virtual_WZ'])
        assert all(a['hook_calls']) and all(b['hook_calls'])
    runtime_comparison = numerical_comparison(runtime_local, runtime_remote, ("omega", "Xf", "Tf_GeV", "SI_pb"))
    by_case = {(r['backend'], r['case'], r['virtual_WZ']): r for r in runtime_local}
    version_differences = [abs(r['omega']/by_case[('6.1.15', r['case'], r['virtual_WZ'])]['omega']-1)
                           for r in runtime_local if r['backend'] == '7.1.4']
    pa, pb = read_json(base / "point_223615-release/report.json"), read_json(base / "manto/point223615-report.json")
    assert pa['status'] == pb['status'] == 'passed'
    for key in ('ewpt_baryo_candidate', 'ewpt_gw_candidate'):
        assert pa['flags'][key] == pb['flags'][key]
    for key in ('ewpt_ew_entry_temperature_GeV', 'ewpt_gw_crit_temperature_GeV'):
        assert abs(pa['flags'][key]-pb['flags'][key]) <= .003
    for key in ('ewpt_ew_entry_jump_over_T', 'ewpt_gw_max_field_jump_over_T'):
        assert math.isclose(pa['flags'][key], pb['flags'][key], rel_tol=1e-3, abs_tol=1e-5)
    for p in (pa, pb):
        assert p['crossing']['high_T_GeV']-p['crossing']['low_T_GeV'] <= .001
    amplitudes = [read_json(base / p) for p in ('dm-amplitudes-final/independent-amplitude-report.json', 'manto/amplitudes.json')]
    assert all(len(values) == 22 and all(r['status'] == 'passed' for r in values) for values in amplitudes)
    report = {
        'status': 'passed', 'physics_version': PHYSICS_VERSION,
        'scope': '100 zero-temperature rows on both hosts; all 22 eligible thermal rows locally at Tmax=300 GeV; point223615 and integrated pilot on both hosts',
        'strata': dict(Counter(r['selection_stratum'] for r in old)),
        'counts': {k: dict(Counter(verdict(r.get(k)) for r in new)) for k in (*flags, 'ewpt_baryo_candidate', 'ewpt_gw_candidate')},
        'historical_transitions': {k: dict(Counter(verdict(a.get(k))+' -> '+verdict(b.get(k)) for a,b in zip(old,new))) for k in old_flags},
        'cross_host_stored100': numerical, 'cross_host_28_runtime_cases': runtime_comparison,
        'maximum_backend_version_relative_omega_difference': max(version_differences),
        'W_below_virtual_on_over_off': by_case[('7.1.4','W_below',True)]['omega_virtual_on_over_off'],
        'point223615': {'local': pa, 'manto': pb},
        'tolerances': {'DM_cross_host_rtol':2e-6, 'DM_cross_host_atol':1e-30, 'HS_rtol':1e-8, 'HS_atol':1e-10,
                       'thermal_temperature_atol_GeV':.003, 'thermal_ratio_rtol':.001,
                       'loop_amplitude_rtol':2e-8, 'on_shell_width_rtol':2e-6, 'crossing_bracket_GeV':.001},
        'input_sha256': {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (
            ROOT/'benchmarks/v2/stored100.tsv',base/'stored100-local-canonical.tsv',base/'manto/stored100-manto-canonical.tsv',base/'stored100-ewpt/assessed100.tsv')},
    }
    for name, path in {'pilot_local':'pilot-integrated-local-release/report.json', 'pilot_manto':'manto/pilot-report.json',
                       'audit_local':'audit-comparison-local.json','audit_manto':'manto/audit-comparison-manto.json'}.items():
        payload = read_json(base/path)
        assert payload['status'] == 'passed', name
        report[name] = payload
    (destination/'comparison-summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    shutil.copy2(base/'stored100-ewpt/assessed100.tsv',destination/'assessed100.tsv')
    for source, target in [('dm-amplitudes-final/independent-amplitude-report.json','amplitudes-local.json'),
                           ('manto/amplitudes.json','amplitudes-manto.json'), ('runtime-local/runtime-comparison.json','runtime-local.json'),
                           ('manto/runtime-comparison.json','runtime-manto.json')]:
        shutil.copy2(base/source,destination/target)
    print(json.dumps({k:v for k,v in report.items() if k in ('status','counts','historical_transitions','maximum_backend_version_relative_omega_difference')},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=Path,default=ROOT/'validation/v2')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'validation/v2/release')
    args=parser.parse_args();run(args.base.resolve(),args.output_dir.resolve())
