"""Calculation coverage, distinct from positive transition evidence."""

STATUS_COLUMNS = (
    "ewpt_transition_strengths",
    "ewpt_calctemps_returncode", "ewpt_minimatracer_returncode",
    "ewpt_execution_status", "ewpt_nlo_stability_status", "ewpt_tracing_status",
    "ewpt_coexistence_status", "ewpt_equilibrium_status", "ewpt_candidate_assessment_reason",
)


def assessment(payload):
    row = payload.get("calctemps") or {}
    if not isinstance(row, dict):
        row = {}
    complete = (row.get("status_nlo_stability") == "success" and
                row.get("status_tracing") == "success" and
                row.get("status_coex_pairs") in ("success", "no_coex_pair"))
    if not complete:
        return False,False
    if row.get("status_coex_pairs")=="no_coex_pair":
        return True,True
    import math
    indices=[key.removeprefix("status_crit_") for key in row if key.startswith("status_crit_")]
    if not indices:return False,False
    def assessed(kind,index):
        status=row.get(f"status_{kind}_{index}")
        if status=="not_met":return True
        if status!="success":return False
        for strength in payload.get("transition_strengths") or []:
            if str(strength.get("transition_index"))!=index or strength.get("temperature_kind")!=kind:continue
            try:
                return float(strength["temperature"])>0 and all(math.isfinite(float(strength[side][field]))
                    for side in ("false_vev","true_vev") for field in ("w1","wx","ws"))
            except (KeyError,ValueError,TypeError):return False
        return False
    critical_complete=all(assessed("crit",i) for i in indices)
    gw_complete=critical_complete and all(assessed(kind,i) for i in indices for kind in ("nucl","perc"))
    return critical_complete, gw_complete


def status_updates(payload):
    import json
    from trsm_inputs import json_safe
    row = payload.get("calctemps") or {}
    minima = payload.get("minimatracer") or {}
    critical, gw = assessment(payload)
    return {
        "ewpt_transition_strengths": json.dumps(json_safe(payload.get("transition_strengths") or []), separators=(",", ":"), allow_nan=False),
        "ewpt_calctemps_returncode": (payload.get("execution") or {}).get("calctemps_returncode"),
        "ewpt_minimatracer_returncode": (payload.get("execution") or {}).get("minimatracer_returncode"),
        "ewpt_execution_status": "completed" if all((payload.get("execution") or {}).get(k)==0 for k in ("calctemps_returncode","minimatracer_returncode")) else "partial_or_unknown",
        "ewpt_nlo_stability_status": row.get("status_nlo_stability", "unknown"),
        "ewpt_tracing_status": row.get("status_tracing", "unknown"),
        "ewpt_coexistence_status": row.get("status_coex_pairs", "unknown"),
        "ewpt_equilibrium_status": minima.get("equilibrium_status", "legacy_interpolated_unverified"),
        "ewpt_candidate_assessment_reason": "complete" if critical and gw else "incomplete_required_calculation",
        "ewpt_status": "assessed" if critical and gw else "incomplete",
    }
