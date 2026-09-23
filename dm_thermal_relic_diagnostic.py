"""Screen fixed-vacuum relic inputs against the sampled global thermal branch.

These observables describe equilibrium minima and zero-temperature mass gaps.
They do not establish the cosmological phase history or validate a relic density.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Mapping


from trsm_inputs import M1 as M1_GEV
VEV_SHIFT_SCREEN = 0.10
WINDOW_LOW_FACTOR = 0.5
WINDOW_HIGH_FACTOR = 2.0

RESONANCE_COLUMNS = (
    "dm_resonance_h1_mass_gap_GeV",
    "dm_resonance_h2_mass_gap_GeV",
    "dm_resonance_h1_abs_gap_over_Tf",
    "dm_resonance_h2_abs_gap_over_Tf",
    "dm_resonance_nearest_mediator",
    "dm_resonance_h1_abs_gap_over_width",
    "dm_resonance_h2_abs_gap_over_width",
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


def resonance_proximity_updates(m2, m3, freezeout_temperature=None, *, m1=M1_GEV, widths=(None, None)):
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
    for name, width in zip(("h1", "h2"), widths):
        width = finite_number(width)
        if width is not None and width > 0:
            updates[f"dm_resonance_{name}_abs_gap_over_width"] = abs(gaps[name]) / width
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
        if temperature is None or temperature < 0:
            return None
        if ew_vev is None or s_vev is None:
            ew_vev = s_vev = math.nan
        phase = point.get("phase_index")
        if phase is None and "phase_index" not in point:
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
        return samples[position][1:3] if samples[position][3] is not None else None
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
    if samples[0][3] is None or not all(math.isfinite(x) for x in zero) or zero[0] <= 1e-9 or zero[1] <= 1e-9:
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
    if full_coverage and all(sample[3] is not None for sample in samples if low <= sample[0] <= high):
        updates["dm_relic_thermal_phase_boundary_bracket_overlaps_window"] = boundary

    endpoint_vevs = [_vevs_at(samples, low), _vevs_at(samples, high)]
    vevs = [sample[1:3] for sample in samples if low <= sample[0] <= high and sample[3] is not None]
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


THERMAL_COUPLING_COLUMNS = (
    "dm_thermal_coupling_status",
    *tuple(f"dm_thermal_{name}_{suffix}" for name in ("K133", "K233") for suffix in
           ("T0_GeV", "Tf_GeV", "max_abs_change_GeV", "max_fractional_change",
            "T0_cancellation_ratio", "minimum_cancellation_ratio", "cancellation_sensitive")),
)


def coupling_sensitivity_updates(payload, point):
    """Fixed-zero-T-mixing proxies, not thermally rediagonalized couplings.

    K133=(lambda_HX*w1*cos(a)-lambda_SX*ws*sin(a))/2, with the analogous
    orthogonal K233. Report absolute changes even if K(0) cancels to zero.
    The cancellation ratio |sum terms|/sum |terms| has a 0.1 screening cut.
    """
    from trsm_inputs import VEV
    result = dict.fromkeys(THERMAL_COUPLING_COLUMNS)
    result['dm_thermal_coupling_status']='unassessed'
    tf=finite_number(point.get('dm_freezeout_temperature_GeV'))
    inputs=[finite_number(point.get(k)) for k in ('lPhiX','lSX','a12','vs')]
    if tf is None or tf<=0 or any(x is None for x in inputs):return result
    hp,sp,angle,vs=inputs; c,s=math.cos(angle),math.sin(angle)
    branch=(payload.get('minimatracer') or {}).get('global_branch') or []
    samples=_global_samples(payload)
    if not samples:return result
    low,high=tf/2,tf*2
    indices=[i for i,row in enumerate(samples) if low<=row[0]<=high]
    # Add interpolation neighbours at both window boundaries and Tf.
    for t in (low,tf,high):
        i=bisect.bisect_left([row[0] for row in samples],t)
        indices.extend(j for j in (i-1,i) if 0<=j<len(samples))
    indices=sorted(set(indices)); by_temp={finite_number(row.get('temp')):row for row in branch}
    if any(samples[i][3] is None for i in indices):return result
    if any(finite_number(by_temp[samples[i][0]].get('wx')) is None or
           abs(float(by_temp[samples[i][0]]['wx']))>=1 for i in indices):
        result['dm_thermal_coupling_status']='incompatible_or_unknown_X_branch';return result
    if any(samples[i][1]<5 for i in indices):
        result['dm_thermal_coupling_status']='outside_EW_broken_branch';return result
    endpoints=[_vevs_at(samples,t) for t in (low,tf,high)]
    if any(v is None for v in endpoints) or not indices:return result
    if len({samples[i][3] for i in indices})!=1:
        result['dm_thermal_coupling_status']='phase_boundary_in_window';return result
    values=[row[1:3] for row in samples if low<=row[0]<=high]+endpoints
    for name,projection in [('K133',(c,-s)),('K233',(s,c))]:
        def calc(v):
            terms=(hp*v[0]*projection[0]/2,sp*v[1]*projection[1]/2)
            norm=sum(abs(x) for x in terms); k=sum(terms)
            return k,abs(k)/norm if norm>0 else None,norm
        k0,ratio0,norm0=calc((VEV,vs)); ktf,_,_=calc(endpoints[1])
        pairs=[calc(v) for v in values];change=max(abs(k-k0) for k,_,_ in pairs)
        ratios=[r for _,r,_ in pairs if r is not None]
        prefix=f'dm_thermal_{name}_'
        result.update({prefix+'T0_GeV':k0,prefix+'Tf_GeV':ktf,prefix+'max_abs_change_GeV':change,
                       prefix+'max_fractional_change':change/abs(k0) if abs(k0)>1e-12*max(norm0,1) else None,
                       prefix+'T0_cancellation_ratio':ratio0,
                       prefix+'minimum_cancellation_ratio':min(ratios) if ratios else None,
                       prefix+'cancellation_sensitive':min(ratios+[ratio0])<.1 if ratios and ratio0 is not None else None})
    result['dm_thermal_coupling_status']='fixed_mixing_screen_assessed'
    return result


def thermal_input_updates(payload,point):
    result=thermal_vev_updates(payload,point.get('dm_freezeout_temperature_GeV'))
    result.update(coupling_sensitivity_updates(payload,point))
    result.update(resonance_proximity_updates(point.get('M2'),point.get('M3'),
        point.get('dm_freezeout_temperature_GeV'),
        widths=(point.get('dm_h1_width_GeV'),point.get('dm_h2_width_GeV'))))
    return result
