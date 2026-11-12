"""Cross-check trsm_sidm.py against CalcHEP (micrOMEGAs) for H3 H3 -> H3 H3.

The CalcHEP driver ``DM/main_sidm.c`` must be compiled in a TRSM micrOMEGAs
project (``make main=main_sidm.c``).  For each benchmark the script

1. writes a micrOMEGAs card and runs ``main_sidm`` for several relative
   velocities, collecting the total cross section and cos(theta) bins;
2. rebuilds the point for trsm_sidm with the CalcHEP widths and the CalcHEP
   value of v, so that both codes see identical inputs;
3. compares CalcHEP with trsm_sidm's tree amplitude evaluated with
   CalcHEP's own propagator prescription (below).  Agreement at the 1e-8
   level tests the couplings, Feynman rules, kinematics, identical-particle
   factor and angular dependence; the script exits non-zero above
   ``--tolerance``;
4. prints, for information, the ratio of the physical trsm_sidm modes
   (``tree-fixed-width`` and the default unitarised ``kmatrix``) to CalcHEP.

CalcHEP propagator prescription (CalcHEP_src/c_source/num/sqme_aux.c,
``prepDen``): every propagator with a non-zero width, s- or t-channel,
carries ``i*sqrt(w2)`` with ``w = M*Gamma`` and ``q = s_or_t - M^2``::

    w2 = w^2                    if q^2 <= B^2 w^2
    w2 = (B^2 + 1) w^2 - q^2    if B^2 w^2 < q^2 <= (B^2 + 1) w^2
    w2 = 0                      otherwise,

with B the Breit-Wigner range (2.7 for processes made with ``newProcess``).
The check runs both with the CalcHEP default and with B = 1e6, where all
propagators keep their fixed widths.

Example::

    python trsm_sidm_calchep_check.py \
        --main-sidm ../runtime-v2/micromegas_7.1.4/TRSM/main_sidm
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

import trsm_sidm as sidm

# CalcHEP converts GeV^-2 to pb with 3.8937966e8 (an older (hbar c)^2); the
# same constant is used here so that the comparison tests the physics.
PB_PER_GEV2 = 3.8937966e8

SM_CARD = {  # trsm_inputs.micromegas_sm_inputs() for config/sm-inputs-v2.json
    "Mh": 125.09, "MW": 80.385, "SW": 0.47211992387942947,
    "EE": 0.30827238978402111, "Mtp": 172.5,
}

# name: (MX, Mh2, vevs, SinT, LX, LHX, LSX)
BENCHMARKS = {
    "data_par_threshold": (50.0, 100.0, 500.0, 0.3, 0.2, 0.1, 0.1),
    "generic_heavy": (200.0, 500.0, 700.0, 0.1, 0.5, 0.3, 0.8),
    "light_dm": (10.0, 40.0, 300.0, -0.2, 0.1, 0.01, 0.5),
    "h1_near_resonant": (60.0, 300.0, 400.0, 0.05, 0.3, 0.05, 0.02),
    "light_mediator_tchannel": (300.0, 20.0, 400.0, 0.02, 0.0, 0.0, 0.3),
    "negative_lx_relativistic_pole": (40.0, 90.0, 250.0, -0.1, -0.2, 0.2, 0.6),
    "near_threshold_dark_h2": (50.0, 100.001, 500.0, 0.001, 0.1, 0.0, 0.2),
}
VELOCITIES = (1.0e-4, 1.0e-3, 6.0e-3, 1.0e-2, 0.1, 0.5)
CALCHEP_DEFAULT_BWRANGE = 2.7


def write_card(path, params):
    MX, Mh2, vevs, SinT, LX, LHX, LSX = params
    card = {"LX": LX, "LHX": LHX, "LSX": LSX, "MX": MX, "vevs": vevs, "SinT": SinT, "Mh2": Mh2}
    card.update(SM_CARD)
    path.write_text("".join(f"{k}\t{v:.17g}\n" for k, v in card.items()), encoding="ascii")


def run_calchep(main_sidm, card, velocities, bins, runtime_dir, bwrange):
    env = dict(os.environ, TRSM_RUNTIME_DIR=str(runtime_dir))
    cmd = [str(main_sidm), str(card), "--bins", str(bins), "--bwrange", f"{bwrange:.17g}",
           *(f"{v:.17g}" for v in velocities)]
    text = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env,
                          cwd=Path(main_sidm).parent).stdout
    inputs = json.loads(re.search(r"^TRSM_SIDM_inputs (\{.*\})$", text, re.M).group(1))
    totals = {float(m.group(1)): float(m.group(2)) for m in
              re.finditer(r"^TRSM_SIDM_cs v=(\S+) pcm=\S+ sigma_pb=(\S+) err=0$", text, re.M)}
    parts = {}
    for m in re.finditer(r"^TRSM_SIDM_dcs v=(\S+) cmin=(\S+) cmax=(\S+) sigma_pb=(\S+) err=0$", text, re.M):
        parts.setdefault(float(m.group(1)), []).append((float(m.group(2)), float(m.group(3)), float(m.group(4))))
    return inputs, totals, parts


def point_from_calchep(inputs):
    K133, K233 = sidm.vxzero_portal_couplings(
        inputs["LHX"], inputs["LSX"], inputs["vevs"], math.asin(inputs["SinT"]), v=inputs["vev"])
    M1, M2, M3 = inputs["Mh"], inputs["Mh2"], inputs["MX"]
    g1 = inputs["width_h1"] - sidm.scalar_to_identical_scalar_width(M3, M1, K133)
    g2 = inputs["width_h2"] - sidm.scalar_to_identical_scalar_width(M3, M2, K233)
    return sidm.SIDMPoint(M3=M3, M2=M2, K133=K133, K233=K233, lX=inputs["LX"],
                          gamma1_other=max(g1, 0.0), gamma2_other=max(g2, 0.0), M1=M1)


def _calchep_width(q, w, bwrange):
    b2 = bwrange * bwrange
    q2, w2 = q * q, w * w
    out = np.where(q2 > b2 * w2, (b2 + 1.0) * w2 - q2, w2)
    return np.sqrt(np.where(q2 > (b2 + 1.0) * w2, 0.0, out))


def calchep_amplitude(point, v, x, bwrange):
    """Tree amplitude with CalcHEP's width prescription in every propagator."""
    v = np.atleast_1d(np.asarray(v, dtype=float))
    _, k2, _ = sidm.kinematics(point.M3, v)
    x = np.atleast_1d(np.asarray(x, dtype=float))[None, :]
    k2 = k2[:, None]
    t = -2.0 * k2 * (1.0 - x)
    u = -2.0 * k2 * (1.0 + x)
    amp = np.full(t.shape, -6.0 * point.lX, dtype=complex)
    for index, (Mi, Ki, _) in enumerate(point.mediators, start=1):
        w = Mi * point.total_width(index)
        qs = sidm._s_minus_mass2(point.M3, v, Mi)[:, None]
        for q in (qs, t - Mi**2, u - Mi**2):
            amp = amp - 4.0 * Ki**2 / (q + 1j * _calchep_width(q, w, bwrange))
    return amp


def calchep_like_binned(point, v, c1, c2, bwrange, n=48):
    """sigma for cos(theta) in [c1, c2] in pb, identical-particle factor included."""
    t, w = np.polynomial.legendre.leggauss(n)
    x = c1 + 0.5 * (c2 - c1) * (t + 1.0)
    s, _, _ = sidm.kinematics(point.M3, np.array([v]))
    dsig = np.abs(calchep_amplitude(point, [v], x, bwrange)[0]) ** 2 / (64.0 * math.pi**2 * s[0])
    return 0.5 * 2.0 * math.pi * 0.5 * (c2 - c1) * np.dot(w, dsig) * PB_PER_GEV2


def compare(main_sidm, bwrange, bins, tmp, verbose=True):
    worst, report = 0.0, {}
    edges = np.linspace(-1.0, 1.0, 9)
    for name, params in BENCHMARKS.items():
        card = tmp / f"{name}.par"
        write_card(card, params)
        runtime = Path(tempfile.mkdtemp(prefix=f"rt_{name}_", dir=tmp))
        inputs, totals, parts = run_calchep(main_sidm, card, VELOCITIES, bins, runtime, bwrange)
        point = point_from_calchep(inputs)
        v = np.array(sorted(totals))
        tree, _ = sidm.cross_sections(point, v, mode="tree-fixed-width")
        unit, _ = sidm.cross_sections(point, v, mode="kmatrix")
        if verbose:
            print(f"\n== {name}: MX={point.M3:g} Mh2={point.M2:g} K133={point.K133:.4g} "
                  f"K233={point.K233:.4g} lX={point.lX:g} w1={inputs['width_h1']:.4g} "
                  f"w2={inputs['width_h2']:.4g} born_h2={point.born_parameter(2):.3g}")
            print(f"   {'v/c':>8} {'CalcHEP [pb]':>14} {'emulation/CH-1':>15} {'bins max dev':>13}"
                  f" {'tree-fixed/CH':>14} {'kmatrix/CH':>12}")
        rows = []
        for vi, t_gev, u_gev in zip(v, tree, unit):
            ch = totals[vi]
            emulated = sum(calchep_like_binned(point, vi, a, b, bwrange) for a, b in zip(edges[:-1], edges[1:]))
            dev_total = emulated / ch - 1.0
            dev_bins = max((abs(calchep_like_binned(point, vi, c1, c2, bwrange) / part - 1.0)
                            for c1, c2, part in parts.get(vi, [])), default=0.0)
            worst = max(worst, abs(dev_total), dev_bins)
            row = {"v": vi, "calchep_pb": ch, "emulation_minus_1": dev_total, "bins_max_dev": dev_bins,
                   "tree_fixed_over_calchep": t_gev * PB_PER_GEV2 / ch,
                   "kmatrix_over_calchep": u_gev * PB_PER_GEV2 / ch}
            rows.append(row)
            if verbose:
                print(f"   {vi:8.1e} {ch:14.6e} {dev_total:15.2e} {dev_bins:13.2e}"
                      f" {row['tree_fixed_over_calchep']:14.6g} {row['kmatrix_over_calchep']:12.5g}")
        report[name] = {"inputs": inputs, "rows": rows}
    return worst, report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--main-sidm", required=True, help="compiled TRSM/main_sidm executable")
    parser.add_argument("--bins", type=int, default=8)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    parser.add_argument("--json", help="also write all numbers to this JSON file")
    args = parser.parse_args(argv)
    main_sidm = Path(args.main_sidm).resolve()

    report, worst = {}, 0.0
    with tempfile.TemporaryDirectory(prefix="trsm_sidm_check_") as tmp:
        for bwrange in (CALCHEP_DEFAULT_BWRANGE, 1.0e6):
            print(f"\n######## CalcHEP Breit-Wigner range B = {bwrange:g}")
            w, r = compare(main_sidm, bwrange, args.bins, Path(tmp))
            worst = max(worst, w)
            report[f"bwrange_{bwrange:g}"] = r
    print(f"\nlargest relative deviation, trsm_sidm amplitude (CalcHEP prescription) vs CalcHEP: {worst:.3e}")
    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0 if worst <= args.tolerance else 1


if __name__ == "__main__":
    sys.exit(main())
