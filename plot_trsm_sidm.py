"""Plots of the H3 self-scattering cross section (needs matplotlib).

Three figures, each written as PNG and PDF:

``landscape``  sigma_V/m against M2/M3 at fixed M3 for several K233, at a
               dwarf and a cluster velocity.  The range covers the light
               mediator region M3 > M2 and the s-channel pole at M2 = 2 M3.
``velocity``   sigma_V/m against v for benchmarks with M3 > M2, comparing
               the Schroedinger phase shifts with the unitarised Born
               approximation.
``scan``       the points of a scan TSV: model value at cluster and dwarf
               velocities against M2/M3, marked by regime and abundance.
``points``     velocity-averaged sigma_V/m against the mean relative velocity
               for a few selected points of a scan, with the lower limit for
               galactic cores and the group and cluster upper limits.

Examples::

    python plot_trsm_sidm.py landscape --m3 60 --k233 1 10 50 120 --output-dir plots/sidm
    python plot_trsm_sidm.py velocity --output-dir plots/sidm
    python plot_trsm_sidm.py scan output/trsm_points_X.dat --output-dir plots/sidm
    python plot_trsm_sidm.py points output/trsm_points_X.dat --output-dir plots/sidm
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import trsm_sidm as sidm

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")   # fixed categorical order
MARKERS = ("o", "s", "^")
C = sidm.C_KMS

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.linewidth": 0.8, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "legend.frameon": False, "legend.labelcolor": INK,
})


def _sci(value):
    """Mathtext for a number, with powers of ten written out."""
    mantissa, exponent = f"{value:.2e}".split("e")
    mantissa, exponent = float(mantissa), int(exponent)
    if -3 <= exponent <= 3:
        return f"{value:g}"
    return f"10^{{{exponent}}}" if mantissa == 1.0 else f"{mantissa:g}\\times 10^{{{exponent}}}"


def _style(ax):
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _save(fig, directory, stem, formats):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for extension in formats:
        path = directory / f"{stem}.{extension}"
        fig.savefig(path, dpi=170)
        paths.append(path)
    plt.close(fig)
    print("wrote", ", ".join(str(p) for p in paths))
    return paths


def yukawa_thresholds(b_max, M3=50.0, m=5.0):
    """b = alpha M3/m at which s-wave zero-energy states appear, up to b_max."""
    k = np.array([1.0e-7 * m])

    def denominator(b):
        return sidm.yukawa_tan_delta(M3, [b * m / M3], [m], k, [0])[1][0, 0]

    out, grid = [], np.geomspace(0.5, max(b_max, 0.6), int(60 * math.log(max(b_max, 0.6) / 0.5)) + 2)
    values = [denominator(b) for b in grid]
    for lo, hi, flo, fhi in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if flo * fhi < 0.0:
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                fmid = denominator(mid)
                lo, hi, flo = (mid, hi, fmid) if flo * fmid > 0.0 else (lo, mid, flo)
            out.append(0.5 * (lo + hi))
    return out


def _clustered(centre, lo, hi, n=60, closest=1.0e-12, widest=0.3):
    offsets = np.geomspace(closest, widest, n)
    points = np.concatenate([centre * (1.0 - offsets), [centre], centre * (1.0 + offsets)])
    return points[(points > lo) & (points < hi)]


def landscape_curve(M3, K233, v_kms, ratio_lo, ratio_hi, *, lX, gamma2_other, n=260, yukawa="auto"):
    """sigma_V/m [cm^2/g] against M2/M3 at fixed relative velocity."""
    v = v_kms / C
    ratios = [np.geomspace(ratio_lo, ratio_hi, n)]
    ratios.append(_clustered(2.0 / math.sqrt(1.0 - 0.25 * v * v), ratio_lo, ratio_hi))     # sqrt(s) = M2
    alpha = K233**2 / (4.0 * math.pi * M3**2)
    for b in yukawa_thresholds(alpha / ratio_lo):
        ratios.append(_clustered(alpha / b, ratio_lo, ratio_hi, closest=1.0e-9))
    ratios = np.unique(np.concatenate(ratios))
    values = np.empty_like(ratios)
    for i, ratio in enumerate(ratios):
        point = sidm.SIDMPoint(M3=M3, M2=ratio * M3, K133=0.0, K233=K233, lX=lX, gamma2_other=gamma2_other)
        _, sigma_v = sidm.cross_sections(point, [v], yukawa=yukawa)
        values[i] = sidm.per_mass_cm2_g(sigma_v[0], M3)
    return ratios, values


def run_landscape(args):
    fig, axes = plt.subplots(1, len(args.velocities), figsize=(5.6 * len(args.velocities), 4.6), sharey=True)
    axes = np.atleast_1d(axes)
    lo, hi = args.ratio_min, args.ratio_max
    for ax, v_kms in zip(axes, args.velocities):
        for color, K233 in zip(SERIES, args.k233):
            ratios, values = landscape_curve(args.m3, K233, v_kms, lo, hi, lX=args.lx,
                                             gamma2_other=args.gamma2_other)
            alpha = K233**2 / (4.0 * math.pi * args.m3**2)
            ax.plot(ratios, values, lw=1.6, color=color, solid_joinstyle="round",
                    label=rf"$K_{{233}}={K233:g}$ GeV ($\alpha={_sci(float(f'{alpha:.2g}'))}$)")
        ceiling = float(sidm.unitarity_limit_cm2_g(args.m3, v_kms))
        ax.axhline(ceiling, color=INK2, lw=1.1, ls=(0, (4, 3)))
        ax.annotate(r"s-wave unitarity $8\pi/(k^2 m)$", (hi, ceiling), xytext=(-4, 4), textcoords="offset points",
                    ha="right", va="bottom", fontsize=8.5, color=INK2)
        ax.axvspan(lo, 1.0, color=INK2, alpha=0.06, lw=0)
        ax.annotate(r"$M_3 > M_2$", (math.sqrt(lo * 1.0), 0.02), xycoords=("data", "axes fraction"),
                    ha="center", fontsize=9, color=INK2)
        ax.axvline(2.0, color=GRID, lw=1.0, zorder=0)
        ax.annotate(r"$M_2 = 2M_3$", (2.0, 0.02), xycoords=("data", "axes fraction"), xytext=(4, 0),
                    textcoords="offset points", ha="left", fontsize=9, color=INK2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_xlabel(r"$M_2/M_3$")
        ax.set_title(rf"$v = {v_kms:g}$ km/s", loc="left", fontsize=10, color=INK)
        _style(ax)
    axes[0].set_ylabel(r"$\sigma_V/m_{\rm DM}$ [cm$^2$/g]")
    axes[0].set_ylim(args.ymin, args.ymax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=9, loc="lower center", ncol=len(handles))
    fig.suptitle(rf"$M_3 = {args.m3:g}$ GeV, $\lambda_X = {args.lx:g}$, $K_{{133}} = 0$, "
                 rf"$\Gamma_{{2,\rm other}} = {_sci(args.gamma2_other)}$ GeV", x=0.01, ha="left", fontsize=10, color=INK2)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    return _save(fig, args.output_dir, f"sidm_landscape_m3_{args.m3:g}", args.formats)


def default_benchmarks(M3=60.0, K233=120.0):
    """M3 > M2 benchmarks: Born regime, then approaching the first bound state."""
    alpha = K233**2 / (4.0 * math.pi * M3**2)
    resonance = yukawa_thresholds(3.0)[0]

    def point(M2, K):
        return sidm.SIDMPoint(M3=M3, M2=M2, K133=0.0, K233=K, lX=0.1)

    at = alpha * M3 / resonance
    return [
        (r"$M_2=20$ GeV, $K_{233}=10$ GeV ($b=0.007$, Born regime)", point(20.0, 10.0)),
        (rf"$M_2=14$ GeV, $K_{{233}}={K233:g}$ GeV ($b={alpha * M3 / 14.0:.2f}$)", point(14.0, K233)),
        (rf"$M_2={at * 1.001:.3f}$ GeV, $K_{{233}}={K233:g}$ GeV ($b$ 0.1% below ${resonance:.3f}$)",
         point(at * 1.001, K233)),
        (rf"$M_2={at:.3f}$ GeV, $K_{{233}}={K233:g}$ GeV ($b={resonance:.3f}$, bound state at threshold)",
         point(at, K233)),
    ]


def run_velocity(args):
    benchmarks = default_benchmarks(args.m3, args.k233)
    v_kms = np.geomspace(args.vmin, args.vmax, args.n)
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for color, (label, point) in zip(SERIES, benchmarks):
        _, exact = sidm.cross_sections(point, v_kms / C)
        _, born = sidm.cross_sections(point, v_kms / C, yukawa="off")
        ax.plot(v_kms, sidm.per_mass_cm2_g(exact, point.M3), lw=2, color=color, label=label)
        ax.plot(v_kms, sidm.per_mass_cm2_g(born, point.M3), lw=1.2, color=color, ls=(0, (2, 2)))
    ax.plot(v_kms, sidm.unitarity_limit_cm2_g(args.m3, v_kms), lw=1.1, color=INK, ls=(0, (4, 3)))
    ax.annotate(r"s-wave unitarity $8\pi/(k^2 m)$", (v_kms[0] * 3.0, float(sidm.unitarity_limit_cm2_g(args.m3, v_kms[0] * 3.0))),
                xytext=(8, 6), textcoords="offset points", ha="left", fontsize=8.5, color=INK2)
    ax.plot([], [], lw=1.2, color=INK2, ls=(0, (2, 2)), label="dotted: same point with unitarised Born exchange")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(v_kms[0], v_kms[-1])
    ax.set_xlabel("relative velocity $v$ [km/s]")
    ax.set_ylabel(r"$\sigma_V/m_{\rm DM}$ [cm$^2$/g]")
    ax.set_title(rf"Light mediator, $M_3 = {args.m3:g}$ GeV $> M_2$ ($\lambda_X = 0.1$, $K_{{133}} = 0$)",
                 loc="left", fontsize=10, color=INK)
    _style(ax)
    ax.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1)
    fig.tight_layout()
    return _save(fig, args.output_dir, f"sidm_velocity_m3_{args.m3:g}", args.formats)


def run_scan(args):
    import trsm_sidm_bounds as bounds

    limits = {limit.id: limit for limit in bounds.load_limits(args.limits, "all")[0]}
    panels = [limits[name] for name in args.panels]
    _, rows = sidm.read_scan(args.input)
    data = {"ratio": [], "regime": [], "xi": [], **{limit.id: [] for limit in panels}}
    for row in rows:
        try:
            point = sidm.SIDMPoint.from_scan_row(row)
            values, _ = bounds.model_values(point, panels)
        except Exception:
            continue
        data["ratio"].append(point.M2 / point.M3)
        data["regime"].append(point.regime())
        data["xi"].append(bounds.abundance_fraction(sidm._row_float(row, args.omega_column)))
        for limit in panels:
            data[limit.id].append(values[limit.id])
    ratio, regime, xi = np.array(data["ratio"]), np.array(data["regime"]), np.array(data["xi"])
    fig, axes = plt.subplots(1, len(panels), figsize=(5.6 * len(panels), 4.6), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, limit in zip(axes, panels):
        values = np.array(data[limit.id])
        for color, marker, name in zip(SERIES, MARKERS, ("perturbative", "nonperturbative", "strong")):
            for filled in (False, True):
                pick = (regime == name) & ((xi >= args.min_fraction) == filled)
                if not pick.any():
                    continue
                ax.scatter(ratio[pick], values[pick], s=30 if filled else 22, marker=marker,
                           facecolor=color if filled else "none", edgecolor=SURFACE if filled else color,
                           linewidth=0.8, zorder=3 if filled else 2)
        ax.axhline(limit.limit_cm2_g, color=INK, lw=1.2)
        ax.annotate(f"limit {limit.limit_cm2_g:g} cm$^2$/g ({limit.confidence_level})", (0.99, limit.limit_cm2_g),
                    xycoords=("axes fraction", "data"), xytext=(0, -5), textcoords="offset points",
                    ha="right", va="top", fontsize=8.5, color=INK)
        ax.axvspan(ratio.min() / 1.3, 1.0, color=INK2, alpha=0.06, lw=0)
        ax.annotate(r"$M_3 > M_2$", (0.02, 0.02), xycoords="axes fraction", fontsize=9, color=INK2)
        ax.axvline(2.0, color=GRID, lw=1.0, zorder=0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(ratio.min() / 1.3, ratio.max() * 1.3)
        ax.set_xlabel(r"$M_2/M_3$")
        where = (rf"$\langle v\rangle = {limit.mean_relative_velocity_kms:g}$ km/s" if limit.kind == "maxwell"
                 else limit.describe_velocity())
        ax.set_title(f"{limit.label}\n{where}", loc="left", fontsize=9.5, color=INK)
        _style(ax)
    axes[0].set_ylabel(r"$\sigma_V/m_{\rm DM}$ [cm$^2$/g], before abundance rescaling")
    handles = [plt.Line2D([], [], ls="none", marker=m, mfc=c, mec=SURFACE, ms=7, label=name)
               for c, m, name in zip(SERIES, MARKERS, ("perturbative", "non-perturbative exchange", "strong coupling (unreliable)"))]
    handles.append(plt.Line2D([], [], ls="none", marker="o", mfc="none", mec=INK2, ms=6,
                              label=rf"open: abundance fraction $< {args.min_fraction:g}$ (no bound applies)"))
    fig.legend(handles=handles, fontsize=8.5, loc="lower center", ncol=2)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    return _save(fig, args.output_dir, Path(args.input).name.split(".")[0] + "_sidm", args.formats)


def important_rows(rows, omega_column="dm_omega", min_fraction=0.9, M1=sidm.DEFAULT_M1):
    """Pick a few rows of a scan that matter for self-interactions.

    Returns a list of (description, row index).  The criteria are: the largest
    sigma_V/m at 100 km/s overall, among reliable points, among points where
    H3 is the dark matter (abundance fraction >= min_fraction) and among
    reliable points with M3 > M2; and the points closest to the s-channel
    poles M2 = 2 M3 and M1 = 2 M3.
    """
    import trsm_sidm_bounds as bounds

    data = []
    for index, row in enumerate(rows):
        try:
            point = sidm.SIDMPoint.from_scan_row(row, M1=M1)
            _, sigma_v = sidm.cross_sections(point, [100.0 / C])
        except Exception:
            continue
        xi = bounds.abundance_fraction(sidm._row_float(row, omega_column))
        data.append((index, point, float(sigma_v[0]) / point.M3, point.regime() != "strong", xi))
    if not data:
        raise SystemExit("no usable rows")

    def best(key, keep=lambda item: True):
        pool = [item for item in data if keep(item)]
        return max(pool, key=key)[0] if pool else None

    criteria = [
        ("largest overall", best(lambda d: d[2])),
        ("largest reliable", best(lambda d: d[2], lambda d: d[3])),
        (rf"largest with $\xi\geq{min_fraction:g}$", best(lambda d: d[2], lambda d: d[3] and d[4] >= min_fraction)),
        (r"closest to $M_2=2M_3$", best(lambda d: -abs(d[1].M2 - 2.0 * d[1].M3))),
        (r"closest to $M_1=2M_3$", best(lambda d: -abs(d[1].M1 - 2.0 * d[1].M3))),
        (r"largest reliable with $M_3>M_2$", best(lambda d: d[2], lambda d: d[3] and d[1].M3 > d[1].M2)),
    ]
    chosen = {}
    for description, index in criteria:
        if index is not None:
            chosen.setdefault(index, []).append(description)
    return [("; ".join(names), index) for index, names in chosen.items()]


def run_points(args):
    import trsm_sidm_bounds as bounds
    import trsm_sidm_lower_limit as lower

    _, rows = sidm.read_scan(args.input)
    if args.rows:
        selection = [(f"row {index}", index) for index in args.rows]
    else:
        selection = important_rows(rows, args.omega_column, args.min_fraction)
    selection = selection[:len(SERIES)]
    mean_v = np.geomspace(args.vmin, args.vmax, args.n)
    dispersions = [sidm.mean_relative_velocity_to_dispersion(v) for v in mean_v]

    def curve(point):
        averages = sidm.maxwell_averages(point, dispersions, [1])
        return sidm.per_mass_cm2_g(np.array([averages[(d, 1)] for d in dispersions]), point.M3)

    fig, ax = plt.subplots(figsize=(8.4, 6.6))
    for color, (description, index) in zip(SERIES, selection):
        point = sidm.SIDMPoint.from_scan_row(rows[index])
        xi = bounds.abundance_fraction(sidm._row_float(rows[index], args.omega_column))
        flag = ", strong coupling" if point.regime() == "strong" else ""
        label = (rf"{description}: $M_2={point.M2:.4g}$, $M_3={point.M3:.4g}$, $K_{{133}}={point.K133:.2g}$, "
                 rf"$K_{{233}}={point.K233:.3g}$ GeV, $\xi={_sci(float(f'{xi:.1g}'))}${flag}")
        ax.plot(mean_v, curve(point), lw=2, color=color, label=label)
    if args.tuned_example:
        M3, K233, binding_keV = 20.0, 5.0, 0.8
        tuned = sidm.SIDMPoint(M3=M3, M2=2.0 * M3 - binding_keV * 1.0e-6, K133=0.0, K233=K233, lX=0.1)
        ax.plot(mean_v, curve(tuned), lw=1.6, color=INK, ls=(0, (5, 2)),
                label=rf"tuned example, not in the scan: $M_3={M3:g}$ GeV, $2M_3-M_2={binding_keV:g}$ keV, "
                      rf"$K_{{233}}={K233:g}$ GeV")

    requirement = lower.load_requirements(args.lower_table)[0][0]
    lo, hi = requirement.velocity_range_kms
    ax.fill_between([lo, hi], requirement.minimum_cm2_g, args.ymax, color=INK2, alpha=0.10, lw=0)
    ax.plot([lo, hi], [requirement.minimum_cm2_g] * 2, color=INK, lw=1.2)
    ax.annotate(rf"galactic cores: $\sigma/m \geq {requirement.minimum_cm2_g:g}$ cm$^2$/g" "\n"
                rf"for $\langle v\rangle = {lo:g}$-${hi:g}$ km/s",
                (math.sqrt(lo * hi), args.ymax), xytext=(0, -6), textcoords="offset points",
                ha="center", va="top", fontsize=8.5, color=INK)
    for limit in bounds.load_limits(args.limits)[0]:
        if limit.kind != "maxwell":
            continue
        ax.annotate("", (limit.mean_relative_velocity_kms, limit.limit_cm2_g / 4.0),
                    (limit.mean_relative_velocity_kms, limit.limit_cm2_g),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
        ax.plot([limit.mean_relative_velocity_kms / 1.18, limit.mean_relative_velocity_kms * 1.18],
                [limit.limit_cm2_g] * 2, color=INK, lw=1.2)
    ax.annotate("group and cluster\nupper limits (95%)", (1500.0, 1.1), xytext=(0, 8), textcoords="offset points",
                ha="center", va="bottom", fontsize=8.5, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(mean_v[0], mean_v[-1])
    ax.set_ylim(args.ymin, args.ymax)
    ax.set_xlabel(r"mean relative velocity $\langle v\rangle$ [km/s]")
    ax.set_ylabel(r"$\langle\sigma_V v\rangle/(\langle v\rangle\, m_{\rm DM})$ [cm$^2$/g], before abundance rescaling")
    _style(ax)
    fig.legend(fontsize=8.5, loc="lower left", bbox_to_anchor=(0.06, 0.0), ncol=1)
    fig.tight_layout(rect=(0, 0.034 * (len(selection) + int(args.tuned_example)) + 0.03, 1, 1))
    return _save(fig, args.output_dir, Path(args.input).name.split(".")[0] + "_sidm_points", args.formats)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p):
        p.add_argument("--output-dir", default="plots/sidm")
        p.add_argument("--formats", nargs="+", default=["png", "pdf"])

    landscape = sub.add_parser("landscape", help="sigma_V/m against M2/M3")
    landscape.add_argument("--m3", type=float, default=60.0)
    landscape.add_argument("--k233", type=float, nargs="+", default=[1.0, 10.0, 50.0, 120.0])
    landscape.add_argument("--velocities", type=float, nargs="+", default=[30.0, 1000.0])
    landscape.add_argument("--lx", type=float, default=0.1)
    landscape.add_argument("--gamma2-other", type=float, default=1.0e-7)
    landscape.add_argument("--ratio-min", type=float, default=0.07)
    landscape.add_argument("--ratio-max", type=float, default=3.5)
    landscape.add_argument("--ymin", type=float, default=1.0e-14)
    landscape.add_argument("--ymax", type=float, default=1.0e3)
    common(landscape)
    velocity = sub.add_parser("velocity", help="sigma_V/m against v for M3 > M2 benchmarks")
    velocity.add_argument("--m3", type=float, default=60.0)
    velocity.add_argument("--k233", type=float, default=120.0)
    velocity.add_argument("--vmin", type=float, default=1.0)
    velocity.add_argument("--vmax", type=float, default=1.0e4)
    velocity.add_argument("--n", type=int, default=240)
    common(velocity)
    scan = sub.add_parser("scan", help="scan points against two limits")
    scan.add_argument("input")
    scan.add_argument("--limits", help="limit table (default datafiles/sidm_limits_v1.json)")
    scan.add_argument("--panels", nargs="+", default=["clusters_sagunski2021", "draco_read2018"],
                      help="limit ids shown, one panel each")
    scan.add_argument("--omega-column", default="dm_omega")
    scan.add_argument("--min-fraction", type=float, default=0.1)
    common(scan)
    points = sub.add_parser("points", help="sigma_V/m against velocity for selected scan points")
    points.add_argument("input")
    points.add_argument("--rows", type=int, nargs="+", help="0-based row indices (default: automatic selection)")
    points.add_argument("--omega-column", default="dm_omega")
    points.add_argument("--min-fraction", type=float, default=0.9,
                        help="abundance fraction above which H3 counts as the dark matter")
    points.add_argument("--no-tuned-example", dest="tuned_example", action="store_false",
                        help="omit the tuned near-threshold example")
    points.add_argument("--limits", help="upper-limit table (default datafiles/sidm_limits_v1.json)")
    points.add_argument("--lower-table", help="lower-limit table (default datafiles/sidm_lower_limits_v1.json)")
    points.add_argument("--vmin", type=float, default=10.0)
    points.add_argument("--vmax", type=float, default=5000.0)
    points.add_argument("--n", type=int, default=30)
    points.add_argument("--ymin", type=float, default=1.0e-10)
    points.add_argument("--ymax", type=float, default=1.0e3)
    common(points)
    args = parser.parse_args(argv)
    if args.command == "landscape" and len(args.k233) > len(SERIES):
        parser.error(f"at most {len(SERIES)} values of --k233")
    {"landscape": run_landscape, "velocity": run_velocity, "scan": run_scan,
     "points": run_points}[args.command](args)


if __name__ == "__main__":
    main()
