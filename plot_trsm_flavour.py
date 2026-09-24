"""Flavour observables and renderers used by the comprehensive plot suite."""

import hashlib
import json
import math

import numpy as np
from matplotlib.colors import SymLogNorm

from flavour import flavour_observables as ups

# Colours follow the suite's existing colour-blind-friendly palette.
STYLES = {
    "outside_coverage": ("Outside coverage", "#BDBDBD", "."),
    "passed": ("Passing", "#009E73", "o"),
    "zero_signal": ("Zero signal", "#56B4E9", "d"),
    "mumu": (r"$\mu\mu$ excluded", "#E69F00", "^"),
    "tautau": (r"$\tau\tau$ excluded", "#D55E00", "v"),
    "both": ("Both channels excluded", "#CC79A7", "X"),
    "excluded": ("Excluded; channel unavailable", "#882255", "P"),
    "unassessed": ("Unassessed", "#333333", "x"),
}


def derive_flavour(strings, passed, available):
    size = len(passed)
    categories = np.full(size, "unassessed", dtype=object)
    arrays = {f"flavour_h2_{name}": np.full(size, np.nan) for name in
              ("k2sq", "mumu_prediction", "tautau_prediction")}
    for i in range(size):
        try:
            details = json.loads(strings["flavour_details"][i])
            if not isinstance(details, dict):
                details = {}
        except (ValueError, TypeError):
            details = {}
        channels_excluded = set()
        for scalar in details.values():
            if not isinstance(scalar, dict):
                continue
            for channel, diagnostic in scalar.get("channels", {}).items():
                if diagnostic.get("excluded") is True:
                    channels_excluded.add(channel)
        if available[i]:
            status = strings["flavour_status"][i]
            if passed[i]:
                categories[i] = status if status in ("outside_coverage", "zero_signal") else "passed"
            else:
                categories[i] = ("both" if {"mumu", "tautau"} <= channels_excluded
                                 else next(iter(channels_excluded)) if channels_excluded else "excluded")
        h2 = details.get("h2", {})
        if not isinstance(h2, dict):
            continue
        values = {"k2sq": h2.get("mixing_squared")}
        for channel in ("mumu", "tautau"):
            values[channel + "_prediction"] = h2.get("channels", {}).get(channel, {}).get("prediction")
        for key, value in values.items():
            if isinstance(value, (int, float)) and math.isfinite(value) and value >= 0:
                arrays[f"flavour_h2_{key}"][i] = value
    arrays["flavour_category"] = categories
    arrays["flavour_available"] = available
    return arrays


def _scatter_status(ax, data, x, y, mask):
    for category, (label, colour, marker) in STYLES.items():
        selected = mask & (data.s("flavour_category") == category)
        if np.any(selected):
            ax.scatter(x[selected], y[selected], c=colour, marker=marker, s=25,
                       alpha=.8, linewidths=.7, rasterized=True,
                       label=f"{label}: {np.count_nonzero(selected):,}")


def _impact_text(data):
    before = data.b("pre_flavour_full_viability")
    removed = before & data.b("flavour_available") & ~data.b("flavour")
    unknown = before & ~data.b("flavour_available")
    return (f"Viable before flavour: {np.count_nonzero(before):,}; "
            f"removed: {np.count_nonzero(removed):,}; "
            f"unassessed: {np.count_nonzero(unknown):,}")


def _curve_provenance_matches(data):
    metadata = data.metadata or {}
    config = metadata.get("flavour") or metadata.get("physics_manifest", {}).get("flavour", {})
    hashes = config.get("sources_sha256", {})
    for filename in ups.LIMIT_FILES:
        expected = hashes.get("flavour/" + filename)
        if expected and hashlib.sha256((ups.LIMIT_DATA_DIR / filename).read_bytes()).hexdigest() != expected:
            return False
    return True


def render_flavour(fig, ax, data, spec, style_mass_axis, unavailable):
    x, y = data.f("M2"), data.f("M3")
    valid = np.isfinite(x) & np.isfinite(y)
    if spec.kind in ("flavour_status", "flavour_zoom"):
        if spec.kind == "flavour_zoom":
            valid &= (x >= 4) & (x <= ups.SEARCH_MAX)
        _scatter_status(ax, data, x, y, valid)
        style_mass_axis(ax, data)
        if spec.kind == "flavour_zoom":
            ax.set_xscale("linear")
            ax.set_xlim(4, ups.SEARCH_MAX)
            positive_y = y[valid & (y > 0)]
            if len(positive_y):
                # Fit the zoom to its own points, retaining the invisible-decay
                # threshold instead of inheriting the full scan's M3 ceiling.
                lower = min(1.8, float(positive_y.min()) * .85)
                upper = max(5., float(positive_y.max()) * 1.15)
                ax.set_yscale("log" if upper / lower > 10 else "linear")
                ax.set_ylim(lower, upper)
        ax.axvline(ups.SEARCH_MAX, color="#666666", ls=":", lw=1)
        ax.set_title(spec.title + "\n" + _impact_text(data), fontsize=10)
        if np.any(valid):
            ax.legend(fontsize=8, loc="best")
        else:
            ax.text(.5, .5, "No stored points in this mass window", transform=ax.transAxes, ha="center")
    elif spec.kind == "flavour_strength":
        ratios = data.f("flavour_max_ratio")
        selected = valid & np.isfinite(ratios) & (ratios >= 0)
        if not np.any(selected):
            raise unavailable("no finite flavour prediction/limit ratios")
        maximum = max(10.0, float(np.max(ratios[selected])))
        norm = SymLogNorm(linthresh=1.0, vmin=0, vmax=maximum, base=10)
        ax.scatter(x[valid & ~selected], y[valid & ~selected], c=".85", s=8,
                   label="Ratio unavailable / outside coverage")
        points = ax.scatter(x[selected], y[selected], c=ratios[selected], norm=norm,
                            cmap="viridis", s=22, rasterized=True)
        failed = selected & (ratios > 1)
        ax.scatter(x[failed], y[failed], s=42, facecolors="none", edgecolors="#D55E00",
                   linewidths=.7, label="Excluded: ratio > 1")
        fig.colorbar(points, ax=ax, label="Maximum prediction / individual limit (excluded > 1)")
        style_mass_axis(ax, data)
        ax.set_title(spec.title + "\nLinear colour scale below one; logarithmic above", fontsize=10)
        ax.legend(fontsize=8)
    elif spec.kind == "flavour_products":
        if not any(np.any(np.isfinite(data.f(f"flavour_h2_{c}_prediction"))) for c in ("mumu", "tautau")):
            raise unavailable("no stored h2 lepton-product predictions")
        if not _curve_provenance_matches(data):
            raise unavailable("local flavour curves differ from saved assessment provenance")
        ax.remove()
        fig.set_size_inches(12, 4.8)
        axes = fig.subplots(1, 2)
        for axis, channel, label, filenames in zip(axes, ("mumu", "tautau"),
                (r"$\mu^+\mu^-$", r"$\tau^+\tau^-$"),
                (ups.LIMIT_FILES[:2], ups.LIMIT_FILES[2:])):
            prediction = data.f(f"flavour_h2_{channel}_prediction")
            finite = np.isfinite(prediction) & np.isfinite(x) & (x >= 4) & (x <= ups.SEARCH_MAX)
            positive = finite & (prediction > 0)
            axis.scatter(x[positive], prediction[positive], s=16, c="#0072B2", alpha=.6,
                         rasterized=True, label="TRSM prediction")
            for filename, experiment, colour in zip(filenames, ("Belle", "BaBar"), ("#009E73", "#D55E00")):
                masses, limits = ups._load_limit_curve(filename)
                keep = (masses >= 4) & (masses <= ups.SEARCH_MAX)
                axis.plot(masses[keep], limits[keep], lw=1, color=colour, label=experiment + " 90% limit")
            axis.set(xlim=(4, ups.SEARCH_MAX), yscale="log", xlabel=r"$M_2$ [GeV]",
                     ylabel=r"BR($\Upsilon\to\gamma h_2$) BR($h_2\to\ell\ell$)", title=label)
            axis.set_title(label + f"\nZero predictions omitted: {np.count_nonzero(finite & (prediction == 0)):,}", fontsize=10)
            axis.grid(alpha=.2)
            axis.legend(fontsize=8)
        fig.suptitle("Upsilon lepton products; prompt, narrow-resonance prescription")
    elif spec.kind == "flavour_mixing":
        k2sq = data.f("flavour_h2_k2sq")
        selected = np.isfinite(x) & np.isfinite(k2sq) & (k2sq > 0)
        if not np.any(selected):
            raise unavailable("no positive assessed h2 mixing values")
        _scatter_status(ax, data, x, k2sq, selected)
        opened = selected & (x > 2 * y)
        if np.any(opened):
            ax.scatter(x[opened], k2sq[opened], facecolors="none", edgecolors="black", s=60,
                       lw=.7, label=r"$h_2\to h_3h_3$ kinematically open")
        ax.set(xlim=(4, ups.SEARCH_MAX), yscale="log", xlabel=r"$M_2$ [GeV]",
               ylabel=r"$|k_2|^2$", title=spec.title)
        ax.set_title(spec.title + f"\nZero mixing omitted: {np.count_nonzero(k2sq == 0):,}", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
