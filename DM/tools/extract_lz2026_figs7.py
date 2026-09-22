#!/usr/bin/env python3
"""Reproduce a temporary SI limit from the vector figures of arXiv:2609.02823v1.

Requires pdfplumber only when regenerating the data, not when running scans.
The complete delta curve is an audit artifact; only its delta=0 endpoint is
used for the elastic SI mass table. The mass dependence is an approximation.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path


SOURCE_URL = "https://arxiv.org/src/2609.02823v1"
FIGURE_HASHES = {
    "FigS7_O1_as_SI.pdf": "1f9624781aeb4e4e174fe532786fda0e3648960da7a5427fce956f8530a9fed1",
    "Fig6_O1_L10_limit_stacked.pdf": "034f64d8eb8206d2379e274c7fb12057b47f9e2987e6a37f597cc9a67856aee4",
}


def observed_upper_curve(page, x_left, x_right):
    # Both figures show the same eight splitting hypotheses, 0..350 keV.
    # The expected median is dashed; the solid lower bound starts at delta>0.
    candidates = [
        curve for curve in page.curves
        if curve["stroke"] and not curve["fill"]
        and curve["stroking_color"] == 0 and not curve.get("dash")
        and math.isclose(curve["linewidth"], 3.0)
        and len(curve["pts"]) == 8
        and math.isclose(curve["pts"][0][0], x_left, abs_tol=1e-5)
        and math.isclose(curve["pts"][-1][0], x_right, abs_tol=1e-5)
    ]
    if len(candidates) != 1:
        raise ValueError("Could not uniquely identify the solid observed upper curve")
    points = candidates[0]["pts"]
    for index, (x, _y) in enumerate(points):
        if not math.isclose(350 * (x - x_left) / (x_right - x_left), 50 * index, abs_tol=1e-5):
            raise ValueError("Unexpected mass-splitting grid in figure")
    return points


def extract(source_dir):
    import pdfplumber

    for name, expected in FIGURE_HASHES.items():
        actual = hashlib.sha256((source_dir / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"{name}: source hash differs from arXiv v1")

    with pdfplumber.open(source_dir / "FigS7_O1_as_SI.pdf") as pdf:
        points = observed_upper_curve(pdf.pages[0], 70.7, 491.18)
        # Calibrated against the labelled log-axis endpoints: 10^-35..10^-47.
        cross_sections = [10 ** (-35 - 12 * (y - 14.7) / (330.06 - 14.7)) for _x, y in points]
    with pdfplumber.open(source_dir / "Fig6_O1_L10_limit_stacked.pdf") as pdf:
        coefficient_points = observed_upper_curve(pdf.pages[0], 66.7, 487.18)
        # Major ticks at 10^-2 and 10^-9 in the top panel (not the frame edges).
        coefficients = [
            10 ** (-2 - 7 * (y - 31.155258) / (269.136452 - 31.155258))
            for _x, y in coefficient_points
        ]
    mu_n = 1000 * 0.939 / (1000 + 0.939)
    conversion = mu_n**2 / (math.pi * 246.2**4) * 3.8937966e-28
    converted = [coefficient * conversion for coefficient in coefficients]
    differences = [sigma6 / sigma7 - 1 for sigma6, sigma7 in zip(converted, cross_sections)]
    if max(abs(value) for value in differences) > 0.01:
        raise ValueError("Fig. 6 / Eq. (6) normalization check disagrees with Fig. S7 by over 1%")

    report = {
        "schema": "lz2026_figure_s7_digitization_v1",
        "source": "https://arxiv.org/abs/2609.02823v1",
        "source_archive": SOURCE_URL,
        "source_files_sha256": FIGURE_HASHES,
        "figure": "S7",
        "DM_mass_GeV": 1000.0,
        "operator": "isoscalar O1, scalar normalization",
        "limit_kind": "observed_upper",
        "confidence_level": 0.9,
        "method": "Extract solid upper-curve vector vertices; invert linear delta and log10 sigma axes",
        "warning": "delta_keV is inelastic mass splitting, not DM mass. Only delta=0 is applicable to elastic TRSM.",
        "figure_s7_axes_pdf_points": {
            "x_at_delta_0": 70.7, "x_at_delta_350": 491.18,
            "y_at_sigma_1e_minus35": 14.7, "y_at_sigma_1e_minus47": 330.06,
            "y_direction": "downwards",
        },
        "normalization_check": {
            "source": "Fig. 6 top panel and supplemental Eq. (6)",
            "higgs_vev_GeV": 246.2,
            "nucleon_mass_GeV": 0.939,
            "GeV_minus2_to_cm2": 3.8937966e-28,
            "max_relative_difference": max(abs(value) for value in differences),
            "note": "A numerical cross-check, not an uncertainty estimate for the mass-scaling approximation",
        },
        "points": [
            {
                "delta_keV": float(50 * index),
                "upper_limit_cm2": sigma,
                "figure_s7_vector_xy": list(points[index]),
                "figure_6_vector_xy": list(coefficient_points[index]),
                "figure_6_squared_coefficient": coefficients[index],
                "figure_6_converted_sigma_cm2": converted[index],
            }
            for index, sigma in enumerate(cross_sections)
        ],
    }
    # Three significant figures avoid suggesting official tabulated precision.
    anchor = float(f"{cross_sections[0]:.3g}")
    table = {
        "schema": "trsm_si_upper_limit_v1",
        "label": "lz2026-figs7-highmass-approx",
        "source": "https://arxiv.org/abs/2609.02823v1, Fig. S7 at delta=0; approximate mass scaling",
        "confidence_level": 0.9,
        "interaction": "elastic_isoscalar_si",
        "limit_kind": "observed_upper",
        "cross_section": "per_nucleon",
        "cross_section_unit": "cm2",
        "provenance": {
            "status": "temporary_high_mass_approximation",
            "official_LZ_mass_table": False,
            "published_anchor_mass_GeV": 1000.0,
            "published_anchor_delta_keV": 0.0,
            "extracted_anchor_sigma_cm2": cross_sections[0],
            "rounded_anchor_sigma_cm2": anchor,
            "assumed_mass_dependence": "sigma_limit(m) = sigma_limit(1000 GeV) * m / (1000 GeV)",
            "validity_range_GeV": [400.0, 4000.0],
            "motivation": "Heavy-WIMP fixed-shape rate scales approximately as sigma/m. The paper notes nearly degenerate recoil shapes above 400 GeV in its look-elsewhere calculation.",
            "limitation": "The mass dependence is assumed, not a published LZ limit curve or a likelihood recast. Its uncertainty and coverage away from the 1 TeV anchor have not been validated.",
            "event_treatment": "Published observed upper endpoint; no event removal and no signal/lower-bound requirement",
            "extraction_report": "figure-s7-o1-1tev-digitized.json",
            "source_files_sha256": FIGURE_HASHES,
        },
        "points": [
            {"mass_GeV": mass, "upper_limit": float(f"{anchor * mass / 1000:.8g}")}
            for mass in (400.0, 1000.0, 4000.0)
        ],
    }
    return report, table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path, help="Extracted arXiv v1 source directory")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data/lz2026")
    parser.add_argument("--check", action="store_true", help="Compare with checked-in JSON without modifying it")
    args = parser.parse_args()
    report, table = extract(args.source_dir)
    if not args.check:
        args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, content in (
        ("figure-s7-o1-1tev-digitized.json", report),
        ("lz2026-figs7-highmass-approx.json", table),
    ):
        text = json.dumps(content, indent=2, allow_nan=False) + "\n"
        path = args.output_dir / name
        if args.check:
            if path.read_text(encoding="utf-8") != text:
                raise ValueError(f"Digitization does not reproduce {path}")
        else:
            path.write_text(text, encoding="utf-8")
        print(path)
    print(f"Extracted sigma(1000 GeV) = {report['points'][0]['upper_limit_cm2']:.8g} cm^2")
    print(f"Eq. (6) relative cross-check: {report['normalization_check']['max_relative_difference']:.3%}")


if __name__ == "__main__":
    main()
