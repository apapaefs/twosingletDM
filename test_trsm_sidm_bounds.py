"""Tests for trsm_sidm_bounds.py.  Run with ``python -m unittest test_trsm_sidm_bounds.py``."""

import csv
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import trsm_sidm as S
import trsm_sidm_bounds as B


def contact_point(lX=0.5, M3=50.0):
    return S.SIDMPoint(M3=M3, M2=300.0, K133=0.0, K233=0.0, lX=lX)


def resonant_light_point(K233=2.0):
    """M3 = 5 GeV with the h2 pole at v = 1900 km/s: large sigma/m at cluster velocities."""
    M3, beta = 5.0, 0.5 * 1900.0 / S.C_KMS
    return S.SIDMPoint(M3=M3, M2=2.0 * M3 / math.sqrt(1.0 - beta**2), K133=0.0, K233=K233,
                       lX=0.1, gamma2_other=1e-9)


class TableTests(unittest.TestCase):
    def test_default_table_loads_and_sets_are_nested(self):
        default, table = B.load_limits()
        stringent, _ = B.load_limits(selected="stringent")
        everything, _ = B.load_limits(selected="all")
        dwarf, _ = B.load_limits(selected="dwarf")
        ids = lambda limits: {limit.id for limit in limits}
        self.assertEqual(table["schema"], "trsm_sidm_limits_v1")
        self.assertTrue(ids(default) < ids(stringent) < ids(everything))
        self.assertTrue(ids(dwarf) <= ids(everything))
        self.assertFalse(ids(dwarf) & ids(stringent))
        self.assertEqual(len(everything), len(table["limits"]))

    def test_published_numbers(self):
        limits = {limit.id: limit for limit in B.load_limits(selected="all")[0]}
        clusters = limits["clusters_sagunski2021"]
        self.assertEqual((clusters.limit_cm2_g, clusters.mean_relative_velocity_kms, clusters.moment), (0.35, 1900.0, 1))
        self.assertEqual(limits["groups_sagunski2021"].limit_cm2_g, 1.1)
        self.assertEqual(limits["clusters_andrade2022"].limit_cm2_g, 0.13)
        self.assertEqual(limits["bullet_randall2008"].confidence_level, "68%")
        self.assertEqual(limits["mergers_wittman2018"].velocity_range_kms, (500.0, 4000.0))
        for limit in limits.values():
            self.assertTrue(limit.reference.startswith("arXiv:"))

    def test_dispersion_from_mean_relative_velocity(self):
        limit = B.load_limits()[0][0]
        self.assertAlmostEqual(limit.dispersion_kms, limit.mean_relative_velocity_kms * math.sqrt(math.pi) / 4.0)

    def test_bad_tables_are_rejected(self):
        with self.assertRaises(ValueError):
            B.load_limits(selected="no-such-set")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "limits.json"
            path.write_text(json.dumps({"schema": "other", "sets": {}, "limits": []}))
            with self.assertRaises(ValueError):
                B.load_limits(path)

    def test_moment_override(self):
        limits, _ = B.load_limits(moment_override=5)
        self.assertTrue(all(limit.moment == 5 for limit in limits if limit.kind == "maxwell"))


class AssessmentTests(unittest.TestCase):
    def setUp(self):
        self.limits, _ = B.load_limits(selected="all")

    def test_constant_cross_section_gives_the_same_value_for_every_limit(self):
        p = contact_point()
        sigma, _ = S.cross_sections(p, [1e-6])
        expected = float(S.per_mass_cm2_g(sigma[0], p.M3))
        out = B.assess_point(p, 0.12, self.limits, set_name="all")
        for limit in self.limits:
            # constant up to the O(v^2) change of s: 36 lX^2/(32 pi s)
            self.assertAlmostEqual(out[f"sidm_{limit.id}_sigma_over_m"] / expected, 1.0, places=3)
            self.assertAlmostEqual(out[f"sidm_{limit.id}_ratio"] * limit.limit_cm2_g / expected, 1.0, places=3)
        self.assertIs(out["sidm_excluded"], False)
        self.assertEqual(out["sidm_bounds_strongest"], "segue1_hayashi2021")   # smallest limit

    def test_resonant_light_dark_matter_is_excluded_by_cluster_cores(self):
        limits, _ = B.load_limits()
        out = B.assess_point(resonant_light_point(), 0.12, limits)
        self.assertEqual(out["sidm_regime"], "perturbative")
        self.assertGreater(out["sidm_clusters_sagunski2021_ratio"], 1.0)
        self.assertIs(out["sidm_excluded"], True)
        self.assertEqual(out["sidm_bounds_strongest"], "clusters_sagunski2021")
        # the smallest value over 500-4000 km/s is far from the pole
        self.assertLess(out["sidm_mergers_wittman2018_ratio"], 1.0)

    def test_abundance_rescaling(self):
        limits, _ = B.load_limits()
        p = resonant_light_point()
        full = B.assess_point(p, 0.12, limits)
        half = B.assess_point(p, 0.06, limits)                          # default: xi^2
        linear = B.assess_point(p, 0.06, limits, abundance_power=1.0)
        self.assertAlmostEqual(half["sidm_abundance_fraction"], 0.5)
        self.assertEqual(half["sidm_abundance_power"], 2.0)
        self.assertAlmostEqual(half["sidm_bounds_max_ratio"] / full["sidm_bounds_max_ratio"], 0.25)
        self.assertAlmostEqual(linear["sidm_bounds_max_ratio"] / full["sidm_bounds_max_ratio"], 0.5)
        self.assertEqual(half["sidm_clusters_sagunski2021_sigma_over_m"],
                         full["sidm_clusters_sagunski2021_sigma_over_m"])   # model value is not rescaled
        over = B.assess_point(p, 0.5, limits)
        self.assertEqual(over["sidm_abundance_fraction"], 1.0)

    def test_subdominant_component_is_not_constrained(self):
        limits, _ = B.load_limits()
        p = resonant_light_point()
        out = B.assess_point(p, 0.12 * 0.05, limits)
        self.assertIs(out["sidm_excluded"], False)
        self.assertEqual(out["sidm_bounds_status"], "not_applicable")
        self.assertIs(out["sidm_bounds_applicable"], False)
        self.assertGreater(out["sidm_clusters_sagunski2021_sigma_over_m"], 0.35)
        no_rescale = B.assess_point(p, 0.12 * 0.05, limits, rescale=False)
        self.assertIs(no_rescale["sidm_excluded"], True)
        self.assertEqual(no_rescale["sidm_abundance_fraction"], 1.0)

    def test_missing_abundance_is_unavailable_unless_rescaling_is_off(self):
        limits, _ = B.load_limits()
        p = resonant_light_point()
        out = B.assess_point(p, math.nan, limits)
        self.assertTrue(math.isnan(out["sidm_excluded"]))
        self.assertEqual(out["sidm_bounds_status"], "unavailable")
        self.assertIs(B.assess_point(p, math.nan, limits, rescale=False)["sidm_excluded"], True)

    def test_strong_coupling_is_unavailable_when_bounds_could_apply(self):
        limits, _ = B.load_limits()
        strong = S.SIDMPoint(M3=5.0, M2=3.0, K133=0.0, K233=30.0, lX=0.1)
        self.assertEqual(strong.regime(), "strong")
        out = B.assess_point(strong, 0.12, limits)
        self.assertGreater(out["sidm_bounds_unitarity_ratio"], 1.0)
        self.assertTrue(math.isnan(out["sidm_excluded"]))
        self.assertEqual(out["sidm_bounds_status"], "unavailable")
        # but a strongly coupled, tiny component is simply not constrained
        self.assertIs(B.assess_point(strong, 1e-6, limits)["sidm_excluded"], False)

    def test_strong_coupling_out_of_reach_by_unitarity_is_not_excluded(self):
        limits, _ = B.load_limits()
        heavy = S.SIDMPoint(M3=60.0, M2=500.0, K133=0.0, K233=400.0, lX=0.1)
        self.assertEqual(heavy.regime(), "strong")
        out = B.assess_point(heavy, 0.12, limits)
        self.assertLess(out["sidm_bounds_unitarity_ratio"], 1.0)
        self.assertIs(out["sidm_excluded"], False)

    def test_unitarity_values_bound_the_model(self):
        limits, _ = B.load_limits(selected="all")
        p = resonant_light_point(K233=20.0)       # broad resonance near unitarity
        values, _ = B.model_values(p, limits)
        ceilings, short_range = B.unitarity_values(p, limits)
        self.assertTrue(short_range)
        for limit in limits:
            self.assertLess(values[limit.id], ceilings[limit.id])
        # Maxwellian average of 8 pi/k^2 against direct quadrature
        limit = next(l for l in limits if l.id == "clusters_sagunski2021")
        disp = limit.dispersion_kms / S.C_KMS
        v = np.linspace(1e-9, 40.0 * disp, 400001)
        f = v**2 * np.exp(-(v**2) / (4.0 * disp**2))
        direct = np.sum(f * v * 32.0 * math.pi / (p.M3 * v) ** 2) / np.sum(f * v)
        self.assertAlmostEqual(ceilings[limit.id] / (B.UNITARITY_WEIGHT * float(S.per_mass_cm2_g(direct, p.M3))),
                               1.0, places=3)

    def test_fixed_range_uses_the_requested_combination(self):
        limits, _ = B.load_limits()
        p = resonant_light_point()
        merger = next(l for l in limits if l.kind == "fixed")
        _, sigma_v = S.cross_sections(p, merger.velocities_kms / S.C_KMS)
        values, _ = B.model_values(p, limits)
        self.assertAlmostEqual(values[merger.id], float(np.min(S.per_mass_cm2_g(sigma_v, p.M3))))


class ScanTests(unittest.TestCase):
    ROW = {"M2": "300", "M3": "50", "vs": "500", "vx": "0", "a12": "0.001", "lX": "0.5",
           "lPhiX": "0.0", "lSX": "0.0", "K133": "0", "K233": "0", "w1": "0.0041", "w2": "0.1",
           "dm_omega": "0.12"}

    def _run(self, rows, *extra):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "scan.dat", Path(tmp) / "bounds.dat"
            with src.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(self.ROW), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            B.main(["scan", str(src), "--output", str(dst), *extra])
            with dst.open() as stream:
                out = list(csv.DictReader(stream, delimiter="\t"))
            meta = json.loads((Path(tmp) / "bounds.dat.sidm-bounds.json").read_text())
            self.assertEqual(src.read_text().splitlines()[1], "\t".join(rows[0].values()))  # source untouched
        return out, meta

    def test_scan_round_trip(self):
        beta = 0.5 * 1900.0 / S.C_KMS
        M2 = 10.0 / math.sqrt(1.0 - beta**2)
        w2 = S.scalar_to_identical_scalar_width(5.0, M2, 2.0) + 1e-9     # physical total width
        resonant = dict(self.ROW, M3="5", M2=repr(M2), K233="2", w2=repr(w2), lX="0.1")
        rows, meta = self._run([self.ROW, resonant, dict(self.ROW, w1="nan"), dict(resonant, dm_omega="1e-4")])
        self.assertEqual([r["sidm_excluded"] for r in rows], ["False", "True", "nan", "False"])
        self.assertEqual([r["sidm_bounds_status"] for r in rows], ["ok", "ok", "failed", "not_applicable"])
        self.assertEqual(rows[1]["sidm_bounds_strongest"], "clusters_sagunski2021")
        self.assertEqual(rows[0]["dm_omega"], "0.12")
        self.assertEqual(meta["counts"], {"excluded": 1, "not_excluded": 2, "unavailable": 1, "not_applicable": 1})
        self.assertEqual(meta["limit_set"], "default")
        self.assertEqual(len(meta["limits"]), 3)
        self.assertEqual(meta["abundance"]["power"], 2.0)

    def test_missing_omega_column_needs_no_rescale(self):
        row = {k: v for k, v in self.ROW.items() if k != "dm_omega"}
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "scan.dat"
            src.write_text("\t".join(row) + "\n" + "\t".join(row.values()) + "\n")
            with self.assertRaises(SystemExit):
                B.main(["scan", str(src), "--output", str(Path(tmp) / "out.dat")])
            B.main(["scan", str(src), "--output", str(Path(tmp) / "out.dat"), "--no-rescale"])
            self.assertIn("sidm_excluded", (Path(tmp) / "out.dat").read_text().splitlines()[0])


if __name__ == "__main__":
    unittest.main()
