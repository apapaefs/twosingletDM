"""Tests for trsm_sidm_lower_limit.py.  Run with ``python -m unittest test_trsm_sidm_lower_limit.py``."""

import csv
import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import trsm_sidm as S
import trsm_sidm_bounds as B
import trsm_sidm_lower_limit as L


def contact_point(lX=0.5, M3=50.0):
    return S.SIDMPoint(M3=M3, M2=300.0, K133=0.0, K233=0.0, lX=lX)


def tuned_point(M3=20.0, K233=5.0, binding_keV=0.8):
    """h2 pole just below the H3 H3 threshold: large, velocity-dependent sigma/m."""
    return S.SIDMPoint(M3=M3, M2=2.0 * M3 - binding_keV * 1.0e-6, K133=0.0, K233=K233, lX=0.1)


class TableTests(unittest.TestCase):
    def test_table_loads(self):
        default, table = L.load_requirements()
        everything, _ = L.load_requirements(selected="all")
        self.assertEqual(table["schema"], "trsm_sidm_lower_limits_v1")
        self.assertEqual([r.id for r in default], ["cores_tulin_yu2018"])
        self.assertTrue({r.id for r in default} < {r.id for r in everything})
        self.assertEqual(len(everything), len(table["requirements"]))

    def test_published_numbers(self):
        reqs = {r.id: r for r in L.load_requirements(selected="all")[0]}
        cores = reqs["cores_tulin_yu2018"]
        self.assertEqual((cores.minimum_cm2_g, cores.velocity_range_kms, cores.moment), (1.0, (30.0, 200.0), 1))
        self.assertEqual(reqs["tbtf_milky_way_tulin_yu2018"].minimum_cm2_g, 0.6)
        self.assertEqual(reqs["tbtf_local_group_tulin_yu2018"].velocity_range_kms, (50.0, 50.0))
        ren = reqs["rotation_curves_ren2019"]
        self.assertTrue(ren.velocity_assumed)
        # <v> = (4/sqrt(pi)) * 0.64 * V_f for V_f = 20 and 300 km/s
        factor = 4.0 / math.sqrt(math.pi) * 0.64
        self.assertAlmostEqual(ren.velocity_range_kms[0] / (factor * 20.0), 1.0, places=1)
        self.assertAlmostEqual(ren.velocity_range_kms[1] / (factor * 300.0), 1.0, places=1)

    def test_velocity_grid(self):
        cores = L.load_requirements()[0][0]
        grid = cores.mean_velocities_kms
        self.assertEqual((grid[0], len(grid)), (30.0, L.GRID_POINTS))
        self.assertAlmostEqual(grid[-1], 200.0)
        self.assertAlmostEqual(cores.dispersions_kms[0], 30.0 * math.sqrt(math.pi) / 4.0)

    def test_bad_tables_are_rejected(self):
        with self.assertRaises(ValueError):
            L.load_requirements(selected="no-such-set")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.json"
            path.write_text(json.dumps({"schema": "other", "sets": {}, "requirements": []}))
            with self.assertRaises(ValueError):
                L.load_requirements(path)


class AssessmentTests(unittest.TestCase):
    def setUp(self):
        self.default, _ = L.load_requirements()
        self.everything, _ = L.load_requirements(selected="all")

    def test_constant_cross_section(self):
        p = contact_point()
        sigma, _ = S.cross_sections(p, [1e-6])
        expected = float(S.per_mass_cm2_g(sigma[0], p.M3))
        out = L.assess_point(p, 0.12, self.everything, set_name="all")
        for r in self.everything:
            self.assertAlmostEqual(out[f"sidm_lower_{r.id}_min"] / expected, 1.0, places=5)
            self.assertAlmostEqual(out[f"sidm_lower_{r.id}_max"] / expected, 1.0, places=5)
            self.assertAlmostEqual(out[f"sidm_lower_{r.id}_ratio"] * r.minimum_cm2_g / expected, 1.0, places=5)
        self.assertIs(out["sidm_lower_satisfied"], False)
        # the two entries that ask for 1 cm^2/g are the hardest to reach
        self.assertIn(out["sidm_lower_weakest"], ("cores_tulin_yu2018", "rotation_curves_ren2019"))
        default = L.assess_point(p, 0.12, self.default)
        self.assertEqual(default["sidm_lower_reason"], "below cores_tulin_yu2018")
        self.assertAlmostEqual(default["sidm_lower_min_ratio"] / expected, 1.0, places=5)

    def test_tuned_pole_satisfies_the_lower_limit_and_the_upper_limits(self):
        p = tuned_point()
        out = L.assess_point(p, 0.12, self.default)
        self.assertEqual(out["sidm_regime"], "perturbative")
        self.assertIs(out["sidm_lower_satisfied"], True)
        self.assertGreater(out["sidm_lower_cores_tulin_yu2018_min"], 1.0)
        # the cross section falls with velocity: larger at 30 than at 200 km/s
        self.assertGreater(out["sidm_lower_cores_tulin_yu2018_max"], 2.0 * out["sidm_lower_cores_tulin_yu2018_min"])
        upper, _ = B.load_limits(selected="stringent")
        self.assertIs(B.assess_point(p, 0.12, upper)["sidm_excluded"], False)
        # the wider SPARC range reaches velocities where it has dropped below 1
        wide = L.assess_point(p, 0.12, self.everything, set_name="all")
        self.assertIs(wide["sidm_lower_satisfied"], False)
        self.assertEqual(wide["sidm_lower_weakest"], "rotation_curves_ren2019")

    def test_minimum_is_taken_over_the_velocity_range(self):
        p = tuned_point()
        cores = self.default[0]
        values, _ = L.model_values(p, self.default)
        averages = S.maxwell_averages(p, cores.dispersions_kms, [1])
        series = S.per_mass_cm2_g(np.array([averages[(d, 1)] for d in cores.dispersions_kms]), p.M3)
        self.assertAlmostEqual(values[cores.id][0], float(series.min()))
        self.assertAlmostEqual(values[cores.id][1], float(series.max()))
        self.assertEqual(int(np.argmin(series)), len(series) - 1)             # smallest at 200 km/s

    def test_abundance(self):
        p = tuned_point()
        full = L.assess_point(p, 0.12, self.default)
        nearly = L.assess_point(p, 0.95 * 0.12, self.default)
        self.assertAlmostEqual(nearly["sidm_lower_min_ratio"] / full["sidm_lower_min_ratio"], 0.95**2)
        self.assertEqual(nearly["sidm_lower_abundance_power"], 2.0)
        self.assertIs(nearly["sidm_lower_satisfied"], True)
        half = L.assess_point(p, 0.06, self.default)
        self.assertIs(half["sidm_lower_satisfied"], False)
        self.assertEqual(half["sidm_lower_status"], "subdominant")
        self.assertIs(L.assess_point(p, 0.06, self.default, rescale=False)["sidm_lower_satisfied"], True)
        missing = L.assess_point(p, math.nan, self.default)
        self.assertTrue(math.isnan(missing["sidm_lower_satisfied"]))
        self.assertEqual(missing["sidm_lower_status"], "unavailable")

    def test_strong_coupling(self):
        light = S.SIDMPoint(M3=5.0, M2=3.0, K133=0.0, K233=30.0, lX=0.1)
        self.assertEqual(light.regime(), "strong")
        out = L.assess_point(light, 0.12, self.default)
        self.assertGreater(out["sidm_lower_unitarity_ratio"], 1.0)
        self.assertTrue(math.isnan(out["sidm_lower_satisfied"]))
        # a heavy strongly coupled particle cannot reach 1 cm^2/g by unitarity
        heavy = S.SIDMPoint(M3=400.0, M2=900.0, K133=0.0, K233=2000.0, lX=0.1)
        self.assertEqual(heavy.regime(), "strong")
        out = L.assess_point(heavy, 0.12, self.default)
        self.assertLess(out["sidm_lower_unitarity_ratio"], 1.0)
        self.assertIs(out["sidm_lower_satisfied"], False)

    def test_unitarity_values_bound_the_model(self):
        p = tuned_point(K233=50.0, binding_keV=1e-6)      # essentially at unitarity
        values, _ = L.model_values(p, self.everything)
        ceilings, short_range = L.unitarity_values(p, self.everything)
        self.assertTrue(short_range)
        for r in self.everything:
            self.assertLess(values[r.id][0], ceilings[r.id])
            self.assertGreater(values[r.id][0], 0.9 * ceilings[r.id] / B.UNITARITY_WEIGHT)   # s-wave saturated


class ScanTests(unittest.TestCase):
    ROW = {"M2": "300", "M3": "50", "vs": "500", "vx": "0", "a12": "0.001", "lX": "0.5",
           "lPhiX": "0.0", "lSX": "0.0", "K133": "0", "K233": "0", "w1": "0.0041", "w2": "0.1",
           "dm_omega": "0.12"}

    def test_scan_round_trip(self):
        tuned = dict(self.ROW, M3="20", M2=repr(40.0 - 0.8e-6), K233="5", w2="1e-7", lX="0.1")
        rows = [self.ROW, tuned, dict(self.ROW, w1="nan"), dict(tuned, dm_omega="0.05")]
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "scan.dat", Path(tmp) / "lower.dat"
            with src.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(self.ROW), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            L.main(["scan", str(src), "--output", str(dst)])
            with dst.open() as stream:
                out = list(csv.DictReader(stream, delimiter="\t"))
            meta = json.loads((Path(tmp) / "lower.dat.sidm-lower.json").read_text())
        self.assertEqual([r["sidm_lower_satisfied"] for r in out], ["False", "True", "nan", "False"])
        self.assertEqual([r["sidm_lower_status"] for r in out], ["ok", "ok", "failed", "subdominant"])
        self.assertEqual(out[0]["dm_omega"], "0.12")
        self.assertEqual(meta["counts"], {"satisfied": 1, "not_satisfied": 2, "unavailable": 1, "subdominant": 1})
        self.assertEqual(meta["abundance"]["min_fraction"], 0.9)
        self.assertEqual([r["id"] for r in meta["requirements"]], ["cores_tulin_yu2018"])


if __name__ == "__main__":
    unittest.main()
