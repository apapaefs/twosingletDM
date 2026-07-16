import unittest
from pathlib import Path


MODEL_ROOT = Path(__file__).resolve().parent / "DM" / "models"


def compact_text(path):
    return "".join(path.read_text(encoding="ascii").split())


def calc_hep_vertices(path):
    vertices = {}
    for line in path.read_text(encoding="ascii").splitlines():
        fields = tuple(field.strip() for field in line.split("|"))
        if len(fields) != 6:
            continue
        key = fields[:4]
        if key in vertices:
            raise AssertionError(f"duplicate CalcHEP vertex {key!r} in {path}")
        vertices[key] = fields[4:]
    return vertices


class TestTrsmDmModelNormalization(unittest.TestCase):
    def test_lanhep_source_uses_canonical_x_normalization(self):
        model = compact_text(MODEL_ROOT / "lanhep_mdl" / "TRSM_mixed.mdl")

        self.assertIn(
            "muX=(MX**2-LHX*(2*MW/EE*SW)**2/2-LSX*vevs**2/2)/2",
            model,
        )
        self.assertIn("-LX/4*('~X'**4)", model)
        self.assertIn(
            "lterm-LSX/2*('~X'**2)*Si**2-LHX/2*('~X'**2)*shd*shD.",
            model,
        )

    def test_generated_mass_relation_has_half_portal_terms(self):
        expected = "muX|(MX^2-LHX*(2*MW/EE*SW)^2/2-LSX*vevs^2/2)/2"

        for variant in ("h4GOn", "h4GOff"):
            with self.subTest(variant=variant):
                functions = compact_text(MODEL_ROOT / variant / "func1.mdl")
                self.assertIn(expected, functions)

    def test_generated_vertices_match_canonical_normalization(self):
        expected = {
            ("h1", "~X", "~X", ""): (
                "1/EE",
                "EE*LSX*SinT*vevs-2*CosT*LHX*MW*SW",
            ),
            ("h2", "~X", "~X", ""): (
                "-1/EE",
                "CosT*EE*LSX*vevs+2*LHX*MW*SW*SinT",
            ),
            ("h1", "h1", "~X", "~X"): (
                "-1",
                "LSX*SinT^2+CosT^2*LHX",
            ),
            ("h1", "h2", "~X", "~X"): (
                "SinTT/2",
                "LSX-LHX",
            ),
            ("h2", "h2", "~X", "~X"): (
                "-1",
                "CosT^2*LSX+LHX*SinT^2",
            ),
            ("~X", "~X", "~X", "~X"): ("-6*LX", "1"),
        }

        for variant in ("h4GOn", "h4GOff"):
            with self.subTest(variant=variant):
                vertices = calc_hep_vertices(MODEL_ROOT / variant / "lgrng1.mdl")
                for particles, coupling in expected.items():
                    self.assertEqual(vertices[particles], coupling)


if __name__ == "__main__":
    unittest.main()
