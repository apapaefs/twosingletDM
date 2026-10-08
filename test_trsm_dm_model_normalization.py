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
    def test_two_higgs_contacts_are_absent(self):
        # The two-Higgs photon/gluon contacts (legacy LAAh, LGGh*RQCDh) are
        # commented out in TRSM_mixed.mdl: no tree-level 2->2 process used by
        # micrOMEGAs (relic, indirect, direct, widths) contains them.
        pairs = (('h1', 'h1'), ('h1', 'h2'), ('h2', 'h2'))
        for variant in ('h4GOn', 'h4GOff'):
            vertices = calc_hep_vertices(MODEL_ROOT / variant / 'lgrng1.mdl')
            for boson in ('A', 'G'):
                for pair in pairs:
                    with self.subTest(variant=variant, boson=boson, pair=pair):
                        self.assertNotIn((boson, boson, *pair), vertices)
            functions = compact_text(MODEL_ROOT / variant / 'func1.mdl')
            for name in ('LAAh|', 'LGGh|', 'aQCDh|', 'RQCDh|'):
                with self.subTest(variant=variant, parameter=name):
                    self.assertNotIn(name, functions)

    def test_two_higgs_four_gluon_chain_is_absent(self):
        # No vertex may use the x1/X1 or G2t/G2T auxiliary fields of the removed
        # h_i h_j + 4 gluon chain. (Unused auxiliary declarations are harmless.)
        for variant in ('h4GOn', 'h4GOff'):
            vertices = calc_hep_vertices(MODEL_ROOT / variant / 'lgrng1.mdl')
            for particles in vertices:
                with self.subTest(variant=variant, particles=particles):
                    self.assertFalse({'x1', 'X1', 'G2t.t', 'G2T.t'} & set(particles))

    def test_single_higgs_loop_coefficients(self):
        # One mixing projection per mediator, loop evaluated at its own mass.
        for variant in ('h4GOn', 'h4GOff'):
            functions = compact_text(MODEL_ROOT / variant / 'func1.mdl')
            vertices = calc_hep_vertices(MODEL_ROOT / variant / 'lgrng1.mdl')
            for index, projection, mass in ((1, 'CosT', 'Mh'), (2, 'SinT', 'Mh2')):
                for boson, channel, prefix in (('A', 22, 'LAA'), ('G', 21, 'LGG')):
                    coefficient = f'{prefix}{index}'
                    with self.subTest(variant=variant, coefficient=coefficient):
                        self.assertIn(f'{coefficient}|{projection}*trsm_loop_abs({mass},{channel})', functions)
                        self.assertEqual(vertices[(boson, boson, f'h{index}', '')], (
                            '-4*' + coefficient, 'p1.p2*m1.m2-m2.p1*m1.p2',
                        ))

    def test_single_higgs_four_gluon_auxiliary_chains_are_present(self):
        on = calc_hep_vertices(MODEL_ROOT / 'h4GOn/lgrng1.mdl')
        off = calc_hep_vertices(MODEL_ROOT / 'h4GOff/lgrng1.mdl')
        particles = compact_text(MODEL_ROOT / 'h4GOn/prtcls1.mdl')
        self.assertIn('G1T|G1t|G1T|', particles)
        # The two sides join through the constant G1t/G1T tensor propagator.
        # These are h+4g operators even though no five-leg row can appear.
        self.assertEqual(on[('G', 'G', 'G1T.t', '')][0], '-GG/2')
        for higgs, coefficient in (('h1', 'LGG1'), ('h2', 'LGG2')):
            self.assertEqual(on[('G', 'G', 'G1t.t', higgs)][0], '-GG*' + coefficient + '/2')
            self.assertEqual(on[('G', 'G', 'G', higgs)][0], '-4*GG*' + coefficient)
            self.assertNotIn(('G', 'G', 'G1t.t', higgs), off)
            self.assertEqual(on[('G', 'G', higgs, '')][0], '-4*' + coefficient)

    def test_lanhep_source_uses_canonical_x_normalization(self):
        # The X sector is written with Xi = X/sqrt(2). This equals the canonical
        # form (dX)^2/2 - muX*X^2/2 - LX/4*X^4 - LSX/2*X^2*Si^2 - LHX/2*X^2*|Phi|^2,
        # so MX^2 = muX + (LHX*v^2 + LSX*vs^2)/2 and both portal terms enter with +.
        model = compact_text(MODEL_ROOT / "lanhep_mdl" / "TRSM_mixed.mdl")

        self.assertIn("muX=MX**2-(LHX*(2*MW*SW/EE)**2+LSX*vevs**2)/2", model)
        self.assertNotIn("(LHX*(2*MW*SW/EE)**2-LSX*vevs**2)", model)
        self.assertIn("letXi='~X'/Sqrt2.", model)
        self.assertIn("lterm(deriv^mu*Xi)*(deriv^mu*Xi)-muX*(Xi**2)-LX*(Xi**4).", model)
        self.assertIn("lterm-LSX*(Xi**2)*Si**2-LHX*(Xi**2)*shd*shD.", model)

    def test_generated_mass_relation_has_half_portal_terms(self):
        expected = "muX|MX^2-(LHX*(2*MW*SW/EE)^2+LSX*vevs^2)/2"

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
                    import sympy as sp
                    functions={}
                    for line in (MODEL_ROOT/variant/'func1.mdl').read_text().splitlines():
                        parts=line.split('|')
                        if len(parts)>1 and parts[0].strip().startswith('B'):
                            functions[sp.Symbol(parts[0].strip())]=sp.sympify(parts[1].strip().replace('^','**'))
                    actual=sp.sympify('*'.join('('+v+')' for v in vertices[particles]).replace('^','**'))
                    for _ in range(10):actual=actual.subs(functions)
                    expected_expr=sp.sympify('*'.join('('+v+')' for v in coupling).replace('^','**'))
                    self.assertEqual(sp.simplify(actual-expected_expr),0)


if __name__ == "__main__":
    unittest.main()
