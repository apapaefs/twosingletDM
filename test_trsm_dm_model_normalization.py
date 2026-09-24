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
    def test_restored_two_higgs_contacts_match_historical_tables(self):
        # Reference vertices from before 34beddb. Keep both the Lorentz
        # structure and the mixed/identical-Higgs combinatorial factors.
        pairs = {
            ('h1', 'h1'): '-2*CosT^2*EE*{coefficient}/(MW*SW)',
            ('h1', 'h2'): '-EE*{coefficient}*SinTT/(MW*SW)',
            ('h2', 'h2'): '-2*EE*{coefficient}*SinT^2/(MW*SW)',
        }
        for variant in ('h4GOn', 'h4GOff'):
            vertices = calc_hep_vertices(MODEL_ROOT / variant / 'lgrng1.mdl')
            for boson, coefficient in (('A', 'LAAh'), ('G', 'LGGh*RQCDh')):
                for pair, factor in pairs.items():
                    with self.subTest(variant=variant, boson=boson, pair=pair):
                        self.assertEqual(vertices[(boson, boson, *pair)], (
                            factor.format(coefficient=coefficient),
                            'p1.p2*m1.m2-m2.p1*m1.p2',
                        ))

    def test_restored_four_gluon_chain_matches_historical_tables(self):
        on = calc_hep_vertices(MODEL_ROOT / 'h4GOn/lgrng1.mdl')
        off = calc_hep_vertices(MODEL_ROOT / 'h4GOff/lgrng1.mdl')
        expected = {
            ('h1', 'h1', 'x1', ''):
                ('CosT^2*EE*LGGh*Maux*RQCDh/(2*MW*SW)', '1'),
            ('h1', 'h2', 'x1', ''):
                ('EE*LGGh*Maux*RQCDh*SinTT/(4*MW*SW)', '1'),
            ('h2', 'h2', 'x1', ''):
                ('EE*LGGh*Maux*RQCDh*SinT^2/(2*MW*SW)', '1'),
            ('G', 'G', 'G2t.t', 'X1'):
                ('-GG*Maux/2', 'm1.m3*m2.M3-m1.M3*m2.m3'),
            ('G', 'G', 'G2T.t', ''):
                ('-GG/2', 'm1.m3*m2.M3-m1.M3*m2.m3'),
        }
        for particles, vertex in expected.items():
            with self.subTest(particles=particles):
                self.assertEqual(on[particles], vertex)
                self.assertNotIn(particles, off)
        on_particles = compact_text(MODEL_ROOT / 'h4GOn/prtcls1.mdl')
        off_particles = compact_text(MODEL_ROOT / 'h4GOff/prtcls1.mdl')
        # No physical mass pole or decay width for either auxiliary pair.
        for declaration in ('G2T|G2t|G2T||2|Maux|0|8|!*|',
                            'x1|x1|X1||0|Maux|0|1|!*|'):
            self.assertIn(declaration, on_particles)
            self.assertNotIn(declaration, off_particles)

    def test_legacy_coefficients_are_separate_from_single_higgs_loops(self):
        for variant in ('h4GOn', 'h4GOff'):
            functions = compact_text(MODEL_ROOT / variant / 'func1.mdl')
            # Preserve the historical signs, h1 scale and QCD prescription
            # only for the restored contacts; no new matching is implied.
            for expression in (
                'LAAh|-cabs(lAAhiggs(Mh,"h1"))',
                'LGGh|-cabs(lGGhiggs(Mh,"h1"))',
                'aQCDh|alphaQCD(Mh)/acos(-1)',
                'RQCDh|sqrt(1+149/12*aQCDh+68.6482*aQCDh^2-212.447*aQCDh^3)',
            ):
                self.assertIn(expression, functions)
            vertices = calc_hep_vertices(MODEL_ROOT / variant / 'lgrng1.mdl')
            for index, projection, mass in ((1, 'CosT', 'Mh'), (2, 'SinT', 'Mh2')):
                for boson, channel, prefix in (('A', 22, 'LAA'), ('G', 21, 'LGG')):
                    coefficient = f'{prefix}{index}'
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
