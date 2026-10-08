import importlib.util, pathlib, sys, unittest
from fractions import Fraction as Q

ROOT = pathlib.Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


cert = load('certificate')
ident = load('check_identities')


class Witness(unittest.TestCase):
    def test_witness_passes(self):
        out = cert.evaluate()
        self.assertEqual(out['kappa'], '249/50000000000')
        self.assertEqual(out['binding'], 'simultaneous butterflies')

    def rejects(self, reason, **kw):
        with self.assertRaises(cert.Failure) as ctx:
            cert.evaluate(**kw)
        self.assertIn(reason, str(ctx.exception))

    def test_kappa_above_margin_rejected(self):
        self.rejects('kappa', kappa=Q(500, 10**11))

    def test_gaussian_row_rejects_eps_two_thirds(self):
        self.rejects('kappa', eps=Q(2, 3))

    def test_precision_rejects_small_x(self):
        self.rejects('precision', x=Q(1, 5), beta=Q(4, 5))

    def test_guard_rejects_old_exponent_without_large_digits(self):
        self.rejects('guard', x=Q(0), eps=Q(666, 1000))

    def test_leaf_rejects_large_beta(self):
        self.rejects('leaf', beta=Q(9, 10))


class Identities(unittest.TestCase):
    def test_small_cases_exact(self):
        for s, t in ((113, 128), (241, 256)):
            out = ident.check(s, t)
            self.assertGreater(out['checks']['segment'], 0)
            self.assertGreater(out['checks']['cross'], 0)
            self.assertGreaterEqual(Q(out['min_cross_exponent']), 2)


if __name__ == '__main__':
    unittest.main()
