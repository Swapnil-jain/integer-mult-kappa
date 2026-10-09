import json, pathlib, sys, unittest
from fractions import Fraction as Q

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import certificate_round8 as c8


class Round8(unittest.TestCase):
    """the round-eight assembly from the frozen histograms (certificates/round8/)"""
    @classmethod
    def setUpClass(cls):
        cls.bit = c8.load(str(ROOT / 'certificates' / 'round8' / 'bit_hist.json.gz'))
        cls.cx = c8.load(str(ROOT / 'certificates' / 'round8' / 'cx_hist.json.gz'))
        cls.r = c8.certify(cls.bit, cls.cx)

    def test_histograms(self):
        for d in (self.bit, self.cx):
            self.assertEqual(sum(w * n for w, n in d['hist'].items()), d['s'])
            self.assertLess(max(d['hist']), d['m'])
            self.assertGreater(d['W'] * d['m'], d['s'])

    def test_claims_pinned(self):
        self.assertEqual(self.r['astar'], Q(self.bit['claim']['a_star']))
        self.assertEqual(self.r['a_c'], Q(self.cx['claim']['a_c']))

    def test_grid_point_is_largest(self):
        for c, a in ((self.r['cb'], self.r['astar']), (self.r['cc'], self.r['a_c'])):
            lu = {w: c8.ln_up_grid(Q(c['m'], w)) for w in c['hist']}
            self.assertLess(c8.F_exp_upper(c, a, lu), 1)
            self.assertFalse(c8.F_exp_upper(c, a + Q(1, c8.DEN), lu) < 1)

    def test_stopped_mix(self):
        r = self.r
        self.assertEqual(r['a_b'], (1 - r['theta']) * r['astar'] + r['theta'] * r['a_old'])
        self.assertLess(r['a_b'], r['theta'])
        self.assertLess(r['a_b'], r['astar'])

    def test_kappa(self):
        k = self.r['k']
        self.assertTrue(k['ok'], k['bad'])
        self.assertGreater(k['kappa'], c8.TARGET)

    def test_lean_input_matches(self):
        with open(ROOT / 'lean' / 'round8-histograms.json') as f:
            d = json.load(f)
        self.assertEqual(Q(*d['bit']['a']), self.r['astar'])
        self.assertEqual(Q(*d['bit']['stop']['ab']), self.r['a_b'])
        self.assertEqual(Q(*d['cx']['a']), self.r['a_c'])
        self.assertEqual(Q(*d['kappa']['kappa']), self.r['k']['kappa'])

    def test_exp_bound_control(self):
        # a saving well above the certified one must fail the same bound
        c = self.r['cc']; lu = {w: c8.ln_up_grid(Q(c['m'], w)) for w in c['hist']}
        self.assertFalse(c8.F_exp_upper(c, self.r['a_c'] * Q(1001, 1000), lu) < 1)


if __name__ == '__main__':
    unittest.main()
