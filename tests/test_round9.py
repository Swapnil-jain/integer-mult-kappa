import json, pathlib, sys, unittest
from fractions import Fraction as Q

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import certificate_round9 as c9


class PR144(unittest.TestCase):
    """PR #144's published values, from its published inventories"""
    def test_reproduction(self):
        r = c9.reproduce_pr144(verbose=False)
        self.assertEqual(r['k']['kappa'], Q(4609169, 10**10))
        self.assertEqual(len(r['k']['strict']), 47)
        self.assertEqual(len(r['k']['margins']), 7)


class Round9(unittest.TestCase):
    """the round-nine assembly from the frozen inventories (certificates/round9/)"""
    @classmethod
    def setUpClass(cls):
        cls.bit = c9.load(str(ROOT / 'certificates' / 'round9' / 'bit_inventory.json.gz'))
        cls.cx = c9.load(str(ROOT / 'certificates' / 'round9' / 'cx_inventory.json.gz'))
        cls.r = c9.certify(cls.bit, cls.cx)

    def test_claims_pinned(self):
        self.assertEqual(self.r['bit']['astar'], Q(self.bit['claim']['astar']))
        self.assertEqual(self.r['bit']['a_bit'], Q(self.bit['claim']['a_bit']))
        self.assertEqual(self.r['cx']['a_c'], Q(self.cx['claim']['a_c']))

    def test_grid_point_is_largest(self):
        b = self.r['bit']
        self.assertGreater(b['gap'], 0)
        self.assertLessEqual(b['gap_next'], 0)

    def test_stopped_mix(self):
        b = self.r['bit']
        self.assertEqual(b['a_bit'], (1 - b['theta']) * b['astar'] + b['theta'] * b['a_old'])
        self.assertLess(b['a_bit'], b['theta'])

    def test_kappa(self):
        k = self.r['k']
        self.assertTrue(k['ok'], k['bad'])
        self.assertTrue(self.r['bit_binds'])
        self.assertLess(k['kappa'], k['minimum'])
        self.assertLess(k['minimum'] - k['kappa'], Q(1, c9.KGRID))
        self.assertGreater(k['kappa'], Q(1, 2**12))

    def test_lean_input_matches(self):
        with open(ROOT / 'lean' / 'round9-histograms.json') as f:
            d = json.load(f)
        self.assertEqual(Q(*d['bit']['a']), self.r['bit']['astar'])
        self.assertEqual(Q(*d['bit']['stop']['ab']), self.r['bit']['a_bit'])
        self.assertEqual(Q(*d['cx']['a']), self.r['cx']['a_c'])
        self.assertEqual(Q(*d['kappa']), self.r['k']['kappa'])
        self.assertEqual([tuple(x) for x in d['bit']['hist']], sorted(self.r['bit']['p']['hist'].items()))

    def test_controls(self):
        # one grid point above a*, a kappa at the minimum margin, and a broken child histogram must all fail
        bit = dict(self.bit, cert=dict(self.bit['cert'], saving=str(self.r['bit']['astar'] + Q(1, 10**12))))
        with self.assertRaises(AssertionError):
            c9.certify(bit, self.cx)
        k = c9.assembly(self.r['k']['a'], self.r['cx']['a_c'], self.r['bridge'], kappa=self.r['k']['minimum'])
        self.assertFalse(k['ok'])
        h = dict(self.bit['child_histogram']); h['2'] += 1
        with self.assertRaises(AssertionError):
            c9.certify(dict(self.bit, child_histogram=h), self.cx)


if __name__ == '__main__':
    unittest.main()
