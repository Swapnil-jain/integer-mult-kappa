import json, pathlib, sys, unittest
from fractions import Fraction as Q

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'tests'))
import certificate_round11 as c11
import test_round10 as t10     # imported as a module, so its own test classes are not collected again here

FROZEN = ROOT / 'certificates' / 'round11'


class Regressions(unittest.TestCase):
    def test_pr144_round9_round10(self):
        self.assertTrue(all(c11.regressions(verbose=False).values()))


@unittest.skipUnless((FROZEN / 'bit_inventory.json.gz').exists(), 'round eleven not frozen yet')
class Round11(t10.Round10):
    """round ten's frozen-assembly tests (claims, grid points, kappa, Lean input, controls) on round eleven's files"""
    frozen, lean_name = FROZEN, 'round11-histograms.json'

    def test_frozen_equals_release_inventories(self):
        # the frozen inventories are the gated release inventories plus the pinned claims, nothing else
        for frozen, release in ((self.bit, 'bit_inventory_p12.json.gz'), (self.cx, 'cx_recycle_inventory_p11.json.gz')):
            raw = c11.load(str(FROZEN / release))
            self.assertEqual({k: v for k, v in frozen.items() if k != 'claim'}, raw)

    def test_summary(self):
        with open(FROZEN / 'summary.json') as f:
            d = json.load(f)
        self.assertEqual(Q(d['kappa']), self.r['k']['kappa'])
        self.assertEqual(d['bit_binds'], self.r['bit_binds'])
        self.assertTrue(d['inherited_reserve'])


if __name__ == '__main__':
    unittest.main()
