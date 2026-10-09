import json, pathlib, sys, unittest
from fractions import Fraction as Q

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'independent' / 'joint-frame-stack'))
import certificate_round10 as c10
import jfdata

FROZEN = ROOT / 'certificates' / 'round10'


def pair_rule(ops, copies, s0):
    """the old rule (a kept copy is any gate on its (copy, source) pair), kept here only for comparison"""
    kept = {c: x for x, cs in copies.items() for c in cs}
    return {i for i, op in enumerate(ops) if op[0] == 'add' and op[1] in kept and op[2] == s0[kept[op[1]]]}


class KeptCopies(unittest.TestCase):
    """kept V copies are identified by op index (the copy role's first op), never by role pair"""
    def test_truth_table(self):
        copies, s0 = {0: [1]}, {0: 0}
        cases = [
            # (ops, expected op indices, or None when the word must be rejected)
            ([('src', 0, 1), ('add', 1, 0)], {1}),
            ([('src', 0, 1), ('add', 1, 0), ('add', 2, 1), ('add', 1, 0)], {1}),   # the pair is gated again later
            ([('src', 0, 1), ('add', 2, 1), ('add', 1, 0)], None),                 # copy role read before its copy
            ([('src', 0, 1), ('add', 1, 2)], None),                                # copy from the wrong source
            ([('src', 1, 1), ('add', 1, 0)], None),                                # copy role starts with a V gate
            ([('src', 0, 1)], None),                                               # copy never made
        ]
        for ops, want in cases:
            with self.subTest(ops=ops):
                if want is None:
                    with self.assertRaises(AssertionError):
                        jfdata.kept_copy_ops(ops, copies, s0)
                else:
                    self.assertEqual(jfdata.kept_copy_ops(ops, copies, s0), want)
        # the old pair rule misreads the later gate as a copy; this is the defect the op-index rule fixes
        regated = cases[1][0]
        self.assertEqual(pair_rule(regated, copies, s0), {1, 3})

    def test_two_copies_of_one_source(self):
        ops = [('src', 0, 1), ('add', 1, 0), ('add', 2, 0), ('add', 2, 0)]
        self.assertEqual(jfdata.kept_copy_ops(ops, {0: [1, 2]}, {0: 0}), {1, 2})

    def test_witness2_unchanged(self):
        """on the frozen witness 2 both rules agree, so rounds seven to nine are unaffected"""
        W = jfdata.load()['W']
        new = jfdata.kept_copy_ops(W['ops'], W['copies'], W['s0'])
        self.assertEqual(new, pair_rule(W['ops'], W['copies'], W['s0']))
        self.assertEqual(len(new), sum(len(cs) for cs in W['copies'].values()))


class Layers(unittest.TestCase):
    """the NDS and reuse layers of the round-ten cover profile"""
    def test_nds_rebuild(self):
        # h = 4, m = 12: a gauge of nullity u = 2, radical 1, t = 4 (4 * 3 <= 12) replaces child 3(h - u) = 6 by one
        # child of rank 12 - 8 = 4 at weight 3/4
        H, W = c10.nds_layer({6: Q(2), 1: Q(5)}, Q(10), 4, 12, {'2,1,4': 2})
        self.assertEqual(H, {1: Q(5), 4: Q(3, 2)})
        self.assertEqual(W, Q(10) + 2 * (Q(3, 4) - 1))

    def test_nds_full_share_drops_the_child(self):
        H, W = c10.nds_layer({6: Q(1), 1: Q(1)}, Q(3), 4, 12, {'2,0,6': 1})
        self.assertEqual(H, {1: Q(1)})
        self.assertEqual(W, Q(3) + Q(3, 6) - 1)

    def test_nds_rejects(self):
        for degrees in ({'2,1,5': 1},      # 5 (2 + 1) > 12: the degree is too large
                        {'2,0,2': 1},      # below 3
                        {'3,0,4': 1},      # no exterior child 3(h - 3) = 3 to replace
                        {'2,0,4': 3}):     # more gauges than exterior children
            with self.subTest(degrees=degrees), self.assertRaises(AssertionError):
                c10.nds_layer({6: Q(2), 1: Q(1)}, Q(4), 4, 12, degrees)

    def _base(self):
        bit, cx = c10.c9.pr144_inputs()
        return c10.plain(cx)

    def test_scale(self):
        inv = self._base()
        k = 7
        scaled = dict(inv, scale=k, W_per_vertex=str(Q(inv['W_per_vertex']) * k),
                      child_histogram={r: str(Q(n) * k) for r, n in inv['child_histogram'].items()})
        self.assertEqual(c10.profile10(scaled, 'complex')['hist'], c10.profile10(inv, 'complex')['hist'])

    def test_reuse_rejects(self):
        inv = self._base()
        W = Q(inv['W_per_vertex'])
        # W may not grow, and may not drop by more than one unit per pair
        for W2, pairs in ((W + 1, 5), (W - 6, 5)):
            bad = dict(inv, reuse=dict(pairs=pairs, W_per_vertex=str(W2), child_histogram=inv['child_histogram']))
            with self.subTest(W2=W2), self.assertRaises(AssertionError):
                c10.profile10(bad, 'complex')


class Regressions(unittest.TestCase):
    def test_pr144_and_round9(self):
        self.assertTrue(all(c10.regressions(verbose=False).values()))

    def test_ladder_switch(self):
        # OFF by default (round nine's kappa is reproduced above); ON, rung k satisfies a* - a_k = theta^k (a* - a_old)
        bit, cx = c10.c9.pr144_inputs()
        off = c10.certify(bit, cx)
        on = c10.certify(dict(bit, cert=dict(bit['cert'], ladder=3)), cx)
        b = on['bit']
        self.assertEqual(off['bit']['ladder'], 1)
        self.assertEqual(b['astar'] - b['a_bit'], b['theta'] ** 3 * (b['astar'] - b['a_old']))
        self.assertGreater(on['k']['kappa'], off['k']['kappa'])


@unittest.skipUnless((FROZEN / 'bit_inventory.json.gz').exists(), 'round ten not frozen yet')
class Round10(unittest.TestCase):
    """the round-ten assembly from the frozen inventories (certificates/round10/); round eleven reuses these tests
    with its own folder and Lean input"""
    frozen, lean_name = FROZEN, 'round10-histograms.json'

    @classmethod
    def setUpClass(cls):
        cls.bit = c10.load(str(cls.frozen / 'bit_inventory.json.gz'))
        cls.cx = c10.load(str(cls.frozen / 'cx_inventory.json.gz'))
        cls.r = c10.certify(cls.bit, cls.cx)

    def test_claims_pinned(self):
        self.assertEqual(self.r['bit']['astar'], Q(self.bit['claim']['astar']))
        self.assertEqual(self.r['bit']['a_bit'], Q(self.bit['claim']['a_bit']))
        self.assertEqual(self.r['cx']['a_c'], Q(self.cx['claim']['a_c']))

    def test_grid_points_are_largest(self):
        for d in (self.r['bit'], self.r['cx']):
            self.assertGreater(d['gap'], 0)
            if d['gap_next'] is not None:
                self.assertLessEqual(d['gap_next'], 0)

    def test_kappa(self):
        k = self.r['k']
        self.assertTrue(k['ok'], k['bad'])
        self.assertLess(k['kappa'], k['minimum'])
        self.assertLess(k['minimum'] - k['kappa'], Q(1, c10.KGRID))

    def test_lean_input_matches(self):
        with open(ROOT / 'lean' / self.lean_name) as f:
            d = json.load(f)
        self.assertEqual(Q(*d['bit']['a']), self.r['bit']['astar'])
        self.assertEqual(Q(*d['bit']['stop']['ab']), self.r['bit']['a_bit'])
        self.assertEqual(Q(*d['cx']['a']), self.r['cx']['a_c'])
        self.assertEqual(Q(*d['kappa']), self.r['k']['kappa'])
        for side, key in (('bit', 'bit'), ('cx', 'cx')):
            p = self.r[key]['p']
            self.assertEqual([tuple(x) for x in d[side]['hist']],
                             [(r, int(n * p['lcm'])) for r, n in sorted(p['hist'].items())])

    def test_controls(self):
        # one grid point above a*, a kappa at the minimum margin, and a broken child histogram must all fail
        cert = {k: v for k, v in self.bit['cert'].items() if k != 'astar'}
        bit = dict(self.bit, cert=dict(cert, saving=str(self.r['bit']['astar'] + Q(1, 10**12))))
        bit.pop('claim', None)
        with self.assertRaises(AssertionError):
            c10.certify(bit, self.cx)
        k = c10.c9.assembly(self.r['k']['a'], self.r['cx']['a_c'], self.r['bridge'], kappa=self.r['k']['minimum'])
        self.assertFalse(k['ok'])
        h = dict(self.bit['child_histogram'])
        r0 = min(h, key=int)
        h[r0] = str(Q(h[r0]) + 1)
        with self.assertRaises(AssertionError):
            c10.certify(dict(self.bit, child_histogram=h), self.cx)


if __name__ == '__main__':
    unittest.main()
