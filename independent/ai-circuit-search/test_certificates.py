"""Adversarial checks for rational enclosures and accepted certificates."""
import copy
import json
import unittest
from fractions import Fraction as Q
from pathlib import Path
import verify as v
from community_assembly import assembly, InvalidAssembly

DATA=json.loads(Path(__file__).with_name('certificate.json').read_text())


class CertificateTests(unittest.TestCase):
    def test_log_series_against_rational_identities(self):
        for x in (Q(1),Q(2),Q(3,2),Q(529,23),Q(575)):
            lo,hi=v.log_bounds(x)
            doubled=v.log_bounds(x*x)
            self.assertLessEqual(doubled[0],2*hi)
            self.assertGreaterEqual(doubled[1],2*lo)
            self.assertLess(hi-lo,Q(1,10**55))

    def test_exp_zero_and_monotonic_enclosure(self):
        self.assertEqual(v.exp_bounds(0),(Q(1),Q(1)))
        for x in (Q(1,10**6),Q(1,1000),Q(1,2)):
            l,u=v.exp_bounds(x)
            l2,u2=v.exp_bounds(x,degree=18)
            self.assertLess(l,l2)
            self.assertLess(u2,u)

    def test_corrupted_counts_rejected(self):
        n=copy.deepcopy(DATA['networks']['jain_bit']);n['hist']['1']+=1
        with self.assertRaises(ValueError):v.normalize_network(n)

    def test_next_grid_saving_rigorously_fails(self):
        for name,n in DATA['networks'].items():
            a=Q(DATA['root_certificates'][name]['upper'])
            self.assertGreater(v.moment_bounds(n,a)[0],1)

    def test_above_fixed_profile_limit_fails_assembly(self):
        p=copy.deepcopy(DATA['jain_assemblies']['power_split_bit_sharpened'])
        p['kappa']=str(Q(p['a_bit']))
        with self.assertRaises(ValueError):v.check_jain_assembly(p)

    def test_community_negative_controls(self):
        p=DATA['community_assemblies']['sharpened'];b=DATA['community_finite_bridge']
        for flag in ('old_guard','old_exposures','original_prefix'):
            with self.assertRaises(InvalidAssembly):
                assembly(b,Q(p['a_bit']),Q(p['kappa']),h=Q(p['h']),a_complex=Q(p['a_complex']),**{flag:True})

    def test_stale_bridge_rejected(self):
        b=copy.deepcopy(DATA['community_finite_bridge']);b['bit']['maxchild']=1
        with self.assertRaises(InvalidAssembly):assembly(b,Q(1,10**5),Q(1,10**6))


if __name__=='__main__':unittest.main(verbosity=2)
