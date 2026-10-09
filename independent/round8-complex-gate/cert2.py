"""Certificate rebuilt from a sreplay histogram pickle {H, W, R, groups, m}: W re-derived as 2N + 2 * groups * R,
rank sum s, deficit, exact moment a_c (certify.py's independent ln/exp bounds), then
certificate_round3.evaluate(6338633/(5*10**10), a_c, 1/1000, 'crude', m_c=576, s_c=s).
Usage: python3 cert2.py hist.pkl"""
import sys, pickle, math, os
from fractions import Fraction as Q
from certify import certify
d = pickle.load(open(sys.argv[1], 'rb')); H = d['H']; m = d['m']; h = int(math.isqrt(m)); v = math.comb(h, 3); N = v * v
W = 2 * N + 2 * d['groups'] * d['R']; assert W == d['W']
s = sum(r * n for r, n in H.items())
print(dict(W=W, s=s, D=W * m - s, maxrank=max(H), guard_ok=max(H) <= m - 1, widths=len(H)))
a, Mup, Mlo = certify(H, W, m)
print(dict(a_c=str(a), a_c_float=float(a), upper_M_at_a=float(Mup), lower_M_next=float(Mlo), certified=Mup < 1, largest=Mlo > 1))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
import certificate_round3 as C
r = C.evaluate(Q(6338633, 5 * 10 ** 10), a, Q(1, 1000), 'crude', m_c=576, s_c=s)
print(dict(ok=r['ok'], kappa=str(r['kappa']), kappa_float=float(r['kappa']), binding=r['binding'], bad=r['bad'], eps=float(r['eps'])))
