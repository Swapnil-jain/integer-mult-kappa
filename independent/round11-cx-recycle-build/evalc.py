"""Exact complex saving for a per-vertex child multiset (Fractions allowed) with W: rigorous 1e-12 grid floor
(cover-bit/coverlib.certify: rational ln upper bounds, exp upper bound) plus float root."""
import sys, os
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coverlib import certify, F_upper, ln_upper
def acert(C, W, m=72):
    C = {r: Q(n) for r, n in C.items() if n}
    a, r = certify(dict(hist=C, m=m, W=Q(W)))
    return a, r
def D_of(C, W, m=72): return Q(W) * m - sum(r * Q(n) for r, n in C.items())
