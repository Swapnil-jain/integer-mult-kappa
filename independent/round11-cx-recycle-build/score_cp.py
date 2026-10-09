"""Fast cube-prefix score of a word config (caps frames, no descent): relaxed arcs, pair-aware gauges, exact birth
matching, chain recount, float root of the moment. Returns 0.0 for an illegal word."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cube_word as cw
from make_word import build_cfg
from relaxed import relaxed_arcs
from resel_br import reselect, with_sel, gfirst_frames, late_tau, candidates, match, recount
def froot(C, W, m):
    items = [(float(n), r / m) for r, n in C.items() if n]; W = float(W); lo, hi = 0.0, 0.05
    for _ in range(60):
        mid = (lo + hi) / 2; s = math.fsum(n * x ** (1 - mid) for n, x in items) / W
        lo, hi = (mid, hi) if s < 1 else (lo, mid)
    return lo
def score(p, cfg, a0=5.6e-4):
    m = 6 * p
    try:
        g = build_cfg(p, cfg); cw.RELAXED = True
        try: prof, wit = cw.compile_frames(g, lambda D: relaxed_arcs(D, verbose=False))
        finally: cw.RELAXED = False
        pr, s = cw.select_gauges(prof, wit)
        W = dict(p=p, h=2 * p, m=m, g=g, prof=prof, wit=wit, pr=pr, s=s)
        F0 = gfirst_frames(W); sel, _, _ = reselect(W, F0, a0, set(), True)
        W2 = with_sel(W, sel); F = gfirst_frames(W2); tau = late_tau(W2, F)
        pairs = match(candidates(W2, F, tau, a0, True, True))
        C, Wv, inf = recount(W2, F, pairs, tau, True)
        return 0.0 if inf['bad'] else froot(C, Wv, m)
    except (AssertionError, TypeError, KeyError, IndexError): return 0.0
