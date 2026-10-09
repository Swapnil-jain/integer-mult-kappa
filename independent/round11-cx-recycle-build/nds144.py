"""Nullity-dependent sharing (nullshare/) on #144's cube word: per gauged role (u = h - d, radical r of sigma),
degree t = floor(m/(u+r)) (degenerate extension when r > 0), weight 3/t per vertex, tail m - t u; ungauged roles t=3."""
import pickle, sys, json
from fractions import Fraction as Q
from collections import Counter
sys.path.insert(0, '.')
from cube_word import build, compile_frames, select_gauges, children, published, basis
from gtypes import classify
from evalc import acert, D_of
h = 24; m = 72
def nds_children(pr, sel, degen=True, cap=None):
    h = pr['h']; m = 3 * h
    C = Counter({r: Q(c) for r, c in children(dict(pr, gauges={})).items()})
    W = Q(2 * pr['v'] + pr['R']); deg = Counter()
    for z in sel:
        u, r, alt = classify(z['A'], h)
        if alt or (r and not degen): t = 3
        else: t = max(3, m // (u + r))
        if cap: t = min(t, cap)
        deg[(u, r, t)] += 1
        W += Q(3, t) - 1
        if m - t * u: C[m - t * u] += Q(3, t)
    return C, W, deg
def tail_nds(d, A, degen=True):
    u, r, alt = classify(A, h)
    t = 3 if (alt or (r and not degen)) else max(3, m // (u + r))
    return [(m - t * u, Q(3, t))] if m - t * u else []
if __name__ == '__main__':
    W0 = pickle.load(open('word144.pkl', 'rb'))
    pr, sel = W0['pr'], W0['sel']
    C = children(pr); W = 2 * pr['v'] + pr['R']
    a, rf = acert(C, W); print('#144 reproduced: D', D_of(C, W), 'a_c', a, float(a), rf, flush=True)
    for degen in (False, True):
        C2, W2, deg = nds_children(pr, sel, degen)
        a, rf = acert(C2, W2); print('NDS on #144 gauges degen=%s: D %s W %s a_c %s = %.10e  degrees %s' % (degen, D_of(C2, W2), W2, a, float(a), dict(deg)), flush=True)
