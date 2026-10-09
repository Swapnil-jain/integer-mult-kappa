"""Independent endpoint check (Z_w / i^9 / P^-1 child), written from the two-stage transfer statement
A = -F y, B = F P^-2 x' + E y (E = F P^-1), x' = Z_w x, B' = B + P^-1 A, Y_out = i^9 Z_w B', X_out = -A.
Dual picture: F, P = C_<w>, E diagonal with exponents q_F(e) = popc(e), q_P(e) = popc((w.e) w) = 9 (w.e), mod 4;
Z_w acts as e -> e + w.  (1) all 2^16 addresses at h = 4 for every weight-9 w = t1 (x) t2 (all 16 triple pairs),
complex vectors of random Gaussian integers; (2) h = 24, random eta and random pairs with literal 576-bit w.
Controls: drop pre Z_w, drop post Z_w, drop i^9, P instead of P^-1, P^-2 -> I.
Usage: python3 endpoint.py"""
import random, itertools
import numpy as np
from gf2 import kron, pc, dot
I4 = np.array([1, 1j, -1, -1j])
def run(h, ws, addrs, ctl=None, rng=np.random.default_rng(3)):
    bad = 0
    for w in ws:
        n = len(addrs); x = rng.integers(-9, 10, n) + 1j * rng.integers(-9, 10, n); y = rng.integers(-9, 10, n) + 1j * rng.integers(-9, 10, n)
        idx = {e: k for k, e in enumerate(addrs)}
        qF = lambda e: pc(e) % 4; qP = lambda e: (9 * dot(w, e)) % 4
        # x' = Z_w x : x'(e) = x(e + w)
        xp = np.array([x[idx[e ^ w]] for e in addrs]) if ctl != 'nopre' else x
        A = np.array([-I4[qF(e)] for e in addrs]) * y
        B = np.array([I4[(qF(e) - (0 if ctl == 'char2' else 2 * qP(e))) % 4] for e in addrs]) * xp + np.array([I4[(qF(e) - qP(e)) % 4] for e in addrs]) * y
        Bp = B + np.array([I4[(qP(e) if ctl == 'P' else -qP(e)) % 4] for e in addrs]) * A
        Yo = I4[0 if ctl == 'noi9' else 1] * np.array([Bp[idx[e ^ w]] if ctl != 'nopost' else Bp[idx[e]] for e in addrs])   # i^9 = i
        Xo = -A
        Fx = np.array([I4[qF(e)] for e in addrs]) * x; Fy = np.array([I4[qF(e)] for e in addrs]) * y
        bad += int((np.abs(Yo - Fx) > 0).sum() + (np.abs(Xo - Fy) > 0).sum())
    return bad
h = 4; T = [sum(1 << p for p in t) for t in itertools.combinations(range(h), 3)]
ws = [kron(a, b, h) for a in T for b in T]; assert all(pc(w) == 9 for w in ws)
full = list(range(1 << 16))
print('h=4 all 2^16 addresses, all 16 w: wrong =', run(4, ws, full))
for c in ('nopre', 'nopost', 'noi9', 'P', 'char2'): print('  control', c, 'wrong =', run(4, ws[:2], full, c))
# h = 24: the identity is pointwise in (e, e+w); check it on random addresses with literal 576-bit w
rnd = random.Random(5); h = 24; T = [sum(1 << p for p in t) for t in itertools.combinations(range(h), 3)]
bad = 0
for _ in range(20000):
    w = kron(rnd.choice(T), rnd.choice(T), h); e = rnd.getrandbits(h * h)
    # Y_out(e) = i^9 B'(e+w) with B'(f) = F(f) P(f)^-2 x'(f) + [E(f) - P(f)^-1 F(f)] y(f), x'(f) = x(f+w)
    f = e ^ w
    cx = (9 + pc(f) - 2 * 9 * dot(w, f)) % 4          # coefficient of x(e) in Y_out(e) must be F(e)
    cy_terms = ((pc(f) - 9 * dot(w, f)) % 4, (2 + pc(f) - 9 * dot(w, f)) % 4)   # E y and P^-1 A = -P^-1 F y cancel
    bad += (cx != pc(e) % 4) + ((cy_terms[0] - cy_terms[1]) % 4 != 2)
print('h=24 random (e, w) pointwise identity: wrong =', bad)
