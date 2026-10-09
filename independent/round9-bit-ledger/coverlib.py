"""Bit cover accounting (own code). Histograms per cell; exact moment with the full rare-class fallback surcharge
(b = 1e-16 bad fraction per edge, 32 m^2 children of width 1 each, nothing displaced), rigorous rational certificate on a
10^-12 grid (rational ln upper bound, exp(u) <= 1+u+u^2/(2(1-u/3))), stopped saving a_bit = (1-theta)a* + theta a_old."""
import json, os, math
from fractions import Fraction as Q
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
THETA, A_OLD, BAD = Q(1, 1000), Q(384599, 10**10), Q(1, 10**16)

def ln_upper(x, n=28):
    x = Q(x); assert x > 1; k = 0
    while x > 2: x /= 2; k += 1
    def ta(z):
        return 2 * (sum(z ** (2 * j + 1) / (2 * j + 1) for j in range(n)) + z ** (2 * n + 1) / ((2 * n + 1) * (1 - z * z)))
    return k * ta(Q(1, 3)) + (ta((x - 1) / (x + 1)) if x > 1 else 0)

def with_fallback(c):
    """add the whole fallback: every edge (child) also spawns, with weight BAD, 32 m^2 width-1 children"""
    E = sum(c['hist'].values()); return dict(c, fb=BAD * 32 * c['m'] ** 2 * E)

def F_float(c, a):
    m = c['m']; t = math.fsum(n * (w / m) ** (1 - a) for w, n in c['hist'].items())
    return (t + float(c.get('fb', 0)) * (1 / m) ** (1 - a)) / c['W']

def F_upper(c, a, lu):
    m = c['m']; tot = Q(0)
    items = list(c['hist'].items()) + ([(1, c['fb'])] if c.get('fb') else [])
    for w, n in items:
        u = a * lu[w]; assert 0 <= u < 1
        tot += Q(n) * Q(w, m) * (1 + u + u * u / (2 * (1 - u / 3)))
    return tot / c['W']

def root(c):
    lo, hi = 0.0, 0.05
    for _ in range(200):
        mid = (lo + hi) / 2; lo, hi = (mid, hi) if F_float(c, mid) < 1 else (lo, mid)
    return lo

def certify(c, den=10 ** 12):
    m = c['m']; assert max(c['hist']) < m and min(c['hist']) >= 1
    assert c['W'] * m > sum(w * n for w, n in c['hist'].items())
    lu = {w: ln_upper(Q(m, w)) for w in set(c['hist']) | {1}}
    r = root(c); lo, hi = max(0, int(r * den) - 50), int(r * den) + 2
    assert F_upper(c, Q(lo, den), lu) < 1 and not F_upper(c, Q(hi, den), lu) < 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if F_upper(c, Q(mid, den), lu) < 1: lo = mid
        else: hi = mid
    return Q(lo, den), r

def stopped(a): return (1 - THETA) * a + THETA * A_OLD

def kappa(b_c, a_bit):
    """PR130/131/137 assembly (audit-137/work/kap.py): binding b = min(complex a_c, bit a_bit / (1-beta))."""
    beta, eta = Q(1, 10**6), Q(1, 10**8)
    b = min(Q(b_c), Q(a_bit) / (1 - beta))
    a = (1 - beta) * b - Q(1, 10**10); q0 = a * (1 - 2 * eta); c0 = q0 * (1 + eta)
    return (1 - eta) / (1 + c0 + q0) * q0

def loadlocal(name):
    """'pr97' (PR #130 input, read as data) or 'w1'/'w2' (our frozen witnesses, localhist.json)"""
    if name == 'pr97':
        P = json.load(open(os.path.join(HERE, 'in', 'pr130_bit_input.json')))
        L = dict(h=P['h'], v=P['v'], R=P['R'], l=P['center_loss'], local=Counter(), ex={int(k): x for k, x in P['exit_nullity_histogram'].items()})
        for cls in ('aux', 'data', 'center'):
            for r, x in P['local_rank_histograms'][cls].items(): L['local'][int(r)] += x
        return L
    J = json.load(open(os.path.join(HERE, 'localhist.json')))[name]
    L = dict(h=J['h'], v=J['v'], R=J['R'], l=J['center_loss'], local=Counter(), ex={int(k): x for k, x in J['exit_nullity'].items()})
    for cls in ('aux', 'data', 'center'):
        for r, x in J[cls].items(): L['local'][int(r)] += x
    return L

def cover(L, k=1, m=None, share=None):
    """Three-stage cover, one cell = k vertices sharing aux banks (k=1: PR #130; k=3, m=3h: PR #137).
    share: optional {u: k_u} per-nullity sharing degree (each k_u divides into the cell); default k for all slots.
    Data: 2v roles per vertex, finishing child of width p = m - (3h-2) per data role when p > 0.
    Aux: per stage, slot with nullity u and degree k_u: k/k_u roles per cell, each with exterior m - k_u u."""
    h, v, R = L['h'], L['v'], L['R']; m0 = 3 * h - 2; m = m or max(m0, k * h); p = m - m0; assert p >= 0
    H = Counter(); Wn = Q(2 * k * v)
    for r, x in L['local'].items(): H[r] += 3 * k * x
    if p: H[p] += 2 * k * v
    for u, n in L['ex'].items():
        ku = (share or {}).get(u, k); assert k % ku == 0 and ku * u <= m
        groups = Q(k, ku) * n * 3; Wn += groups
        if m - ku * u: H[m - ku * u] += groups
    assert all(Q(x).denominator == 1 for x in H.values()) and Wn.denominator == 1
    H = {r: int(x) for r, x in H.items() if x}
    s = sum(r * x for r, x in H.items()); W = int(Wn)
    return dict(m=m, W=W, s=s, D=W * m - s, hist=H, k=k)

def lam(c): return sum(n * w * math.log(c['m'] / w) for w, n in c['hist'].items())

def cover_var(L, m, tmax=None):
    """Nullity-dependent orbit sharing (our lever): per stage and vertex, a slot of nullity u gets 1/t_u roles, where
    t_u = min(floor(m/u), tmax) vertices share one role through an order-t_u free right action; exterior m - t_u*u.
    Histograms are per t-lcm cell; returned with counts scaled to integers by the lcm."""
    h, v = L['h'], L['v']; m0 = 3 * h - 2; p = m - m0; assert p >= 0
    t = {u: min(m // u, tmax or 10**9) for u in L['ex']}
    from math import lcm
    K = 1
    for x in t.values(): K = lcm(K, x)
    H = Counter(); W = 2 * K * v
    for r, x in L['local'].items(): H[r] += 3 * K * x
    if p: H[p] += 2 * K * v
    for u, n in L['ex'].items():
        g = K // t[u] * n * 3; W += g
        if m - t[u] * u: H[m - t[u] * u] += g
    H = {r: x for r, x in H.items() if x}; s = sum(r * x for r, x in H.items())
    return dict(m=m, W=W, s=s, D=W * m - s, hist=H, k=K, t=t)
