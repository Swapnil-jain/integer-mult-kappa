"""R5 census on the released word (pf_G37_rx): every ungauged role that could still take a gauge (untouched by the
centre cone, has a first op), its admissible gauge dim d (sigma = (ann[first op] + its targets' current limits)^perp,
the same rule as reselect, evaluated against the FINAL selection), and an optimistic bound on the gain if every role
with d >= DMIN were gauged AND paired for free: no target-split cost, a donor with the best last-frame dim e <= d always
available. Each such role then changes 3 x child r -> r - d (its first step), 3 x child (h - e) -> (d - e) (the
donor's exit) and W -> W - 1, which keeps the deficit. The bound is the float root of the edited histogram."""
import sys, os, pickle, math
from collections import Counter
from resel_br import *
import make_word
from score_cp import froot
DMIN = int(os.environ.get('DMIN', 6))
W = make_word.make(11, 'pf_G37_rx'); PL = pickle.load(open('gated_plan2_pf_G37_rx_p11.pkl', 'rb'))
sel = PL['sel']; pairs = PL['pairs']; tau = PL['tau']
W2 = with_sel(W, sel); F = dict(gfirst_frames(W2), frames=PL['frames'])
C, Wv, inf = recount(W2, F, pairs, tau, True); h, m = W['h'], W['m']; a_now = froot(C, Wv, m)
G = W['wit']['g']; s = W['s']; ann = W['wit']['ann']; first = s['first']; touched = s['touched']; R = W['pr']['R']
roots = G['roots']; inputs = G['inputs']; v = W['pr']['v']; ops = s['ops']
co = [0] * R
for r, sr in zip(roots, s['rootroles']): co[sr] |= sum(1 << t for t in r['targets'])
for a, b, x in reversed(ops): co[b] |= co[a]
limit = [None] * v
for r in roots:
    if r['kind'] == 'center': continue
    A = basis(inputs[t] for t in r['targets'])
    for t in r['targets']:
        if limit[t] is None: limit[t] = A
for t in range(v):
    if limit[t] is None: limit[t] = (inputs[t],)
for z in sel:
    for t in z['targets']: limit[t] = tuple(z['A'])
gauged = {z['role'] for z in sel}; rootroles = set(s['rootroles'])
cands = [x for x in range(R) if x not in touched and first[x] is not None and x not in gauged]
D = Counter(); rows = []; NEG = []
for x in cands:
    A = tuple(ann[first[x]]); bits = co[x]
    while bits:
        low = bits & -bits; t = low.bit_length() - 1; bits ^= low; A = basis(A + tuple(limit[t]))
        if len(A) == h: break
    d = h - len(A); r = h - len(ann[first[x]]); D[d] += 1
    tg = [t for t in range(v) if co[x] >> t & 1]
    exa = lambda t: t * math.expm1(a_now * math.log(m / t)) if t > 0 else 0.0
    split = 3 * sum(exa(d) + exa(h - len(limit[t]) - d) - exa(h - len(limit[t])) for t in tg)
    e = max(range(0, d + 1), key=lambda e: exa(h - e) - exa(d - e))
    gain = 3 * (exa(r) - exa(r - d)) + 3 * (exa(h - e) - exa(d - e)) - split
    if d >= DMIN: rows.append((r, d)); NEG.append(gain > 0)
print('ungauged gaugeable roles: %d (touched by centre cone or never used: %d); gauge dim histogram %s' %
      (len(cands), R - len(cands) - len(gauged), dict(sorted(D.items()))))
print('roles with d >= %d: %d; first-frame ranks %s' % (DMIN, len(rows), dict(Counter(r for r, _ in rows))))
C2 = Counter({k: float(n) for k, n in C.items()}); W2v = float(Wv)
ex = lambda t, a=a_now: t * math.expm1(a * math.log(m / t)) if t > 0 else 0.0
for r, d in rows:
    e = max(range(0, d + 1), key=lambda e: ex(h - e) - ex(d - e))
    for k, dlt in ((r, -3), (r - d, 3), (h - e, -3), (d - e, 3)):
        if k > 0: C2[k] += dlt
    W2v -= 1
a_bound = froot(C2, W2v, m)
print('of the %d roles, %d have a positive gain after their target-split cost (reselect rule, free donor)' % (len(rows), sum(NEG)))
print('a now %.6e; optimistic bound with all %d gauged+paired free: %.6e (%+.3f%%)' % (a_now, len(rows), a_bound, 100 * (a_bound / a_now - 1)))
