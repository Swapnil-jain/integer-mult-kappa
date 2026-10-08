"""Round-seven witness 2, part C (stdlib): B-defer + V frames exactly over Q, late copies, side lemma, certificate.

 Q1. V-leaf starts recomputed EXACTLY from their definition s_i = N_i cap s_{i+1} along the X_S time order (kept
     copies' first frames intersected into the source role's start; fallback <t_S> if G-degenerate), equal to
     the frozen V starts; t_S <= s_i <= N_i; X_S chain <t_S>, early V, deferred V (by dim), F nested exactly.
 Q2. the frozen sigma_u (canonical integer bases) have dimension f_u; exactly: sigma_u <= F0(u), sigma_u <= t_T^perp for every target, G-nondegenerate; deferred readout levels on
     every y_T nested in TIME order (sorted by dim sigma), strict as a set of levels.
 Q3. late copies: every gate at a late-copy frame C (an exact intersection of use frames, part A): C strictly inside
     the copy role's next frame (its use) and the pivot's next frame (it keeps climbing), exactly.
 Q4. side lemma (generic north-east pivots of U^T (P_B - P_A) U^-T and V^T (..) V^-T, mod 2^31-1 at
     check_lifted's integer point) on every distinct step of rank r with 2r > h that the certificate charges.
 Q5. certificate rebuilt from the per-role chains in time (the walk of part B): rank sum s, moment.certify, kappa with the
     round-6 complex side; staircase variant.
 N.  negative controls.
Usage: python3 check_design.py"""
import time, sys, os, random
from fractions import Fraction as Fr
from math import isqrt, comb
from collections import Counter, defaultdict
from xq import rref, ann, inside_ann, red, rref_mod, nondeg, prim
T0 = time.time()
def log(*a): print('[%5.0fs]' % (time.time() - T0), *a, flush=True)
OK = {}
def report(name, ok):
    OK[name] = bool(ok); log('%-104s %s' % (name, 'PASS' if ok else 'FAIL'))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, '..', 'deferred-readout'), os.path.join(HERE, '..', 'two-stage-bit'),
                os.path.join(HERE, '..', '..', 'scripts')]
import jfdata
K = jfdata.load(); EXa, eqida = jfdata.load_frames()
X = K['X']; W = K['W']; E = dict(EX=EXa, eqid=eqida)
ro = jfdata.release_orders(K); B_ = dict(sel=sorted(ro['sel']), late_V=ro['late_V'], defer_order=ro['defer_order'])
B_['walk'] = jfdata.walks(K, eqida, ro)
h = 23; P = 67108859; G = X['G']; trip = [tuple(t) for t in G['trip']]; v = len(trip); R = X['R']; args = G['args']
EX = E['EX']; eqid = E['eqid']; chains = X['chains']; fnode = X['fnode']; f = W['f']; sel = set(B_['sel'])
FULL = tuple(tuple(int(i == j) for j in range(h)) for i in range(h))
tv = lambda T: tuple(int(q in T) for q in range(h))
_Z = {}
def Zof(B):
    if B not in _Z: _Z[B] = ann(B, h)
    return _Z[B]
inside = lambda A, B: len(B) == h or inside_ann(A, Zof(B))
def meet(A, B): return ann(rref(list(Zof(A)) + list(Zof(B)), h), h) if len(B) < h else A
tperp = [ann(rref([[9 * x - 3 for x in tv(T)]], h), h) for T in trip]

# ------------------------------------------------------------------ Q1. V-leaf starts
leafroles = {s for s in range(R) if fnode[s] is not None and args[fnode[s]] is None}
copies = W['copies']; s0 = W['s0']; xorder = W['xorder']; srcop = W['srcop']
def N(s): return EX[chains[s][1]] if len(chains[s]) > 1 else FULL
VS = {}; nfall = 0
for x, order in xorder.items():
    cur = None
    for s in reversed(order):
        B = N(s)
        if s == s0[x]:
            for c in copies.get(x, ()): B = meet(B, N(c))
        cur = B if cur is None else meet(B, cur)
        if not nondeg(cur): cur = rref([tv(trip[x - 1])], h); nfall += 1
        VS[s] = cur
    for c in copies.get(x, ()): VS[c] = VS[s0[x]]
report('Q1. exact V starts for all %d leaf roles (fallbacks %d) == frozen V starts' % (len(leafroles), nfall),
       set(VS) == leafroles == set(K['VS']) and all(VS[s] == K['VS'][s] for s in leafroles))
late = set(B_['late_V'])
badv = sum(1 for s in leafroles if not inside((tv(trip[fnode[s] - 1]),), VS[s]) or not inside(VS[s], N(s)) or not nondeg(VS[s]))
report('Q1. exactly: t_S <= s_i <= N_i (first frame of the use) and s_i G-nondegenerate (all leaf roles)', badv == 0)
badx = 0; xsteps = []; xdata = []
for x, order in xorder.items():
    early = [s for s in order if s not in late]; lt = sorted((s for s in order if s in late), key=lambda s: (len(VS[s]), s))
    seq = [rref([tv(trip[x - 1])], h)] + [VS[s] for s in early + lt] + [FULL]
    for A, B in zip(seq, seq[1:]):
        if len(A) > len(B) or not inside(A, B) or (len(A) == len(B) and A != B): badx += 1
        if len(B) > len(A): xsteps.append((A, B))
    ds = [len(B) for B in seq]; ds = [d for i, d in enumerate(ds) if i == 0 or d != ds[i - 1]]
    xdata.append([b - a for a, b in zip(ds, ds[1:])])
report('Q1. X_S chains nested exactly in TIME order (<t_S>, early V, deferred V by dim, F) on %d sources; xdata == frozen' % len(xorder),
       badx == 0 and sorted(map(tuple, xdata)) == sorted(map(tuple, W['xdata'])))

# ------------------------------------------------------------------ Q2. sigma
def F0ex(s): return VS[s] if s in leafroles else EX[chains[s][0]]
SIG = dict(K['SIG'])
report('Q2. frozen sigma frames: one per deferred role, canonical bases, dims == f (%d roles)' % len(sel),
       set(SIG) == sel and all(rref(SIG[s], h) == SIG[s] and len(SIG[s]) == f[s] for s in sel))
log('   max |entry| in sigma bases %d' % max((abs(x) for B in SIG.values() for r in B for x in r), default=0))
F2 = W['F2']; tg = {s: [t for t in range(v) if F2[s] >> t & 1] for s in sel}
report('Q2. sigma_u <= F0(u) exactly', all(inside(SIG[s], F0ex(s)) for s in SIG))
report('Q2. sigma_u <= t_T^perp exactly for every target T of u', all(inside(SIG[s], tperp[t]) for s in SIG for t in tg[s]))
report('Q2. G nondegenerate on every sigma frame (exact)', all(nondeg(SIG[s]) for s in SIG))
def ychains(order):
    byT = defaultdict(list)
    for s in order:
        for t in tg[s]: byT[t].append(s)
    bad = 0; lev = {}
    for t, ss in byT.items():
        seq = [SIG[s] for s in ss]
        for A, B in zip(seq, seq[1:]):
            if len(A) > len(B) or not inside(A, B) or (len(A) == len(B) and A != B): bad += 1
        lv = []
        for B in seq:
            if not lv or lv[-1] != B: lv.append(B)
        lev[t] = lv
    return bad, lev
bY, lev = ychains(B_['defer_order'])
report('Q2. deferred readout levels on every y_T nested in TIME order (sorted by dim sigma), exactly (%d targets)' % len(lev), bY == 0)
yd = []
for t in range(v):
    ds = sorted(set([0] + [len(B) for B in lev.get(t, [])] + [h - 1])); yd.append([b - a for a, b in zip(ds, ds[1:])])
report('Q2. y_T chains recorded in the frozen schedule == these exact chains', yd == W['yd'])
bYr, _ = ychains(B_['defer_order'][::-1])
report('N. NEG reversed readout order breaks the y_T nesting (%d bad)' % bYr, bYr > 0)
rng = random.Random(5); pert = tot = 0
for s in [s for s in sorted(SIG) if len(F0ex(s)) < h][:300]:
    Bq = [list(r) for r in SIG[s]]; Bq[0] = [x + rng.randint(1, 3) for x in Bq[0]]; tot += 1; pert += not inside(rref(Bq, h), F0ex(s))
report('N. NEG perturbed sigma frames fail sigma <= F0 (%d / %d)' % (pert, tot), tot and pert >= 0.9 * tot)

# ------------------------------------------------------------------ Q3. late copies
ops = W['ops']; opfr = X['opfr']; walk = B_['walk']
# frame ids defined as late-copy intersections (not also a region, root or use frame id; same precedence as part A)
_other = set(X['Mid'].values()) | set(X['rootid']['out'].values()) | set(X['rootid']['ret'].values()) | set(X['dfr'].values())
capid = {i for i in X['capdef'] if i not in _other}
def key_frame(k):
    if isinstance(k, tuple) and k and k[0] == 'sig': return SIG[k[1]]
    if isinstance(k, tuple) and k and k[0] == 'v': return VS[k[1]]
    return EX[k]
def after(s, k):
    w = walk[s]; i = w.index(k); return w[i + 1] if i + 1 < len(w) else None
nlate = 0; badl = Counter(); eqc = 0
for i, op in enumerate(ops):
    if op[0] != 'add': continue
    if opfr[i] not in capid: continue
    # a late copy: gate (copy role r <- pivot) at C = an intersection of use frames (exact, part A). Needed: the copy
    # leaves C upward into its use (C <= next frame of r), the pivot keeps climbing (C <= its next frame); recorded
    # too: whether C is exactly (copy's next) cap (pivot's next), i.e. C_i = D(l_i) cap C_{i+1}.
    nlate += 1; r, piv = op[1], op[2]; Ck = eqid[opfr[i]]
    nr = after(r, Ck); npv = after(piv, Ck)
    if nr is None or npv is None: badl['missing next frame'] += 1; continue
    Dr, Dp = key_frame(nr), key_frame(npv)
    if not (len(EX[Ck]) < len(Dr) and inside(EX[Ck], Dr)): badl['C not strictly inside the copy use frame'] += 1
    if not (len(EX[Ck]) < len(Dp) and inside(EX[Ck], Dp)): badl['pivot does not climb from C'] += 1
    eqc += meet(Dr, Dp) == EX[Ck]
report('Q3. %d late copies at intersection frames: C < copy use frame, C < pivot next frame (exact) %s; C == use cap pivot-next on %d' % (nlate, dict(badl) or '', eqc), nlate > 0 and not badl)
negc = 0
for i, op in enumerate(ops):
    if op[0] == 'add' and opfr[i] in capid:
        comps = [EX[c] for c in X['capdef'][opfr[i]]]
        r = op[1]; Dr = key_frame(after(r, eqid[opfr[i]]))
        rest = [B for B in comps if B != Dr]
        if rest:
            Cn = rest[0]
            for B in rest[1:]: Cn = meet(Cn, B)
            negc += not inside(Cn, Dr)
report('N. NEG copy taken at C_{i+1} instead of C_i leaves its use frame (%d copies)' % negc, negc > 0)

# ------------------------------------------------------------------ chains with dims, certificate data
def kdim(k): return len(key_frame(k))
rk = Counter(); corners = Counter(); steps = {}
def add_step(A, B, kind):
    r = len(B) - len(A)
    if 2 * r > h: steps.setdefault((A, B), (A, B, r, kind))
ZERO = ()
for s in range(R):
    fr = [key_frame(k) for k in walk[s]] + [FULL]
    if s not in sel: fr = [ZERO] + fr
    ds = [len(B) for B in fr]
    assert all(b >= a for a, b in zip(ds, ds[1:]))
    for A, B in zip(fr, fr[1:]):
        if len(B) > len(A): rk[len(B) - len(A)] += 1; add_step(A, B, 'chain')
    corners[h - f[s]] += 1
    add_step(SIG[s] if s in sel else ZERO, FULL, 'corner')
for t, lv in lev.items():
    q = [ZERO] + lv + [tperp[t]]
    for A, B in zip(q, q[1:]): add_step(A, B, 'ydata')
for A, B in xsteps: add_step(A, B, 'xdata')
for c in range(h): add_step(ZERO, ann(rref([[1 - 3 * (q == c) for q in range(h)]], h), h), 'centre')
report('Q5. per-role chain rank == corner rank (sum %d = h R)' % sum(r * n for r, n in rk.items()),
       sum(r * n for r, n in rk.items()) == sum(r * n for r, n in corners.items()) == sum(h - f[s] for s in range(R)))

# ------------------------------------------------------------------ Q4. side lemma
from linalg import inv_mod_matrix, rightmost_pivots, predicted, Q31
from check_lifted import point_UV
q = Q31; U, V = point_UV(h, 1)
conj = []
for A_ in (U, V):
    Ai = inv_mod_matrix(A_, q); assert Ai is not None; conj.append(([list(c) for c in zip(*A_)], Ai))
ID = [[int(i == j) for j in range(h)] for i in range(h)]
inv9 = pow(9, q - 2, q)
def proj_conj(B, cj):
    if not B: return [[[0] * h for _ in range(h)]] * 2
    if len(B) == h: return [ID, ID]
    B = [[x % q for x in r] for r in B]
    BG = [[(x - sum(b) * inv9) % q for x in b] for b in B]
    Gr = [[sum(x * y for x, y in zip(bg, c)) % q for c in B] for bg in BG]
    Gi = inv_mod_matrix(Gr, q); assert Gi is not None, 'Gram singular mod q'
    out = []
    for Acols, Ai in cj:
        LU = [[sum(x * y for x, y in zip(ac, b)) % q for b in B] for ac in Acols]
        RU = [[sum(x * y for x, y in zip(bg, ai)) % q for ai in Ai] for bg in BG]
        Mi = [[sum(Gi[i][k] * RU[k][j] for k in range(len(B))) % q for j in range(h)] for i in range(len(B))]
        out.append([[sum(LU[i][k] * Mi[k][j] for k in range(len(B))) % q for j in range(h)] for i in range(h)])
    return out
_pc = {}
def PC(B, control=False):
    k = (B, control)
    if k not in _pc:
        if len(_pc) > 4000: _pc.clear()
        _pc[k] = proj_conj(B, [([list(c) for c in zip(*ID)], ID)] * 2 if control else conj)
    return _pc[k]
pred = predicted(h)
def lemma(items, control=False):
    bad = 0
    for A, B, r, _ in items:
        QA, QB = PC(A, control), PC(B, control)
        for t in range(2):
            Mx = [[(QB[t][i][j] - QA[t][i][j]) % q for j in range(h)] for i in range(h)]
            if rightmost_pivots(Mx, q) != pred[r]: bad += 1
    return bad
items = sorted(steps.values(), key=lambda x: (x[0], x[1]))
dimq = sum(1 for A, B, _, _ in items for M in (A, B) if M and rank_ok(M) is False) if False else 0
log('   distinct high-rank steps %d %s' % (len(items), dict(Counter(k for *_, k in items))))
badL = 0
for k0 in range(0, len(items), 1000):
    badL += lemma(items[k0:k0 + 1000]); log('   side lemma %d/%d, failures %d' % (min(k0 + 1000, len(items)), len(items), badL))
report('Q4. side lemma on all %d distinct high-rank steps (both conjugations, mod 2^31-1, exact frames)' % len(items), badL == 0)
badc = lemma(items[:300], control=True)
report('N. NEG side lemma with U = V = I fails (%d / 600)' % badc, badc > 0)

# ------------------------------------------------------------------ Q5. certificate
import moment
from certificate_round3 import evaluate
m = h * h
def inner(r): return [1] * (h - r) + [2 * r - h] if 2 * r > h else [1] * r
def histogram(stair=False):
    Nn = v * v; A = v * R; H = Counter()
    for _ in range(2):
        for r, c in corners.items():
            H[m - 2 * r] += v * c
            for w in inner(r): H[w] += v * c
        for r, c in rk.items():
            for w in inner(r): H[w] += v * c
        for _c in range(h):
            for w in inner(h - 1): H[w] += v
        for st in yd:
            for x_ in st:
                for w in inner(x_): H[w] += v
        for xs in xdata:
            for x_ in xs:
                for w in inner(x_): H[w] += v
    H[m - 4 * h + 2] += 2 * Nn
    for w in ([h - 2, h - 5] + [1] * 6 if stair else [h - 2] + [1] * (h + 1)): H[w] += 2 * Nn
    H[1] += Nn
    H = {w: c for w, c in H.items() if c}
    Wd = 2 * Nn + 2 * A; L = 2 * v * h * (h - 1); s_ = Wd * m - Nn + L
    return dict(h=h, m=m, v=v, N=Nn, W=Wd, L=L, s=s_, R=R, hist=H)
c = histogram()
rs = sum(w * n for w, n in c['hist'].items())
report('Q5. rank sum == s = W m - N + L = %d, max child < m' % c['s'], rs == c['s'] and max(c['hist']) < m and c['W'] * m > c['s'])
import deferred as dr
S7 = jfdata.Shape(K); rkC, ylevC, xdC, _ = jfdata.accounting(K, eqida)
cC = dr.histogram(S7, rkC, ylevC, xdC)
report('Q5. rebuilt histogram == the certificate histogram (scripts/certificate_round7.py) of this schedule', c['hist'] == cC['hist'] and c['s'] == cC['s'])
a, _ = moment.certify(c)
from fractions import Fraction as Q_
report('Q5. moment.certify: a_b = %s (%.6e), claim %s' % (a, float(a), K['claimed']['plain']['a_b']), a == Q_(K['claimed']['plain']['a_b']))
k = evaluate(a, Q_(36926111, 5 * 10**11), Q_(1, 1000), 'crude', m_c=576, s_c=119453132304)
report('Q5. kappa = %s (%.7e) ok=%s > 2^-14 = %.7e' % (k['kappa'], float(k['kappa']), k['ok'], 2**-14),
       k['ok'] and k['kappa'] > Q_(1, 2**14) and k['kappa'] == Q_(K['claimed']['plain']['kappa']))
cs = histogram(True); rs2 = sum(w * n for w, n in cs['hist'].items())
a2, _ = moment.certify(cs)
k2 = evaluate(a2, Q_(36926111, 5 * 10**11), Q_(1, 1000), 'crude', m_c=576, s_c=119453132304)
report('Q5. staircase entrance [h-2, h-5, 1^6]: rank sum == s; a_b = %s; kappa = %s (%.7e) ok=%s' % (a2, k2['kappa'], float(k2['kappa']), k2['ok']),
       rs2 == cs['s'] and k2['ok'] and (a2, k2['kappa']) == (Q_(K['claimed']['staircase']['a_b']), Q_(K['claimed']['staircase']['kappa'])))
H2 = dict(c['hist']); H2[1] += 1
report('N. NEG perturbed histogram breaks the rank sum', sum(w * n for w, n in H2.items()) != c['s'])
log('ALL %s' % ('PASS' if all(OK.values()) else 'FAIL: ' + str([k for k, x in OK.items() if not x])))
sys.exit(0 if all(OK.values()) else 1)
