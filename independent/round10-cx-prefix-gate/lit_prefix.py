"""Items 1, 2, 4(complex analogue), 6: own end-to-end replay of #144's paired-cube local word with arbitrary dirty scratch.
Nothing from the PR runs; the circuit is rebuilt here from the note's definitions (our own naive balanced-sum modules
instead of the PR117 restrictions, so the ledger is not theirs; the mechanism is).

Geometry: h = 2p address bits, coordinate 2i+b = pair i selector b; port S = one coordinate from each of 3 distinct pairs,
q_S = chi_S. Signed DAG: per cube E (12), A (6), G (6, signed), F (1); per target cube J the 35 root calls (F / face /
edge aggregates, sums over the right source cubes) and the h coordinate-star centres S_c read through copied centres at
frame 0 with coefficients 1/3 (c in T) and -1/6. K: two 4x4 Hadamard/2 blocks per cube on the original source registers.

Frames (arbitrary-subspace interface of the note): T_U = K_G C_{E_r} K_G^{-1}, K_G = D_M P_G (P_G: |a> -> |G a>, D_M the
quadratic phase i^{q_M(a)}, M = I + (G G^T)^{-1}), C = sqrt(X) = (1+i)/2 [[1,-i],[-i,1]] (#130), G = basis of U completed.
A register 'at frame U' holds T_U applied to its logical value. A transition U -> V is applied LITERALLY as T_V T_U^{-1};
a gate requires both registers at the IDENTICAL frame (canonical subspace) and is a scalar Gaussian-dyadic combination
(coefficients with denominators 2 and 3: all data and scratch are multiples of 3, so 1/3 is exact).
Node frame Phi_x = (annihilator of everything x reaches)^perp; roles climb monotonically (asserted), ranks counted with
d(L_U, L_V) = dim U + dim V - 2 dim(U cap V) and, in lit mode, MEASURED as log2 column support on sample columns.

Schedule (the note's): t=0 old-value reads of gauge-0 roles at frame 0 (exact adjoint coefficients); every root use of a
source injected from the source at its rank-one frame <q_S>; phase 1 = centre closure, centre copies -> 0, scatter reads;
phase 2 by frame dimension (gauged old reads at sigma first, then computes, then root reads); K: sources -> U_par (rank 2),
Hadamard/2 mix at the identical representative, -> q_T^perp (rank h-4), added to their opposite-parity target at that
identical representative, -> F (rank 1), inverse blocks at F; all roles -> F; inverse word at F; sources subtracted at F.
Checks: targets end at T_{q_T^perp}(y + x_T) (i.e. y += (K+H+B)x = x), sources at F x, every aux role at F T_sigma^{-1} w0
(w0 its arbitrary physical start), every gate at an identical frame, every chain nested, histograms.
Modes: lit (physical 2^h registers, exact) | scal (frame-tagged exact Fraction scalars, k parallel trials).
Gauges: 'gauge' arg selects entrance gauges (sigma = first frame cap target limits, PR-style greedy, all accepted);
'omit=K' then moves K of them (every 3rd) back to the frame-0 prelude (gauge omission, complex analogue of item 4).
Controls (must FAIL): noKundo, badgroup, rep, noprelude, lateread, skipclone.
Usage: python3 cube_lit.py p lit|scal SEED [gauge] [omit] [control=...]"""
import sys, itertools, random, time
from fractions import Fraction as Q
from collections import Counter, defaultdict
import numpy as np

p = int(sys.argv[1]); MODE = sys.argv[2]; SEED = int(sys.argv[3])
FLAGS = set(a for a in sys.argv[4:] if '=' not in a); KV = dict(a.split('=') for a in sys.argv[4:] if '=' in a)
CTRL = KV.get('control', '')
rnd = random.Random(SEED); h = 2 * p; N = 1 << h; t0 = time.time()
OUT = {}
def check(name, ok):
    OUT[name] = bool(ok); print(('PASS ' if ok else 'FAIL ') + name, flush=True)

# ---------------- F2 subspaces ----------------
def red(vs):
    b = {}
    for x in vs:
        for p_ in sorted(b, reverse=True):
            if x >> p_ & 1: x ^= b[p_]
        if x:
            p_ = x.bit_length() - 1
            for k in list(b):
                if b[k] >> p_ & 1: b[k] ^= x
            b[p_] = x
    return b
def canon(vs): return tuple(sorted(red(vs).values()))
def dim(U): return len(U)
def dot(a, b): return bin(a & b).count('1') & 1
def inside(A, B):
    bb = red(B)
    for x in A:
        for p_ in sorted(bb, reverse=True):
            if x >> p_ & 1: x ^= bb[p_]
        if x: return False
    return True
def perp(vs):
    piv = red(vs); out = []
    for f in range(h):
        if f in piv: continue
        x = 1 << f
        for c, y in piv.items():
            if y >> f & 1: x |= 1 << c
        out.append(x)
    return canon(out)
def meet(U, V): return canon(perp(list(perp(U)) + list(perp(V))))
def dL(U, V): return len(U) + len(V) - 2 * len(meet(U, V))
Z0 = (); FULL = canon([1 << j for j in range(h)])

# ---------------- ports, cubes ----------------
cubes = list(itertools.combinations(range(p), 3))
bl = list(itertools.product(range(2), repeat=3))
ports = [(I, b) for I in cubes for b in bl]; v = len(ports); idx = {pt: k for k, pt in enumerate(ports)}
coords = [frozenset(2 * i + b for i, b in zip(*pt)) for pt in ports]
q = [sum(1 << c for c in cs) for cs in coords]
par = [sum(b) % 2 for _, b in ports]
def src(I, b): return idx[(I, tuple(b))]

# ---------------- signed DAG ----------------
args = []; sgn = []; sup = []; intern = {}
def node_in(s):
    args.append(None); sgn.append(1); sup.append(frozenset([s])); return len(args) - 1
INP = [node_in(s) for s in range(v)]
def add(a, b, s=1):
    if a is None: assert s == 1; return b
    if b is None: return a
    key = (a, b, s)
    if key not in intern:
        assert not (sup[a] & sup[b]); args.append((a, b)); sgn.append(s); sup.append(sup[a] | sup[b]); intern[key] = len(args) - 1
    return intern[key]
def tsum(xs):
    xs = [x for x in xs if x is not None]
    if not xs: return None
    while len(xs) > 1: xs = [add(xs[i], xs[i + 1]) if i + 1 < len(xs) else xs[i] for i in range(0, len(xs), 2)]
    return xs[0]
Fn, An, Gn = {}, {}, {}
for I in cubes:
    E = {}
    for k, l in itertools.combinations(range(3), 2):
        m_ = 3 - k - l
        for a, b in itertools.product(range(2), repeat=2):
            b0 = [0] * 3; b0[k], b0[l] = a, b; b1 = list(b0); b1[m_] = 1
            E[k, l, a, b] = add(INP[src(I, b0)], INP[src(I, b1)])
        Gn[I, I[k], I[l], 0] = add(E[k, l, 0, 0], E[k, l, 1, 1], -1)
        Gn[I, I[k], I[l], 1] = add(E[k, l, 0, 1], E[k, l, 1, 0], -1)
    for k in range(3):
        l = 1 if k == 0 else 0
        for a in range(2):
            key = (lambda b: (k, l, a, b)) if k < l else (lambda b: (l, k, b, a))
            An[I, I[k], a] = add(E[key(0)], E[key(1)])
    Fn[I] = add(An[I, I[0], 0], An[I, I[0], 1])
roots = []      # dict(node, targets, coef(list), kind, frame)
# cube-prefix flags: 'prefix' builds every edge aggregate from one shared nested-prefix all-but-one module per
# (pair, mode) (restrict.allbut_prefix, #168's module), instead of a fresh balanced tsum per root; 'merge' reads
# face2 + (2b0-1) edge02 as one signed node per (b0, b2), coefficient 1/2, shared by the two ports b1 = 0, 1.
QPRE = {}
def allbut_prefix(n):
    """Nested-prefix all-but-one module (#168's construction, own code): prefix sums p_{k+1} = p_k + x_k, suffix sums
    s_k = x_k + s_{k+1}; y_0 = s_1, y_{n-1} = p_{n-1}, y_{i+1} = p_i + (x_i + s_{i+2}). 3n - 6 additions; each prefix
    p_i (i >= 2) feeds both p_{i+1} and y_{i+1}, with nested supports, which is what lets more carrier arcs in."""
    args = [None] * n; sup = [1 << i for i in range(n)]; by = {x: i for i, x in enumerate(sup)}
    def add(a, b):
        s = sup[a] | sup[b]; assert not sup[a] & sup[b]
        if s not in by: by[s] = len(args); args.append([a, b]); sup.append(s)
        return by[s]
    P = {1: 0}
    for k in range(1, n - 1): P[k + 1] = add(P[k], k)
    S = {n - 1: n - 1}
    for k in range(n - 2, 0, -1): S[k] = add(k, S[k + 1])
    roots = [None] * n; roots[0] = S[1]; roots[n - 1] = P[n - 1]
    for i in range(n - 2):
        inner = add(i, S[i + 2])
        roots[i + 1] = add(P[i], inner) if i >= 1 else inner
    full = (1 << n) - 1
    assert all(sup[r] == full ^ (1 << i) for i, r in enumerate(roots))
    return dict(input_count=n, args=args, roots=roots)
if 'prefix' in FLAGS:
    am = allbut_prefix(p - 2)
    for i, j in itertools.combinations(range(p), 2):
        others = [c for c in range(p) if c not in (i, j)]
        for mode in range(2):
            img = [Gn[tuple(sorted((i, j, c))), i, j, mode] for c in others]
            for a_, b_ in am['args'][len(others):]: img.append(add(img[a_], img[b_]))
            for c, r_ in zip(others, am['roots']): QPRE[i, j, c, mode] = img[r_]
for J in cubes:
    allt = [src(J, b) for b in bl]
    D = tsum(Fn[I] for I in cubes if not set(I) & set(J))
    groups = [('F', allt, D, Q(1, 2))]
    def grp(keybits):
        for fix in itertools.product(range(2), repeat=len(keybits)):
            yield fix, [src(J, b) for b in bl if all(b[kk] == f for kk, f in zip(keybits, fix))]
    def face(k, kb):
        for fix, tg in grp(kb):
            ek = fix[kb.index(k)]
            groups.append(('face%d' % k, tg, tsum(An[I, J[k], 1 - ek] for I in cubes if set(I) & set(J) == {J[k]}), Q(1, 2)))
    def edge(k, l, kb):
        for fix, tg in grp(kb):
            ek, el = fix[kb.index(k)], fix[kb.index(l)]
            third = [c for c in J if c not in (J[k], J[l])][0]
            nd_ = QPRE[J[k], J[l], third, ek ^ el] if QPRE else tsum(Gn[I, J[k], J[l], ek ^ el] for I in cubes if set(I) & set(J) == {J[k], J[l]})
            groups.append(('edge%d%d' % (k, l), tg, nd_, Q(2 * ek - 1, 2)))
    face(0, (0,)); face(1, (0, 1)); edge(0, 1, (0, 1)); face(2, (0, 1, 2)); edge(0, 2, (0, 1, 2)); edge(1, 2, (0, 1, 2))
    if 'merge' in FLAGS:
        f2 = {tuple(tg): nd for nm, tg, nd, co in groups if nm == 'face2'}
        e02 = {tuple(tg): (nd, co) for nm, tg, nd, co in groups if nm == 'edge02'}
        rest_ = [g_ for g_ in groups if g_[0] not in ('face2', 'edge02')]
        merged = []
        for tg in f2:
            nd, co = e02[tg]; sg = 1 if co > 0 else -1
            merged.append(('merge', list(tg), add(f2[tg], nd, sg), Q(1, 2)))
        k12 = [i_ for i_, g_ in enumerate(rest_) if g_[0] == 'edge12'][0]
        groups[:] = rest_[:k12] + merged + rest_[k12:]
    for name, tg, nd, co in groups:
        if nd is None: continue
        roots.append(dict(node=nd, targets=tg, coef=[co] * len(tg), kind='side', name=name, frame=perp([q[t] for t in tg])))
centres = []
for c in range(h):
    i, e = c // 2, c % 2
    nd = tsum(An[I, i, e] for I in cubes if i in I); Uc = canon([q[s] for s in sup[nd]])
    roots.append(dict(node=nd, targets=list(range(v)), coef=[Q(1, 3) if c in coords[t] else Q(-1, 6) for t in range(v)],
                      kind='centre', name='centre', frame=Uc))
nadd = sum(1 for a in args if a is not None)
print('active nodes', 'computed below')
print('p=%d h=%d v=%d: %d additions (%d signed), %d roots (%d side + %d centres)' % (
    p, h, v, nadd, sum(1 for s in sgn if s < 0), len(roots), len(roots) - h, h), flush=True)

# ---------------- frames ----------------
nn = len(args); succ = [[] for _ in range(nn)]
for x, a in enumerate(args):
    if a: succ[a[0]].append(x); succ[a[1]].append(x)
annv = [set() for _ in range(nn)]
for r in roots:
    annv[r['node']] |= set(perp(r['frame'])) if True else set()
for x in reversed(range(nn)):
    for y in succ[x]: annv[x] |= annv[y]
active = set(); st_ = [r['node'] for r in roots]
while st_:
    x = st_.pop()
    if x in active: continue
    active.add(x)
    if args[x]: st_.extend(args[x])
succ = [[y for y in sl if y in active] for sl in succ]
Phi = [None] * nn; badspan = 0
for x in range(nn):
    if x not in active: continue
    Phi[x] = perp(list(canon(list(annv[x]))))
    badspan += not inside([q[s] for s in sup[x]], Phi[x])
check('every node: span of its source addresses lies in its frame (orthogonal to all targets it reaches)', badspan == 0)
check('every centre frame is exactly its star span U_c (dim h-2), so the copy pays l = h(h-2) = %d' % (h * (h - 2)),
      all(Phi[r['node']] == r['frame'] and len(r['frame']) == h - 2 for r in roots if r['kind'] == 'centre'))
for x in range(nn):
    if args[x] and Phi[x] is not None:
        assert all(inside(Phi[y], Phi[x]) for y in args[x]), 'operand frame not nested'

# ---------------- cube-stack: frame descent (flag 'desc') ----------------
# Role chains are simulated exactly as the role loop below assigns them; every add node's frame is then moved between
# lo = sum of the frames preceding it on any chain (plus its source support) and hi = meet of the frames following it,
# minimising the local cover Lambda sum f(step), f(r) = r ln(3h/r). Centre-root nodes keep U_c. Input frames unused.
if 'desc' in FLAGS:
    import math as _m
    _m3 = 3 * h
    _f = lambda r: r * _m.log(_m3 / r) if r > 0 else 0.0
    _uses = defaultdict(list)
    for x in range(nn):
        if args[x] and Phi[x] is not None:
            for j, y in enumerate(args[x]): _uses[y].append(('n', x, j))
    for k, r in enumerate(roots): _uses[r['node']].append(('r', k))
    _ch = []; _ro = {}
    for x in range(nn):
        if Phi[x] is None: continue
        us = _uses[x]
        if args[x] is None:
            for u in us: _ro[x, u] = len(_ch); _ch.append([('F', Z0), ('F', canon([q[x]]))])
            continue
        a, b = args[x]; d_, c_ = _ro[a, ('n', x, 0)], _ro[b, ('n', x, 1)]
        _ch[d_].append(('N', x)); _ch[c_].append(('N', x)); _ro[x, us[0]] = d_
        for u in us[1:]: _ro[x, u] = len(_ch); _ch.append([('F', Z0), ('N', x)])
    fixed = set()
    for k, r in enumerate(roots):
        z = _ro[r['node'], ('r', k)]; _ch[z].append(('F', r['frame']))
        if r['kind'] == 'centre': fixed.add(r['node'])
    for c in _ch: c.append(('F', FULL))
    occ = defaultdict(list)
    for ci, c in enumerate(_ch):
        for i, e in enumerate(c):
            if e[0] == 'N': occ[e[1]].append((ci, i))
    U_ = {x: Phi[x] for x in occ}
    def fr(e): return e[1] if e[0] == 'F' else U_[e[1]]
    def local(x, dx):
        tot = 0.0
        for ci, i in occ[x]:
            p_, n_ = _ch[ci][i - 1], _ch[ci][i + 1]
            if not (p_[0] == 'N' and p_[1] == x): tot += _f(dx - len(fr(p_)))
            if not (n_[0] == 'N' and n_[1] == x): tot += _f(len(fr(n_)) - dx)
        return tot
    lam0 = sum(local(x, len(U_[x])) for x in occ)
    for sweep in range(6):
        moved = 0
        for x in (sorted(occ, reverse=True) if sweep % 2 == 0 else sorted(occ)):
            if x in fixed: continue
            lo = [q[s] for s in sup[x]]; hi = None
            for ci, i in occ[x]:
                p_, n_ = _ch[ci][i - 1], _ch[ci][i + 1]
                if not (p_[0] == 'N' and p_[1] == x): lo += list(fr(p_))
                if not (n_[0] == 'N' and n_[1] == x): hi = fr(n_) if hi is None else meet(hi, fr(n_))
            lo = canon(lo); hi = FULL if hi is None else hi
            if not inside(lo, hi): continue
            comp = []; cur = list(lo)
            for z in hi:
                if len(canon(cur + [z])) > len(canon(cur)): comp.append(z); cur.append(z)
            base = local(x, len(U_[x])); best = (0.0, None)
            for k in range(len(comp) + 1):
                d = local(x, len(lo) + k) - base
                if d < best[0] - 1e-9: best = (d, k)
            if best[1] is not None: U_[x] = canon(list(lo) + comp[:best[1]]); moved += 1
        if not moved: break
    lam1 = sum(local(x, len(U_[x])) for x in occ)
    ndiff = sum(1 for x in occ if U_[x] != Phi[x])
    for x in occ: Phi[x] = U_[x]
    print('frame descent: %d of %d node frames moved; local Lambda (double-counted per step) %.1f -> %.1f' % (ndiff, len(occ), lam0, lam1), flush=True)
    for x in range(nn):
        if args[x] and Phi[x] is not None:
            assert all(inside(Phi[y], Phi[x]) for y in args[x] if args[y]), 'operand frame not nested after descent'
            assert inside([q[s] for s in sup[x]], Phi[x]), 'support not in frame after descent'

# ---------------- roles and the forward word ----------------
# uses of node x, in a fixed order: consumers (addition nodes, operand position) then roots
uses = defaultdict(list)
for x in range(nn):
    if args[x] and x in active:
        for j, y in enumerate(args[x]): uses[y].append(('n', x, j))
for k, r in enumerate(roots): uses[r['node']].append(('r', k))
R = 0; role_of = {}; first_frame = {}; injections = []   # (role, source)
ops = []           # ('add', dst, ctl, coef, frame, node) / ('copy', dst, src, frame, node) in topological order
home = {}
for x in range(nn):
    if Phi[x] is None: continue
    us = uses[x]
    if args[x] is None:
        s = x                                   # input node index == source index
        for u in us:
            role_of[x, u] = R; injections.append((R, s)); first_frame[R] = canon([q[s]]); R += 1
        continue
    a, b = args[x]; dst, ctl = role_of[a, ('n', x, 0)], role_of[b, ('n', x, 1)]
    ops.append(('add', dst, ctl, Q(sgn[x]), Phi[x], x)); home[x] = dst
    role_of[x, us[0]] = dst
    for u in us[1:]:
        role_of[x, u] = R; first_frame[R] = Phi[x]; ops.append(('copy', R, dst, Phi[x], x)); R += 1
print('roles R = %d, ops %d, injections %d (%.1fs)' % (R, len(ops), len(injections), time.time() - t0), flush=True)
rootrole = [role_of[r['node'], ('r', k)] for k, r in enumerate(roots)]

# phase 1 = centre closure: every op that precedes a centre's completion on its roles (RAW/WAW/WAR closure)
prev = {}; pred = []
for i, op in enumerate(ops):
    rr = (op[1], op[2]); pred.append([prev[z] for z in rr if z in prev])
    for z in rr: prev[z] = i
lastw = {}
for i, op in enumerate(ops): lastw[op[1]] = i
phase1 = set(); st = [lastw[rootrole[k]] for k, r in enumerate(roots) if r['kind'] == 'centre' and rootrole[k] in lastw]
while st:
    i = st.pop()
    if i in phase1: continue
    phase1.add(i); st.extend(pred[i])
touched1 = set(z for i in phase1 for z in ops[i][1:3]) | set(r for r, s in injections)

# adjoint: coefficient of each role's OLD value in each target's reads (reverse sweep over the forward word)
FAST = 'fast' in FLAGS and MODE == 'scal'     # large p: centre-read part of the frame-0 prelude computed by a dirty-only run
co = [defaultdict(Q) for _ in range(R)]
for k, r in enumerate(roots):
    if FAST and r['kind'] == 'centre': continue    # untouched (gaugeable) roles never reach a centre read
    for t, c in zip(r['targets'], r['coef']): co[rootrole[k]][t] += c
for op in reversed(ops):
    if op[0] == 'add':
        _, d_, c_, cf, _, _ = op
        for t, val in list(co[d_].items()): co[c_][t] += cf * val
    else:
        _, d_, s_, _, _ = op
        for t, val in list(co[d_].items()): co[s_][t] += val
co = [{t: c for t, c in d_.items() if c} for d_ in co]

# ---------------- entrance gauges (optional) ----------------
sigma = {z: Z0 for z in range(R)}
if 'gauge' in FLAGS:
    limit = {}
    for k, r in enumerate(roots):
        if r['kind'] != 'side': continue
        for t in r['targets']:
            if t not in limit or len(r['frame']) < len(limit[t]): limit[t] = r['frame']
    cand = sorted((z for z in range(R) if z not in touched1 and co[z]), key=lambda z: (-len(first_frame[z]), len(co[z]), z))
    sel = []
    for z in cand:
        S_ = first_frame[z]
        for t in co[z]: S_ = meet(S_, limit[t])
        if len(S_) == 0: continue
        sigma[z] = S_; sel.append(z)
        for t in co[z]: limit[t] = S_
    print('gauges selected: %d, dims %s' % (len(sel), dict(Counter(len(sigma[z]) for z in sel))), flush=True)
    if 'g20' in FLAGS:                    # keep only the rank-(h-4) gauges (#144 keeps 4,840 of rank 20 at h = 24)
        for z in sel:
            if len(sigma[z]) != h - 4: sigma[z] = Z0
        sel = [z for z in sel if sigma[z]]
        print('kept only rank h-4 gauges: %d (per target cube %.2f)' % (len(sel), len(sel) / len(cubes)), flush=True)
    if 'omit' in FLAGS:
        om = sel[::3]
        for z in om: sigma[z] = Z0
        print('gauge omission: %d gauges moved to the frame-0 prelude' % len(om), flush=True)
gauged = [z for z in range(R) if sigma[z]]
OMITTED = om if ('gauge' in FLAGS and 'omit' in FLAGS) else []
NOPRE = set(OMITTED) if OMITTED else set(z for z in range(R) if z not in touched1 and co[z])

# ---------------- registers ----------------
if MODE == 'lit':
    bits = ((np.arange(N)[:, None] >> np.arange(h)[None, :]) & 1).astype(np.int64)
    class Reg:
        __slots__ = ('re', 'im', 'e')
        def __init__(s, re, im, e=0): s.re, s.im, s.e = re, im, e
        def copy(s): return Reg(s.re.copy(), s.im.copy(), s.e)
    def rand_reg(): return Reg(3 * np.array([rnd.randint(-7, 7) for _ in range(N)], dtype=np.int64),
                               3 * np.array([rnd.randint(-7, 7) for _ in range(N)], dtype=np.int64))
    def zero_reg(): return Reg(np.zeros(N, np.int64), np.zeros(N, np.int64))
    def norm(r):
        while r.e > 0 and not (r.re & 1).any() and not (r.im & 1).any(): r.re >>= 1; r.im >>= 1; r.e -= 1
        assert max(abs(r.re).max(), abs(r.im).max()) < (1 << 52), 'overflow guard'
    def align(a, b):
        if a.e < b.e: a.re = a.re << (b.e - a.e); a.im = a.im << (b.e - a.e); a.e = b.e
    def axpy(dst, sr, c):                 # dst += c * src, c rational with denominator 2^a 3^b
        s_ = sr.copy(); c = Q(c); num, den = c.numerator, c.denominator; e2 = 0
        while den % 2 == 0: den //= 2; e2 += 1
        assert den in (1, 3)
        re, im = s_.re * num, s_.im * num
        if den == 3:
            assert not (re % 3).any() and not (im % 3).any(), 'non-integral /3'
            re //= 3; im //= 3
        s_ = Reg(re, im, s_.e + e2)
        if dst.e < s_.e: align(dst, s_)
        else: align(s_, dst)
        dst.re = dst.re + s_.re; dst.im = dst.im + s_.im; norm(dst)
    def eqreg(a, b):
        a, b = a.copy(), b.copy()
        if a.e < b.e: align(a, b)
        else: align(b, a)
        return bool((a.re == b.re).all() and (a.im == b.im).all())
    def sx(r, k, inv):
        sh = (N >> (k + 1), 2, 1 << k); Rr = r.re.reshape(sh); M = r.im.reshape(sh)
        ar, br, ai, bi = Rr[:, 0, :].copy(), Rr[:, 1, :].copy(), M[:, 0, :].copy(), M[:, 1, :].copy()
        if not inv: ur, ui, vr, vi, fr, fi = ar + bi, ai - br, br + ai, bi - ar, 1, 1
        else: ur, ui, vr, vi, fr, fi = ar - bi, ai + br, br - ai, bi + ar, 1, -1
        Rr[:, 0, :] = fr * ur - fi * ui; M[:, 0, :] = fr * ui + fi * ur
        Rr[:, 1, :] = fr * vr - fi * vi; M[:, 1, :] = fr * vi + fi * vr
        r.re = Rr.reshape(N); r.im = M.reshape(N); r.e += 1
    def rot(r, k):                        # multiply entrywise by i^k (k array mod 4)
        k = k % 4; c = np.array([1, 0, -1, 0])[k]; s_ = np.array([0, 1, 0, -1])[k]
        r.re, r.im = c * r.re - s_ * r.im, s_ * r.re + c * r.im
    def perm(r, idx_):
        re = np.empty(N, np.int64); im = np.empty(N, np.int64); re[idx_] = r.re; im[idx_] = r.im; r.re, r.im = re, im
    def f2inv(cols):                      # inverse of the h x h F2 matrix given by column ints
        n_ = len(cols); rows = [[(cols[j] >> i) & 1 for j in range(n_)] + [int(i == k) for k in range(n_)] for i in range(n_)]
        for c in range(n_):
            piv = next(i for i in range(c, n_) if rows[i][c]); rows[c], rows[piv] = rows[piv], rows[c]
            for i in range(n_):
                if i != c and rows[i][c]: rows[i] = [x ^ y for x, y in zip(rows[i], rows[c])]
        inv_rows = [r_[n_:] for r_ in rows]
        return [sum(inv_rows[i][j] << i for i in range(n_)) for j in range(n_)]
    def apply_cols(cols, a):              # vectorized G a for all a
        out = np.zeros(N, np.int64)
        for k, c in enumerate(cols): out ^= bits[:, k] * c
        return out
    FR = {}
    def frame_ops(U, phase_twist=False):
        key = (U, phase_twist)
        if key in FR: return FR[key]
        cols = list(U)
        for j in range(h):
            if len(cols) == h: break
            if not inside([1 << j], cols): cols.append(1 << j)
        assert len(red(cols)) == h
        Gi = f2inv(cols)
        # M = I + (G G^T)^{-1}
        GT = [sum(((cols[j] >> i) & 1) << j for j in range(h)) for i in range(h)]       # columns of G^T
        GGt = [apply_cols_int(cols, c) for c in GT]
        Mc = f2inv(GGt); Mc = [c ^ (1 << j) for j, c in enumerate(Mc)]
        Mm = np.array([[(Mc[j] >> i) & 1 for j in range(h)] for i in range(h)], dtype=np.int64)
        assert (Mm == Mm.T).all()
        up = np.triu(Mm, 1)
        qv = (bits @ np.diag(Mm) + 2 * ((bits @ up) * bits).sum(1)) % 4
        if phase_twist: qv = (qv + 2 * bits[:, 0]) % 4          # control: same label, different rank-zero factor (Z_0)
        FR[key] = (apply_cols(cols, None), apply_cols(Gi, None), qv, len(U))
        return FR[key]
    def apply_cols_int(cols, x):
        y = 0; k = 0
        while x:
            if x & 1: y ^= cols[k]
            x >>= 1; k += 1
        return y
    def T(r, U, inv=False, twist=False):  # r <- T_U r  (or T_U^{-1} r)
        if len(U) == 0: return
        Gidx, Giidx, qv, rr = frame_ops(U, twist)
        rot(r, -qv); perm(r, Giidx)
        for k in range(rr): sx(r, k, inv)
        perm(r, Gidx); rot(r, qv); norm(r)
    def colrank(U, V):                    # measured Fourier rank of T_V T_U^{-1}: log2 of column support, 3 columns
        out = set()
        for _ in range(3):
            r = zero_reg(); r.re[rnd.randrange(N)] = 1; T(r, U, inv=True); T(r, V)
            nz = int(((r.re != 0) | (r.im != 0)).sum()); t = nz.bit_length() - 1; out.add(t if nz == 1 << t else -1)
        return out
else:
    K_ = 3
    class Reg:
        __slots__ = ('val',)
        def __init__(s, val): s.val = val
        def copy(s): return Reg(list(s.val))
    def rand_reg(): return Reg([Q(rnd.randint(-50, 50), rnd.choice([1, 2, 3, 4])) for _ in range(K_)])
    def zero_reg(): return Reg([Q(0)] * K_)
    def axpy(dst, sr, c): dst.val = [a + Q(c) * b for a, b in zip(dst.val, sr.val)]
    def eqreg(a, b): return a.val == b.val
    def T(r, U, inv=False, twist=False): pass
regs = {}; frame = {}; hist = Counter(); badgate = [0]; badnest = [0]; measured = defaultdict(set)
def move(name, V, kind, twist=False):
    U = frame[name]
    if U == V and not twist: return
    if not (inside(U, V) or (kind == 'copy' and inside(V, U))): badnest[0] += 1
    d_ = dL(U, V); hist[kind, d_] += 1
    if MODE == 'lit':
        if len(measured[kind, d_]) < 2 and len(measured) < 60: measured[kind, d_] |= colrank(U, V)
        T(regs[name], U, inv=True); T(regs[name], V, twist=twist)
    frame[name] = V
def gate(dst, sr, c):
    if frame[dst] != frame[sr]: badgate[0] += 1
    axpy(regs[dst], regs[sr], c)

# initial state: sources at T_<q>, targets at frame 0, aux roles arbitrary at their gauge frame
x0 = {}; y0 = {}; w0 = {}
for s in range(v):
    xl = rand_reg(); x0[s] = xl.copy(); r_ = xl.copy(); T(r_, canon([q[s]])); regs['x', s] = r_; frame['x', s] = canon([q[s]])
for t in range(v):
    y0[t] = rand_reg(); regs['y', t] = y0[t].copy(); frame['y', t] = Z0
for z in range(R):
    w0[z] = rand_reg(); regs['z', z] = w0[z].copy(); frame['z', z] = sigma[z]
print('registers ready (%.1fs)' % (time.time() - t0), flush=True)

# ---- t = 0: old-value reads of every gauge-0 role at frame 0 (exact adjoint coefficients) ----
dropped = None
if FAST:
    # y_T -= (J M z_old)_T restricted to gauge-0 roles: exact by linearity (forward word on the old values alone).
    # Every gauge-0 role and every target is at frame 0 here, so each such read is a legal frame-0 gate.
    assert all(frame['y', t] == Z0 for t in range(v)) and all(frame['z', z] == Z0 for z in range(R) if not sigma[z])
    val = {z: (regs['z', z].copy() if not sigma[z] else zero_reg()) for z in range(R)}
    for op in ops:
        if op[0] == 'add': axpy(val[op[1]], val[op[2]], op[3])
        else: axpy(val[op[1]], val[op[2]], 1)
    corr = {t: zero_reg() for t in range(v)}
    for k, r in enumerate(roots):
        for t, c in zip(r['targets'], r['coef']): axpy(corr[t], val[rootrole[k]], c)
    for t in range(v): axpy(regs['y', t], corr[t], -1)
    del val, corr
for z in range(R):
    if sigma[z] or FAST: continue
    for t, c in co[z].items():
        if CTRL == 'noprelude' and (dropped is None or dropped == z) and z in NOPRE: dropped = z; continue
        gate(('y', t), ('z', z), -c)
# ---- injections at the sources' rank-one entrance frames ----
for z, s in injections:
    move(('z', z), canon([q[s]]), 'aux'); gate(('z', z), ('x', s), 1)
EXEC = []                                 # executed op order: with shared registers the cleanup must be its exact reverse
def do_op(op):
    EXEC.append(op)
    if op[0] == 'add':
        _, d_, c_, cf, F_, _ = op; move(('z', d_), F_, 'aux'); move(('z', c_), F_, 'aux'); gate(('z', d_), ('z', c_), cf)
    else:
        _, d_, s_, F_, _ = op; move(('z', d_), F_, 'aux'); move(('z', s_), F_, 'aux'); gate(('z', d_), ('z', s_), 1)
# ---- phase 1: centre closure, copied centres at frame 0 ----
for i, op in enumerate(ops):
    if i in phase1: do_op(op)
ell = 0
for k, r in enumerate(roots):
    if r['kind'] != 'centre': continue
    z = rootrole[k]; move(('z', z), r['frame'], 'aux')
    if CTRL == 'skipclone':
        continue
    regs['c', k] = regs['z', z].copy(); frame['c', k] = frame['z', z]
    move(('c', k), Z0, 'copy'); ell += len(r['frame'])
    for t, c in zip(r['targets'], r['coef']): gate(('y', t), ('c', k), c)
    del regs['c', k]
# ---- phase 2 events by frame dimension ----
ev = []
firstop = {}
for i, op in enumerate(ops):
    for z in op[1:3]: firstop.setdefault(z, i)
for z in gauged:
    ev.append((len(sigma[z]), 0, firstop.get(z, 0), 'old', z))
for i, op in enumerate(ops):
    if i not in phase1: ev.append((len(op[4] if op[0] == 'add' else op[3]), 1, i, 'op', i))
for k, r in enumerate(roots):
    if r['kind'] == 'side': ev.append((len(r['frame']), 2, k, 'read', k))
ev.sort()
if CTRL == 'lateread':                    # one gauged old-value read moved after that role's first write
    j = next(i for i, e in enumerate(ev) if e[3] == 'old' and e[4] in firstop); e = ev.pop(j)
    k_ = next(i for i, e2 in enumerate(ev) if e2[3] == 'op' and e2[4] == firstop[e[4]]); ev.insert(k_ + 1, e)
# ---- cube-birth: birth-read reuse (flag 'birth') ----
# A gauged role B is untouched until its old-value read at sigma_B, which cancels whatever its register holds; so B
# may be born on the register of a dead donor A (ungauged, not a root role, last op executed before B's read, last
# frame E_A inside sigma_B). The donor register moves E_A -> sigma_B LITERALLY and becomes B's register.
BIRTH = {}; DONOR = {}
if 'birth' in FLAGS:
    P1 = sorted(phase1); seq = [('op', i) for i in P1] + [(e[3], e[4]) for e in ev]
    lastp, lastf = {}, {}
    for kk, (kind, i) in enumerate(seq):
        if kind != 'op': continue
        op = ops[i]; F_ = op[4] if op[0] == 'add' else op[3]
        for z in op[1:3]: lastp[z] = kk; lastf[z] = F_
    oldp = {i: kk for kk, (kind, i) in enumerate(seq) if kind == 'old'}
    roots_set = set(rootrole); gset = set(gauged)
    donors = sorted((z for z in lastp if z not in roots_set and z not in gset), key=lambda z: lastp[z])
    used = set()
    for z in sorted(gauged, key=lambda z: oldp[z]):
        best = None
        for A in donors:
            if lastp[A] >= oldp[z]: break
            if A in used: continue
            ok = inside(lastf[A], sigma[z])
            if CTRL == 'birthframe' and not BIRTH: ok = not ok          # control: one hand-off into a sigma missing E_A
            if ok and (best is None or len(lastf[A]) > len(lastf[best])): best = A
        if best is not None: BIRTH[z] = best; DONOR[best] = z; used.add(best)
    _t = sum(1 for z in gauged for A in donors if lastp[A] < oldp[z]); _f = sum(1 for z in gauged for A in donors if inside(lastf[A], sigma[z]))
    print('diag: donors %d, timely (z,A) %d, frame-compatible (z,A) %d, sigma dims %s, donor last dims %s' % (len(donors), _t, _f,
          dict(Counter(len(sigma[z]) for z in gauged)), dict(Counter(len(lastf[A]) for A in donors))), flush=True)
    print('birth pairs: %d of %d gauged (donor last-frame dims %s)' % (len(BIRTH), len(gauged),
          dict(Counter(len(lastf[A]) for A in DONOR))), flush=True)
    if CTRL == 'birthlate':               # one recipient's old-value read moved after its first op (read after first touch)
        zl = next(z for z in sorted(BIRTH, key=lambda z: oldp[z]) if z in firstop)
        j = next(i for i, e in enumerate(ev) if e[3] == 'old' and e[4] == zl); e = ev.pop(j)
        k_ = next(i for i, e2 in enumerate(ev) if e2[3] == 'op' and e2[4] == firstop[zl]); ev.insert(k_ + 1, e)
DROP = next(iter(sorted(BIRTH)), None) if CTRL == 'dropbirth' else None
for e in ev:
    if e[3] == 'old':
        z = e[4]
        if z in BIRTH:
            A = BIRTH[z]; move(('z', A), sigma[z], 'aux')                 # hand-off E_A -> sigma_B, literal
            regs['z', z] = regs['z', A]; frame['z', z] = frame['z', A]       # B lives on A's register from here
            if z == DROP: continue                                           # control: the birth read is dropped
        for t, c in co[z].items():
            move(('y', t), sigma[z], 'tgt'); gate(('y', t), ('z', z), -c)
    elif e[3] == 'op': do_op(ops[e[4]])
    else:
        k = e[4]; r = roots[k]; z = rootrole[k]; move(('z', z), r['frame'], 'aux')
        for t, c in zip(r['targets'], r['coef']): move(('y', t), r['frame'], 'tgt'); gate(('y', t), ('z', z), c)
# ---- K on the original source registers ----
Kmat = {}
for c_, I in enumerate(cubes):
    for pa in (0, 1):
        cls = [8 * c_ + k for k in range(8) if par[8 * c_ + k] == pa]; opp = [8 * c_ + k for k in range(8) if par[8 * c_ + k] != pa]
        if CTRL == 'badgroup' and pa == 0:  # WRONG grouping: parity-0 half-cubes of cubes 2j and 2j+1 exchange targets
            c2 = c_ + 1 if c_ % 2 == 0 else c_ - 1
            if c2 < len(cubes): opp = [8 * c2 + k for k in range(8) if par[8 * c2 + k] != pa]
        U = canon([q[s] for s in cls])
        # K_{T,S} = -(B_block)_{T,S} for T != S: (P - A)/2 inside the cube
        def Kts(T_, S_):
            dd = sum(a != b for a, b in zip(ports[T_][1], ports[S_][1]))
            return Q(1, 2) if dd == 3 else Q(-1, 2) if dd == 1 else Q(0)
        Mblk = [[Kts(T_, S_) for S_ in cls] for T_ in opp]      # register i (holding source cls[i]) <- row opp[i]
        for i_, s in enumerate(cls):
            move(('x', s), U, 'src', twist=(CTRL == 'rep' and i_ == 0 and c_ == 0))
        new = []
        for i_ in range(4):
            r_ = zero_reg(); frame['tmp'] = U; regs['tmp'] = r_
            for j_, s in enumerate(cls): gate('tmp', ('x', s), Mblk[i_][j_])
            new.append(regs.pop('tmp'))
        for i_, s in enumerate(cls): regs['x', s] = new[i_]
        Kmat[c_, pa] = (cls, opp, Mblk)
for (c_, pa), (cls, opp, Mblk) in Kmat.items():
    for i_, s in enumerate(cls):
        t = opp[i_]; qp = perp([q[t]])
        move(('x', s), qp, 'src'); move(('y', t), qp, 'tgt'); gate(('y', t), ('x', s), 1)
for (c_, pa), (cls, opp, Mblk) in Kmat.items():
    for s in cls: move(('x', s), FULL, 'src')
    if CTRL == 'noKundo': continue
    Mi = [[Q(0)] * 4 for _ in range(4)]
    A_ = [list(r_) + [Q(int(i == j)) for j in range(4)] for i, r_ in enumerate(Mblk)]
    for c in range(4):
        pv = next(i for i in range(c, 4) if A_[i][c] != 0); A_[c], A_[pv] = A_[pv], A_[c]
        A_[c] = [x / A_[c][c] for x in A_[c]]
        for i in range(4):
            if i != c and A_[i][c] != 0: A_[i] = [x - A_[i][c] * y for x, y in zip(A_[i], A_[c])]
    Minv = [r_[4:] for r_ in A_]
    new = []
    for i_ in range(4):
        r_ = zero_reg(); frame['tmp'] = FULL; regs['tmp'] = r_
        for j_, s in enumerate(cls): gate('tmp', ('x', s), Minv[i_][j_])
        new.append(regs.pop('tmp'))
    for i_, s in enumerate(cls): regs['x', s] = new[i_]
# ---- every role to the full frame, inverse word at F, sources subtracted at F ----
for z in range(R):
    if z in DONOR: continue
    move(('z', z), FULL, 'aux')
for A in DONOR: frame['z', A] = FULL
for op in reversed(EXEC if 'birth' in FLAGS else ops):
    if op[0] == 'add': gate(('z', op[1]), ('z', op[2]), -op[3])
    else: gate(('z', op[1]), ('z', op[2]), -1)
for z, s in injections: gate(('z', z), ('x', s), -1)
print('word done (%.1fs)' % (time.time() - t0), flush=True)

# ---------------- checks ----------------
check('every gate between registers at the IDENTICAL frame (%d violations)' % badgate[0], badgate[0] == 0)
check('every transition nested and ascending (copied-centre copies descend) (%d violations)' % badnest[0], badnest[0] == 0)
bad_t = 0
for t in range(v):
    want = y0[t].copy(); axpy(want, x0[t], 1); T(want, perp([q[t]]))
    bad_t += not (frame['y', t] == perp([q[t]]) and eqreg(regs['y', t], want))
check('targets end at T_{q_T^perp}(y + x_T): y += (K+H+B)x = x  (%d/%d wrong)' % (bad_t, v), bad_t == 0)
bad_x = 0
for s in range(v):
    want = x0[s].copy(); T(want, FULL); bad_x += not eqreg(regs['x', s], want)
check('sources end at F x_S (K undone at the full frame) (%d/%d wrong)' % (bad_x, v), bad_x == 0)
bad_z = 0
for z in range(R):
    if z in BIRTH: continue
    want = w0[z].copy(); T(want, sigma[z], inv=True); T(want, FULL); bad_z += not eqreg(regs['z', z], want)
check('aux roles restored: end at F T_sigma^{-1} w0 for arbitrary dirty w0 (%d/%d wrong)' % (bad_z, R), bad_z == 0)
src_h = Counter(d_ for (k, d_), c in hist.items() if k == 'src' for _ in range(c))
print('source transitions (rank: count):', dict(sorted(Counter({d_: c for (k, d_), c in hist.items() if k == 'src'}).items())))
print('target transitions (rank: count):', dict(sorted(Counter({d_: c for (k, d_), c in hist.items() if k == 'tgt'}).items())))
print('copied-centre transforms: l =', ell, ' (h(h-2) =', h * (h - 2), ')')
aux_mass = sum(d_ * c for (k, d_), c in hist.items() if k == 'aux')
want_mass = R * h - sum(len(sigma[z]) for z in range(R)) - sum(h - len(sigma[z]) for z in BIRTH)
check('aux ledger: sum of aux transition ranks = R h - sum dim sigma - sum_pairs (h - dim sigma_B)  (%d vs %d), slots %d' % (aux_mass, want_mass, R - len(BIRTH)),
      aux_mass == want_mass)
srcc = Counter({d_: c for (k, d_), c in hist.items() if k == 'src'})
check('source histogram [2, h-4, 1] per port: %s' % dict(srcc), srcc == Counter({2: v, h - 4: v, 1: v}) or CTRL)
if MODE == 'lit':
    okm = all(ms == {d_} for (k, d_), ms in measured.items())
    check('measured column-support Fourier ranks == d(L_U,L_V) on every sampled transition type %s' %
          sorted((k, d_, sorted(ms)) for (k, d_), ms in measured.items()), okm)
import math as _mm
print('aux Lambda (sum f over aux transitions incl. copies) = %.2f' % sum(c * d * _mm.log(3 * h / d) for (k_, d), c in hist.items() if k_ in ('aux', 'copy') and d > 0), flush=True)
print('RESULT p=%d %s seed=%d flags=%s control=%s: %s (%.1fs)' % (p, MODE, SEED, sorted(FLAGS), CTRL or '-',
      'ALL PASS' if all(OUT.values()) else 'SOME FAIL: ' + '; '.join(k for k, o in OUT.items() if not o), time.time() - t0), flush=True)
