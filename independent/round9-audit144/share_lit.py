"""Item 3, own code: #144's SEQUENTIAL CROSS-STAGE reuse of one auxiliary bank in three orthogonal h-blocks, literal and
exact over Z[i][1/2] on all 2^m addresses, m = 3h (pattern of audit-137/work/sim.py, frames from the note's
arbitrary-subspace interface so degenerate gauges occur).

Ambient F2^m, an INDEPENDENT orthogonal partition P_1 + P_2 + P_3 = Omega(E_1) + Omega(E_2) + Omega(E_3) (E_j coordinate
h-blocks, Omega a random orthogonal map), H_j: A = E_1 -> P_j orthogonal, a random vertex g; the role (g, r) is used in stage
j by the core at vertex k_j = g H_j, whose active block is k_j A = g P_j. Stage offsets O_1 = 0, O_2, O_3 = random subspaces
of the inactive coordinates (transported by k_j), so offset cancellation is computed, not assumed.
Frames at a vertex: T = P_k (T^A_U (x) T^rest_O) P_k^{-1}, T^A_U = K_G C_{E_r} K_G^{-1} on the active coordinate block
(K_G = D_M P_G, M = I + (G G^T)^{-1}, C = sqrt X), i.e. block-factored representatives.
Local chain U_0 < U_1 < ... < U_{h-1} = A from a basis G e_0 = q (odd), G e_1 = w (even, w.q = 0: U_1 is DEGENERATE), rest
random. Toy transparent core (audit-137's): role r joins at step j_r (gauge sigma_r = U_{j_r}; j_r = -1: gauge 0):
dirty pre-read y -= d z, z += c x, all climb together, y += d z, z -= c x at A. Stage 1 and 3 forward (Y += X), stage 2
the literal inverse word with roles exchanged (reversed: its aux residual is the inverse).
Checks: each stage's data on the SHARED bank == the same core on private random scratch; bank == U_3 U_2^{-1} U_1 a0 with
U_j = T_{last} T_{first}^{-1} of that core; tail = F (U_3 U_2^{-1} U_1)^{-1} restores F a0 exactly; tail column support
2^(3 dim sigma) with unit*((1+i)/2)^t entries (one width-3 dim sigma child).
Controls (must FAIL): lockstep (round-robin interleave of the three cores on the bank), shorttail (tail without the stage-3
residual), overlap (P_3 = P_1: blocks not orthogonal -> tail rank != 3 dim sigma), global (non-factored representatives
with the same labels: reported, tail rank).
Usage: python3 share_lit.py h SEED [control]"""
import sys, random, itertools
import numpy as np
h = int(sys.argv[1]); SEED = int(sys.argv[2]); CTRL = sys.argv[3] if len(sys.argv) > 3 else ''
rnd = random.Random(SEED * 7919 + h); m = 3 * h; N = 1 << m
OUT = {}
def check(name, ok): OUT[name] = bool(ok); print(('PASS ' if ok else 'FAIL ') + name, flush=True)
bits = ((np.arange(N)[:, None] >> np.arange(m)[None, :]) & 1).astype(np.int64)
dot = lambda a, b: bin(a & b).count('1') & 1
e = lambda i: 1 << i
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
def applyM(cols, x):
    y = 0; k = 0
    while x:
        if x & 1: y ^= cols[k]
        x >>= 1; k += 1
    return y
def compose(A, B): return [applyM(A, c) for c in B]
def f2inv(cols):
    n_ = len(cols); rows = [[(cols[j] >> i) & 1 for j in range(n_)] + [int(i == k) for k in range(n_)] for i in range(n_)]
    for c in range(n_):
        piv = next(i for i in range(c, n_) if rows[i][c]); rows[c], rows[piv] = rows[piv], rows[c]
        for i in range(n_):
            if i != c and rows[i][c]: rows[i] = [x ^ y for x, y in zip(rows[i], rows[c])]
    return [sum(rows[i][n_ + j] << i for i in range(n_)) for j in range(n_)]
I_m = [e(i) for i in range(m)]
def is_orth(M): return all(dot(M[i], M[j]) == (i == j) for i in range(m) for j in range(m))
def rand_orth(steps=60):
    g = list(I_m)
    for _ in range(steps):
        if rnd.random() < 0.5:
            p_ = list(range(m)); rnd.shuffle(p_); g = compose([e(p_[k]) for k in range(m)], g)
        else:
            while True:
                vv = rnd.randrange(1, N)
                if bin(vv).count('1') % 2 == 0: break
            g = compose([x ^ (vv if dot(x, vv) else 0) for x in I_m], g)
    assert is_orth(g); return g

# ---------------- exact registers ----------------
class Reg:
    __slots__ = ('re', 'im', 'e')
    def __init__(s, re, im, e_=0): s.re, s.im, s.e = re, im, e_
    def copy(s): return Reg(s.re.copy(), s.im.copy(), s.e)
def rand_reg(): return Reg(np.array([rnd.randint(-9, 9) for _ in range(N)], np.int64), np.array([rnd.randint(-9, 9) for _ in range(N)], np.int64))
def zero_reg(): return Reg(np.zeros(N, np.int64), np.zeros(N, np.int64))
def norm(r):
    while r.e > 0 and not (r.re & 1).any() and not (r.im & 1).any(): r.re >>= 1; r.im >>= 1; r.e -= 1
    assert max(abs(r.re).max(), abs(r.im).max()) < (1 << 52)
def align(a, b):
    if a.e < b.e: a.re = a.re << (b.e - a.e); a.im = a.im << (b.e - a.e); a.e = b.e
def eqreg(a, b):
    a, b = a.copy(), b.copy()
    if a.e < b.e: align(a, b)
    else: align(b, a)
    return bool((a.re == b.re).all() and (a.im == b.im).all())
def axpy(dst, sr, c):        # dst += c*src, c Gaussian integer
    s_ = sr.copy()
    if dst.e < s_.e: align(dst, s_)
    else: align(s_, dst)
    cr, ci = c; dst.re = dst.re + cr * s_.re - ci * s_.im; dst.im = dst.im + cr * s_.im + ci * s_.re; norm(dst)
def sx(r, k, inv):
    sh = (N >> (k + 1), 2, 1 << k); Rr = r.re.reshape(sh); M = r.im.reshape(sh)
    ar, br, ai, bi = Rr[:, 0, :].copy(), Rr[:, 1, :].copy(), M[:, 0, :].copy(), M[:, 1, :].copy()
    if not inv: ur, ui, vr, vi, fr, fi = ar + bi, ai - br, br + ai, bi - ar, 1, 1
    else: ur, ui, vr, vi, fr, fi = ar - bi, ai + br, br - ai, bi + ar, 1, -1
    Rr[:, 0, :] = fr * ur - fi * ui; M[:, 0, :] = fr * ui + fi * ur
    Rr[:, 1, :] = fr * vr - fi * vi; M[:, 1, :] = fr * vi + fi * vr
    r.re = Rr.reshape(N); r.im = M.reshape(N); r.e += 1
def rot(r, k):
    k = k % 4; c = np.array([1, 0, -1, 0])[k]; s_ = np.array([0, 1, 0, -1])[k]
    r.re, r.im = c * r.re - s_ * r.im, s_ * r.re + c * r.im
def perm(r, idx_):
    re = np.empty(N, np.int64); im = np.empty(N, np.int64); re[idx_] = r.re; im[idx_] = r.im; r.re, r.im = re, im
def vec_apply(cols):
    out = np.zeros(N, np.int64)
    for k, c in enumerate(cols): out ^= bits[:, k] * c
    return out
_cache = {}
def Kops(cols):
    """(P_G idx, P_G^-1 idx, q_M phase array) for the full m x m matrix G (column ints)."""
    key = tuple(cols)
    if key in _cache: return _cache[key]
    Gi = f2inv(cols)
    GT = [sum(((cols[j] >> i) & 1) << j for j in range(m)) for i in range(m)]
    Mc = f2inv([applyM(cols, c) for c in GT]); Mc = [c ^ (1 << j) for j, c in enumerate(Mc)]
    Mm = np.array([[(Mc[j] >> i) & 1 for j in range(m)] for i in range(m)], dtype=np.int64); assert (Mm == Mm.T).all()
    qv = (bits @ np.diag(Mm) + 2 * ((bits @ np.triu(Mm, 1)) * bits).sum(1)) % 4
    _cache[key] = (vec_apply(cols), vec_apply(Gi), qv); return _cache[key]
# An operator = list of primitive ops ('perm', idx) | ('rot', q) | ('sx', k, inv); applied in order.
def op_apply(ops, r):
    for o in ops:
        if o[0] == 'perm': perm(r, o[1])
        elif o[0] == 'rot': rot(r, o[1])
        else: sx(r, o[1], o[2]); norm(r)       # normalise after every factor (int64 headroom)
    norm(r)
def op_inv(ops):
    out = []
    for o in reversed(ops):
        if o[0] == 'perm':
            inv = np.empty(N, np.int64); inv[o[1]] = np.arange(N); out.append(('perm', inv))
        elif o[0] == 'rot': out.append(('rot', -o[1]))
        else: out.append(('sx', o[1], not o[2]))
    return out
def KCK(cols, coords):
    """K_G C_coords K_G^{-1} as primitive ops (applied right to left: D^-1, P_G^-1, C, P_G, D)."""
    Gi_, Gii, qv = Kops(cols)
    return [('rot', -qv), ('perm', Gii)] + [('sx', c, False) for c in coords] + [('perm', Gi_), ('rot', qv)]

# ---------------- geometry ----------------
A = list(range(h))                         # local active coordinate block
def complete(basis, coords):               # complete a basis inside the coordinate set to a basis of it
    out = list(basis)
    for c in coords:
        if len(out) == len(coords): break
        if len(red(out + [e(c)])) > len(out): out.append(e(c))
    return out
# local chain basis: q odd, w even with w.q = 0 (U_1 degenerate), rest random
if h == 3: qv_, wv_ = e(0), e(1) | e(2)
else: qv_, wv_ = e(0) | e(1) | e(2), e(0) | e(1)
assert dot(qv_, qv_) == 1 and dot(wv_, wv_) == 0 and dot(qv_, wv_) == 0
while True:
    Gloc = [qv_, wv_] + [rnd.randrange(1, 1 << h) for _ in range(h - 2)]
    if len(red(Gloc)) == h: break
def gram_rank(vs):
    rows = [sum(dot(a, b) << j for j, b in enumerate(vs)) for a in vs]; return len(red(rows))
print('local chain dims/gram ranks:', [(j + 1, gram_rank(Gloc[:j + 1])) for j in range(h)], flush=True)
JOINS = [-1, 0, 1, h - 1]                  # gauge dims 0, 1, 2 (degenerate), h
NR = len(JOINS)
_gcache = {}
def frame_ops(k, j, offset, glob=False):
    """frame at vertex k (orthogonal m x m), local chain step j (U_j = span Gloc[:j+1]; j=-1: 0), stage offset
    (list of local inactive vectors). Block-factored: G = Gloc (+) completion on A, offset basis completed on the rest."""
    rest = list(range(h, m))
    GA = complete(Gloc[:j + 1], A) if j >= 0 else [e(c) for c in A]
    Go = complete(offset, rest) if offset else [e(c) for c in rest]
    cols = GA + Go
    coords = list(range(j + 1)) + list(range(h, h + len(offset)))
    if glob and (j, tuple(offset)) in _gcache: cols, coords = _gcache[j, tuple(offset)]
    elif glob:                             # same subspace, generic NON-factored completion, fixed per subspace (control)
        U = (Gloc[:j + 1] if j >= 0 else []) + list(offset)
        cols = complete(U, list(range(m)))
        cols = U + [c for c in cols[len(U):]]
        # scramble the completion with random combinations of all vectors (still a basis, first len(U) = U)
        tail = cols[len(U):]
        for _ in range(3 * m):
            a, b = rnd.randrange(len(tail)), rnd.randrange(len(cols))
            if b >= len(U) and b - len(U) == a: continue
            nt = tail[a] ^ cols[b]
            if len(red(cols[:len(U)] + tail[:a] + [nt] + tail[a + 1:])) == m: tail[a] = nt
        cols = cols[:len(U)] + tail
        coords = list(range(len(U))); _gcache[j, tuple(offset)] = (cols, coords)
    kidx = vec_apply(k); kinv = np.empty(N, np.int64); kinv[kidx] = np.arange(N)
    return [('perm', kinv)] + KCK(cols, coords) + [('perm', kidx)]
def core_word(k, offset, coefs, srcn, dstn, glob=False):
    Fr = lambda j: frame_ops(k, j, offset, glob)
    trans = lambda j: op_inv(Fr(j - 1)) + Fr(j)
    W = []
    for r, j in enumerate(JOINS):
        if j == -1: W.append(('T', ('z', r), trans(0)))
    active = []
    for step in range(h):
        if step > 0:
            for reg in [srcn, dstn] + [('z', r) for r in active]: W.append(('T', reg, trans(step)))
        for r, j in enumerate(JOINS):
            if max(j, 0) == step:
                c, d = coefs[r]
                W.append(('ax', dstn, ('z', r), (-d[0], -d[1]))); W.append(('ax', ('z', r), srcn, c)); active.append(r)
    for r in active:
        c, d = coefs[r]; W.append(('ax', dstn, ('z', r), d)); W.append(('ax', ('z', r), srcn, (-c[0], -c[1])))
    return W
def word_inverse(W):
    out = []
    for w in reversed(W):
        out.append(('T', w[1], op_inv(w[2])) if w[0] == 'T' else ('ax', w[1], w[2], (-w[3][0], -w[3][1])))
    return out
def run(word, regs, aux):
    for w in word:
        Rg = (lambda n: aux[n[1]] if n[0] == 'z' else regs[n])
        if w[0] == 'T': op_apply(w[2], Rg(w[1]))
        else: axpy(Rg(w[1]), Rg(w[2]), w[3])
def residual(word, r): return [op for w in word if w[0] == 'T' and w[1] == ('z', r) for op in w[2]]

# ---------------- cross-stage sharing ----------------
Omega = rand_orth(); g = rand_orth()
blockperm = lambda j: [e(((c // h + j) % 3) * h + c % h) for c in range(m)]      # E_1 -> E_{1+j} (coordinate shift)
Hj = {s: compose(Omega, blockperm(s - 1)) for s in (1, 2, 3)}
if CTRL == 'overlap': Hj[3] = Hj[1]                                               # P_3 = P_1: blocks overlap
k = {s: compose(g, Hj[s]) for s in (1, 2, 3)}
Pblk = {s: [applyM(k[s], e(c)) for c in A] for s in (1, 2, 3)}
if CTRL != 'overlap':
    assert all(dot(a, b) == 0 for s in (1, 2, 3) for t in (1, 2, 3) if s < t for a in Pblk[s] for b in Pblk[t])
    assert len(red(sum(Pblk.values(), []))) == m
offs = {1: [], 2: [e(h), e(h + 1) ^ e(2 * h)], 3: [e(2 * h - 1) ^ e(m - 1), e(h + 2) ^ e(h) if h > 3 else e(h + 2)]}
for s in (2, 3): assert len(red(offs[s])) == len(offs[s])
coefs = [((rnd.randint(1, 3), rnd.randint(-2, 2)), (rnd.randint(1, 3), rnd.randint(-2, 2))) for _ in range(NR)]
GLOB = CTRL == 'global'
words = {1: core_word(k[1], offs[1], coefs, 'X', 'Y', GLOB),
         2: word_inverse(core_word(k[2], offs[2], coefs, 'Y', 'X', GLOB)),
         3: core_word(k[3], offs[3], coefs, 'X', 'Y', GLOB)}
din = {s: {'X': rand_reg(), 'Y': rand_reg()} for s in (1, 2, 3)}
a0 = [rand_reg() for _ in range(NR)]
ref = {}
for s in (1, 2, 3):
    regs = {n: x.copy() for n, x in din[s].items()}; run(words[s], regs, [rand_reg() for _ in range(NR)]); ref[s] = regs
bank = [x.copy() for x in a0]; outs = {}
if CTRL == 'lockstep':
    regsL = {s: {n: x.copy() for n, x in din[s].items()} for s in (1, 2, 3)}
    L = max(len(w) for w in words.values())
    for i in range(L):
        for s in (1, 2, 3):
            if i < len(words[s]): run([words[s][i]], regsL[s], bank)
    outs = regsL
else:
    for s in (1, 2, 3):                    # complete stage-one core, then stage two, then stage three
        regs = {n: x.copy() for n, x in din[s].items()}; run(words[s], regs, bank); outs[s] = regs
bad = sum(not eqreg(outs[s][n], ref[s][n]) for s in (1, 2, 3) for n in 'XY')
check('cross-stage sequential sharing: every stage data == private-scratch run (%d/6 registers wrong)' % bad, bad == 0)
full = [('sx', c, False) for c in range(m)]
for r in range(NR):
    Rops = [op for s in (1, 2, 3) for op in residual(words[s], r)]
    exp = a0[r].copy(); op_apply(Rops, exp)
    okb = eqreg(bank[r], exp)
    tail = op_inv(Rops if CTRL != 'shorttail' else [op for s in (1, 2) for op in residual(words[s], r)]) + full
    z = bank[r].copy(); op_apply(tail, z); fa = a0[r].copy(); op_apply(full, fa)
    okr = eqreg(z, fa)
    rk = set(); unit = True
    for _ in range(4):
        vv = zero_reg(); vv.re[rnd.randrange(N)] = 1; op_apply(tail, vv)
        nz = np.nonzero((vv.re != 0) | (vv.im != 0))[0]; t = len(nz).bit_length() - 1; rk.add(t if len(nz) == 1 << t else -1)
        for i in nz:
            a_, b_ = int(vv.re[i]), int(vv.im[i]); ee = vv.e
            if (a_ * a_ + b_ * b_) << t != 1 << (2 * ee): unit = False
    want = 3 * (JOINS[r] + 1)
    check('role %d (gauge dim %d%s): bank == U3 U2^-1 U1 a0 %s; tail restores F a0 %s; tail rank %s (want 3 dim sigma = %d), unit entries %s' % (
        r, JOINS[r] + 1, ', degenerate' if JOINS[r] == 1 else '', okb, okr, sorted(rk), want, unit), okb and okr and rk == {want} and unit)
print('RESULT h=%d m=%d seed=%d control=%s: %s' % (h, m, SEED, CTRL or '-', 'ALL PASS' if all(OUT.values()) else
      'SOME FAIL: ' + '; '.join(k_ for k_, o in OUT.items() if not o)), flush=True)
