"""Independent checker of a FROZEN paired-cube complex word with birth-read reuse (cube-prefix release gate).

Reads only the word JSON written by export_word.py; imports no construction code. Checks:
  1. decoder identity mod P = 2^61-1 on random x: (sum of root coefficients x DAG values) + K x = x, K = the cube-local
     Hadamard/2 blocks of #144 (antipodal port +1/2, neighbour port -1/2);
  2. legality on actual F2 subspaces: every op frame contains its node's span; every role chain
     start -> op frames -> root frame -> [recipient gauge -> its op frames -> ...] -> full is nested; every pair has a
     non-root donor with a last op before the recipient's read clock, the read no later than the recipient's first op,
     and the donor's last frame inside the recipient gauge; gauge reads keep every target chain nested in read order;
  3. the ledger: child histogram and W recomputed from the chains (#144 normalisation: 3 x positive chain steps,
     3 x source [1, 2, h-4], 3 x target steps, 2v width-two finals, plain 3d tail for any unpaired gauge), the
     telescoping deficit W m - s = 2v - 3 loss, equality with the claim;
  4. the exact saving on the 1e-12 grid (rational ln upper bound, exp upper bound), next grid point rejected, equal
     to the claim;
  5. exact aliased dirty-scratch replay mod P, two seeds (every slot restored, y - y0 = the DAG's root sum), and three
     controls that must fail: a recipient's read omitted, a read after the recipient's first op, a read before the
     donor's last write.
Usage: python3 check_word.py WORD.json.gz"""
import sys, gzip, json, math, random, time
from fractions import Fraction as Q
from collections import Counter, defaultdict

P = (1 << 61) - 1
t0 = time.time()
def log(*a): print('[%4.0fs]' % (time.time() - t0), *a, flush=True)
J = json.load(gzip.open(sys.argv[1], 'rt'))
h, v, m, R = J['h'], J['v'], J['m'], J['R']; args = J['args']; signs = J['signs']; inputs = J['inputs']
roots = J['roots']; ops = [tuple(o) for o in J['ops']]; sched = J['sched']; cut = J['cut']
sources = {int(x): s for x, s in J['sources'].items()}; rootroles = J['rootroles']
gauges = J['gauges']; tau = {int(b): t for b, t in J['tau'].items()}; pairs = [tuple(q) for q in J['pairs']]
frames = [tuple(f) for f in J['frames']]
OK = {}
def check(name, ok): OK[name] = bool(ok); log(('PASS ' if ok else 'FAIL ') + name)

# ---------------- F2 subspaces ----------------
def basis(rows):
    b = {}
    for x in rows:
        for q in sorted(b, reverse=True):
            if x >> q & 1: x ^= b[q]
        if x:
            q = x.bit_length() - 1
            for k in list(b):
                if b[k] >> q & 1: b[k] ^= x
            b[q] = x
    return tuple(b[q] for q in sorted(b, reverse=True))
def inside(A, B):
    bb = basis(B); piv = {y.bit_length() - 1: y for y in bb}
    for x in A:
        for q in sorted(piv, reverse=True):
            if x >> q & 1: x ^= piv[q]
        if x: return False
    return True
def perp(rows):
    rows = basis(rows); piv = {r.bit_length() - 1: r for r in rows}; out = []
    for j in range(h):
        if j in piv: continue
        x = 1 << j
        for q, r in piv.items():
            if r >> j & 1: x |= 1 << q
        out.append(x)
    return basis(out)
FULL = tuple(1 << j for j in range(h - 1, -1, -1))
md = lambda q: Q(q).numerator % P * pow(Q(q).denominator % P, P - 2, P) % P
n = len(args); assert args[0] is None and all(args[x] is None for x in range(1, v + 1))
spans = [()] * n
for x in range(1, n): spans[x] = (inputs[x - 1],) if args[x] is None else basis(spans[args[x][0]] + spans[args[x][1]])
def dag(x):
    val = [0] * n
    for k in range(1, n):
        val[k] = x[k - 1] if args[k] is None else (val[args[k][0]] + signs[k] * val[args[k][1]]) % P
    return val

# ---------------- 1. decoder identity ----------------
rng = random.Random(7); x = [rng.randrange(P) for _ in range(v)]; val = dag(x); y = [0] * v
for r in roots:
    for t, c in zip(r['targets'], r['coefficients']): y[t] = (y[t] + md(c) * val[r['node']]) % P
half = md(Q(1, 2)); bad = 0
for t in range(v):
    k = t - t % 8; acc = 0
    for s in range(k, k + 8):
        w = bin(inputs[t] ^ inputs[s]).count('1')
        acc += x[s] if w == 6 else -x[s] if w == 2 else 0
    bad += (y[t] + acc * half - x[t]) % P != 0
check('decoder identity H x + K x = x on random x mod P (%d wrong of %d)' % (bad, v), bad == 0)

# ---------------- 2. roles, chains, legality ----------------
pos = {i: k for k, i in enumerate(sched)}; assert sorted(sched) == list(range(len(ops)))
rops = [[] for _ in range(R)]
for i in sched:
    a, b, _ = ops[i]; assert a != b; rops[a].append(i); rops[b].append(i)
rootframe, rootkind = {}, {}
for r, sr in zip(roots, rootroles):
    assert sr not in rootframe
    if r['kind'] == 'center': rootframe[sr] = spans[r['node']]; rootkind[sr] = 'center'
    else: rootframe[sr] = perp([inputs[t] for t in r['targets']]); rootkind[sr] = 'side'
start = [()] * R
for xx, sr in sources.items(): start[sr] = (inputs[xx - 1],)
gsig = {}
for z in gauges: gsig[z['role']] = perp(z['A']); start[z['role']] = gsig[z['role']]
donor = {a: b for a, b in pairs}; recip = {b: a for a, b in pairs}
check('pairs are a matching (%d)' % len(pairs), len(donor) == len(recip) == len(pairs))
bad = Counter()
for i, (_, _, xx) in enumerate(ops):
    if not inside(spans[xx], frames[i]): bad['span'] += 1
for a, b in pairs:
    if a in rootframe or not rops[a]: bad['donor kind'] += 1
    if b not in gsig: bad['recipient kind'] += 1
    if not pos[rops[a][-1]] < tau[b]: bad['chronology'] += 1
    if rops[b] and not tau[b] <= pos[rops[b][0]]: bad['read after first touch'] += 1
    if not inside(frames[rops[a][-1]], gsig[b]): bad['donor frame outside gauge'] += 1
for z in gauges:
    b = z['role']
    if rops[b] and not tau[b] <= pos[rops[b][0]]: bad['gauge read after first touch'] += 1
def chain(sr): return [start[sr]] + [frames[i] for i in rops[sr]] + ([rootframe[sr]] if sr in rootframe else [])
local = Counter()
for sr in range(R):
    if sr in recip: continue
    seq = []; cur = sr
    while cur is not None: seq += chain(cur); cur = donor.get(cur)
    seq.append(FULL)
    for A_, B_ in zip(seq, seq[1:]):
        if not inside(A_, B_): bad['nesting'] += 1
    if len(start[sr]) == 1 and sr not in gsig: local[1] += 1
    d = [len(X) for X in seq]
    for d0, d1 in zip(d, d[1:]):
        if d1 > d0: local[d1 - d0] += 1
    cur = sr
    while cur is not None:
        if rootkind.get(cur) == 'center':
            if len(rootframe[cur]) != h - 2: bad['centre frame'] += 1
            local[len(rootframe[cur])] += 1
        cur = donor.get(cur)
order = sorted(range(len(gauges)), key=lambda k: (tau[gauges[k]['role']], len(gauges) - 1 - k))
curA = [FULL] * v; target = Counter()
for k in order:
    A = basis(gauges[k]['A'])
    for t in gauges[k]['targets']:
        if not inside(A, curA[t]): bad['target order'] += 1
        target[len(curA[t]) - len(A)] += 1; curA[t] = A
for r in roots:
    if r['kind'] == 'side':
        A = basis([inputs[t] for t in r['targets']])
        for t in r['targets']:
            if not inside(A, curA[t]): bad['side order'] += 1
            target[len(curA[t]) - len(A)] += 1; curA[t] = A
for t in range(v):
    if not inside((inputs[t],), curA[t]): bad['input'] += 1
    target[len(curA[t]) - 1] += 1
check('legality: spans, chains, pairs, chronology, target order (violations %s)' % dict(bad), not bad)

# ---------------- 3. ledger ----------------
C = Counter()
for hist in (local, {1: v, 2: v, h - 4: v}, target):
    for r_, k in hist.items():
        if r_: C[r_] += Q(3 * k)
C[2] += 2 * v
Wv = Q(2 * v + R - len(pairs))
for z in gauges:
    if z['role'] not in recip: C[3 * (h - len(basis(z['A'])))] += 1
C = {r_: k for r_, k in C.items() if k}
s_ = sum(r_ * k for r_, k in C.items()); loss = J['loss']
check('telescoping deficit W m - s = 2v - 3 loss (%s = %d)' % (Wv * m - s_, 2 * v - 3 * loss), Wv * m - s_ == 2 * v - 3 * loss)
check('every child narrower than m (max %d < %d)' % (max(C), m), max(C) < m and min(C) >= 1)
claimC = {int(r_): Q(k) for r_, k in J['claim']['C'].items()}
check('child histogram and W equal the claim', C == claimC and Wv == Q(J['claim']['W']))

# ---------------- 4. certified saving ----------------
def ln_upper(z, N=28):
    z = Q(z); k = 0
    while z > 2: z /= 2; k += 1
    def ta(u): return 2 * (sum(u ** (2 * j + 1) / (2 * j + 1) for j in range(N)) + u ** (2 * N + 1) / ((2 * N + 1) * (1 - u * u)))
    return k * ta(Q(1, 3)) + (ta((z - 1) / (z + 1)) if z > 1 else 0)
lu = {w: ln_upper(Q(m, w)) for w in C}
def Fup(a):
    tot = Q(0)
    for w, k in C.items():
        u = a * lu[w]; assert 0 <= u < 1
        tot += k * Q(w, m) * (1 + u + u * u / (2 * (1 - u / 3)))
    return tot / Wv
a = Q(J['claim']['a_c']); den = 10 ** 12
check('a_c = %s: grid point accepted (F < 1), denominator divides 1e12' % a, (a * den).denominator == 1 and Fup(a) < 1)
check('next grid point %s rejected (F >= 1)' % (a + Q(1, den)), not Fup(a + Q(1, den)) < 1)

# ---------------- 5. exact dirty-scratch replay ----------------
hold = {sr: xx for xx, sr in sources.items()}; coef = {}
for i in sched:
    d, c, xx = ops[i]
    if hold.get(d) is None: assert hold[c] == xx; coef[i] = (1, 1)
    else:
        aa, bb = args[xx]; sg = signs[xx]
        if hold[d] == aa: assert hold[c] == bb; coef[i] = (1, sg)
        else: assert hold[d] == bb and hold[c] == aa; coef[i] = (sg, 1)
    hold[d] = xx
seeds = [(sr, r['kind'], {t: md(c) for t, c in zip(r['targets'], r['coefficients'])}, r['node']) for r, sr in zip(roots, rootroles)]
cf = [dict() for _ in range(R)]
def addto(dst, src, k):
    for t, u in src.items(): dst[t] = (dst.get(t, 0) + k * u) % P
for sr, kind, sd, _ in seeds:
    if kind == 'side': addto(cf[sr], sd, 1)
for k in range(len(sched) - 1, -1, -1):
    i = sched[k]; d, c, _ = ops[i]; ca, cb = coef[i]
    if k == cut - 1:
        for sr, kind, sd, _ in seeds:
            if kind == 'center': addto(cf[sr], sd, 1)
    addto(cf[c], cf[d], cb)
    if ca != 1: cf[d] = {t: ca * u % P for t, u in cf[d].items()}
if cut == 0:
    for sr, kind, sd, _ in seeds:
        if kind == 'center': addto(cf[sr], sd, 1)
gauged = {z['role'] for z in gauges}
def run(seed, omit=None, after_first=None, before_death=None):
    rnd = random.Random(seed)
    def rs(sr):
        while sr in recip: sr = recip[sr]
        return sr
    slot = {sr: rs(sr) for sr in range(R)}; phys = sorted(set(slot.values()))
    z = {q: rnd.randrange(P) for q in phys}; a_ = dict(z)
    x = [rnd.randrange(P) for _ in range(v)]; y0 = [rnd.randrange(P) for _ in range(v)]; y = list(y0)
    def read(sr, sg, vec):
        val_ = a_[slot[sr]]
        for t, u in vec.items(): y[t] = (y[t] + sg * u * val_) % P
    rt = defaultdict(list)
    for b in gauged:
        t = tau[b]
        if b == after_first: t = pos[rops[b][0]] + 1
        if b == before_death: t = max(pos[i] for i in rops[recip[b]] if ops[i][0] == recip[b])
        rt[t].append(b)
    for sr in range(R):
        if sr not in gauged: read(sr, -1, cf[sr])
    for xx, sr in sources.items(): a_[slot[sr]] = (a_[slot[sr]] + x[xx - 1]) % P
    chron = []
    for k, i in enumerate(sched):
        if k == cut:
            for sr, kind, sd, _ in seeds:
                if kind == 'center': read(sr, +1, sd)
        for b in rt.get(k, ()):
            if b != omit: read(b, -1, cf[b])
        d, c, _ = ops[i]; ca, cb = coef[i]; dd, cc = slot[d], slot[c]
        if dd == cc: return (False, False)
        a_[dd] = (ca * a_[dd] + cb * a_[cc]) % P; chron.append((dd, cc, ca, cb))
    if cut == len(sched):
        for sr, kind, sd, _ in seeds:
            if kind == 'center': read(sr, +1, sd)
    for b in rt.get(len(sched), ()):
        if b != omit: read(b, -1, cf[b])
    for sr, kind, sd, _ in seeds:
        if kind == 'side': read(sr, +1, sd)
    for dd, cc, ca, cb in reversed(chron): a_[dd] = (a_[dd] - cb * a_[cc]) * pow(ca % P, P - 2, P) % P
    for xx, sr in sources.items(): a_[slot[sr]] = (a_[slot[sr]] - x[xx - 1]) % P
    vv = dag(x); want = list(y0)
    for sr, kind, sd, node in seeds:
        for t, u in sd.items(): want[t] = (want[t] + u * vv[node]) % P
    return a_ == z, y == want
for sd in (1, 2): check('dirty-scratch replay seed %d (%d slots): slots restored and y correct' % (sd, R - len(pairs)), run(sd) == (True, True))
if pairs:
    check('control omitted recipient read rejected', run(3, omit=pairs[0][1]) != (True, True))
    b = next(b for _, b in pairs if rops[b])
    check('control read after first touch rejected', run(4, after_first=b) != (True, True))
    b = next(b for a_, b in reversed(pairs) if any(ops[i][0] == a_ for i in rops[a_]))
    check('control read before donor last write rejected', run(5, before_death=b) != (True, True))
log('RESULT check_word %s p=%d a_c=%s: %s' % (J['name'], J['p'], a, 'ALL PASS' if all(OK.values()) else 'FAIL'))
sys.exit(0 if all(OK.values()) else 1)
