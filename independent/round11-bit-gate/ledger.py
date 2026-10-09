"""[round 11: the round-10 gate ledger.py (all seven of its controls, --expect, .gz words) extended by absorb/aledger.py's terminal-sink absorption (D["absorb"]):
sinks have no register, prelude read, gates or root read; their writes land on the pivot target at the placed
times (move Y_c to the control's current frame, Y_c ^= control); group sinks: pre-shear Y_t ^= Y_c at the pre key
(same frame), post-shear at the post key with every member of T at the root cap U. Extra controls: drop_post,
drop_pre, late_write (a write after its control's next write), stray_write (a write on a non-pivot member).
A pivot swap is NOT a control: any member of T can be the pivot.]
Independent physical ledger of a frozen paired-cube bit word (own code; imports NO construction code: it reads only
the JSON written by export.py and w2-cover-gate/cert.py for the certificate). Modelled on bit-recompile/ledger_c.py.

Registers X_t (data), Y_t (targets), aux s; formal basis 2v + R, every bit a separate symbol (arbitrary x, y AND
arbitrary dirty scratch). Every frame is re-canonicalised exactly over Q (Fraction RREF); a move requires exact
nesting (strict, positive rank), a gate requires both registers at the same exact frame; data of every X/aux register
at every gate and read must lie in the frame (each source label chi_S, exact integer annihilators).

Stage-1 word in time order:
  prelude     every role that is not a kept gauge reads its garbage at frame 0 with the F2 adjoint rows (recomputed
              here from the ops and roots by symbolic propagation);
  V           leaf role -> <chi_x>, z ^= X_x;  partner mixing step 1-2: X_c, X_d -> M, X_c ^= X_d;
  phase 1     its ops at their frames;  centre copies: role -> star span, macro read into Y at 0;
  phase 2     its ops, each kept gauge read at sigma at its exported time (before its first touch, Y_t moves to sigma
              for every response target, reads on the F2 row);
  roots       side roots in order (role and targets -> common cap, Y += role); after the receivers' edge12 root the
              carrier X_c -> C and Y_a, Y_b += X_c;
  cleanup     everything to F, ops reversed, X_c ^= X_d, V removal; Y ends at cap(T) = ker(3 chi_T - 1).
Checks: forward F2 shear (y += x, x and scratch restored), reflected word (reverse time, complemented frames, swapped
banks: continuity and the opposite shear), G = I - J/9 nondegenerate over Q on every frame used, event histogram ==
the exported accounting, #144 child histogram rebuilt from events, D = 2v - 3h(h-2), certificate (1e-12 grid, next
point rejected), stopped saving at theta = 1e-3.  Usage: python3 ledger.py WORD.json [--ctl NAME] [--out OUT.json]"""
import sys, os, json, time, argparse, math
from fractions import Fraction as Q
from collections import Counter, defaultdict
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument('word'); ap.add_argument('--ctl', default=''); ap.add_argument('--out', default='')
ap.add_argument('--theta', default='1/1000'); ap.add_argument('--expect', default='')
A = ap.parse_args(); T0 = time.time()
def log(*a): print('[%5.0fs]' % (time.time() - T0), *a, flush=True)
OK = {}
def report(name, ok):
    OK[name] = bool(ok); log('%-96s %s' % (name, 'PASS' if ok else 'FAIL'))
class Fail(Exception): pass
def need(c, msg):
    if not c: raise Fail(msg)

import gzip; D = json.load(gzip.open(A.word) if A.word.endswith('.gz') else open(A.word)); h, v, R = D['h'], D['v'], D['R']; labels = [tuple(l) for l in D['labels']]
m = 3 * h; chi = [[int(c in T) for c in range(h)] for T in labels]

# ------------------------------------------------------------------ exact frames over Q
def rref(rows):
    M = [[Q(x) for x in r] for r in rows]; piv = []; r = 0
    for c in range(h):
        k = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]; inv = 1 / M[r][c]; M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    return tuple(tuple(x) for x in M[:r]), piv
def nullspace(B):
    """integer basis of {u : b.u = 0 for b in B}"""
    Rr, piv = rref(B) if B else ((), [])
    out = []
    for f in [c for c in range(h) if c not in piv]:
        u = [Q(0)] * h; u[f] = Q(1)
        for row, c in zip(Rr, piv): u[c] = -row[f]
        den = 1
        for x in u: den = den * x.denominator // math.gcd(den, x.denominator)
        ints = [int(x * den) for x in u]; g = 0
        for x in ints: g = math.gcd(g, abs(x))
        out.append([x // g for x in ints])
    return out
key_of = {}; fb = []; fann = []; fdim = []
def canon(rows):
    Rr, _ = rref(rows) if rows else ((), [])
    if Rr not in key_of:
        key_of[Rr] = len(fb); fb.append([[int(x * math.lcm(*[y.denominator for y in row])) for x in row] for row in Rr] if Rr else [])
        fdim.append(len(Rr)); fann.append(nullspace(fb[-1]) if Rr else [[int(i == j) for j in range(h)] for i in range(h)])
    return key_of[Rr]
kZ = canon([]); kF = canon([[int(i == j) for j in range(h)] for i in range(h)])
EX = {int(k): canon(r) for k, r in D['frames'].items()}
log('frames: %d exported -> %d canonical' % (len(D['frames']), len(fb)))
kcap = [canon(nullspace([[3 * x - 1 for x in chi[t]]])) for t in range(v)]
kline = [canon([chi[t]]) for t in range(v)]
dim = lambda k: fdim[k]
_in = {}
def inside(a, b):
    if a == b or b == kF or a == kZ: return True
    if (a, b) not in _in:
        _in[a, b] = dim(a) <= dim(b) and all(sum(x * y for x, y in zip(z, u)) == 0 for z in fann[b] for u in fb[a])
    return _in[a, b]
def bareiss_nonzero(M):
    M = [row[:] for row in M]; n = len(M); prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None: return False
            M[k], M[sw] = M[sw], M[k]
        for i in range(k + 1, n):
            for j in range(k + 1, n): M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return M[n - 1][n - 1] != 0
_nd = {}
def nondeg(k):
    if k not in _nd:
        B = fb[k]
        _nd[k] = True if not B else bareiss_nonzero([[9 * sum(a * b for a, b in zip(u, w)) - sum(u) * sum(w) for w in B] for u in B])
    return _nd[k]
CHI = np.array(chi, dtype=np.int64).T
_allow = {}
def allow(k):
    if k not in _allow:
        if dim(k) == h: _allow[k] = (1 << v) - 1
        elif dim(k) == 0: _allow[k] = 0
        else:
            big = max(abs(x) for z in fann[k] for x in z) >= 1 << 58
            Z = np.array(fann[k], dtype=object if big else np.int64)
            ok = np.all((Z.dot(CHI.astype(object) if big else CHI)) == 0, axis=0)
            _allow[k] = sum(1 << t for t in range(v) if ok[t])
    return _allow[k]

# ------------------------------------------------------------------ word data
ops = [(a, b, EX[f]) for a, b, f in D['ops']]; n1 = D['n_phase1']
sig = {int(s): EX[f] for s, f in D['sigma'].items()}; kept = set(sig); tdeg = {int(s): t for s, t in D['tdeg'].items()}
reads = [(tau, d, s) for tau, d, s in D['reads']]; gt = {int(s): t for s, t in D['gtargets'].items()}
leaf = [(s, x, EX[f]) for s, x, f in D['leafrole']]
roots = D['roots']; mixes = D['mix']
pairs = [tuple(p) for p in D.get('pairs', [])]; donor_of = {b: a for a, b in pairs}
phys = list(range(R))
for s_ in range(R):                       # chains A -> B -> C resolve to the first donor's register
    r_ = s_
    while r_ in donor_of: r_ = donor_of[r_]
    phys[s_] = r_
recips = set(donor_of); owner = list(range(R)); remaining = Counter()
for a_, b_, f_ in ops: remaining[a_] += 1; remaining[b_] += 1
need(not (set(donor_of.values()) & {r['role'] for r in roots}), 'a donor is a root role')
for s, x, f in leaf: need(f == kline[x], 'leaf line frame')

# F2 adjoint rows by symbolic propagation of the garbage symbols through the ops and root reads
cont = [1 << s for s in range(R)]
for a, b, f in ops: cont[a] ^= cont[b]
yacc = [0] * v
for r in roots:
    for t in r['targets']: yacc[t] ^= cont[r['role']]
row = [[] for _ in range(R)]
for t in range(v):
    g = yacc[t]
    while g:
        lb = g & -g; row[lb.bit_length() - 1].append(t); g ^= lb
del cont, yacc
report('gauge read rows lie inside the exported response targets', all(set(row[s]) <= set(gt[s]) for s in kept))

# ------------------------------------------------------------------ forward word
size = 2 * v + R; AUX = 2 * v; MASK = (1 << v) - 1
bits = [1 << i for i in range(size)]
cur = [kline[t] for t in range(v)] + [kZ] * v + [sig[s] if s in sig and s not in recips else kZ for s in range(R)]; init = list(cur)
ev = []; rh = {'aux': Counter(), 'src': Counter(), 'tgt': Counter(), 'center': Counter()}; bad = Counter()
def RG(s): return AUX + phys[s]
def kindof(r): return 'aux' if r >= AUX else 'src' if r < v else 'tgt'
def move(r, dst):
    if cur[r] == dst: return
    need(dim(dst) > dim(cur[r]) and inside(cur[r], dst), ('move not strictly nested', r, dim(cur[r]), dim(dst)))
    rk = dim(dst) - dim(cur[r]); ev.append((0, r, cur[r], dst, rk)); rh[kindof(r)][rk] += 1; cur[r] = dst
def xor(t, s, chk=True):
    need(cur[t] == cur[s], ('gate frames differ', t, s))
    if chk:
        mk = allow(cur[t])
        for r in (t, s):
            if not (v <= r < 2 * v) and (bits[r] & MASK) & ~mk: bad[kindof(r)] += 1
    ev.append((1, t, s, cur[t], 0)); bits[t] ^= bits[s]
ctl = A.ctl; rng = __import__('random').Random(11)
nonkept = [s for s in range(R) if s not in kept]
drop = nonkept[len(nonkept) // 3] if ctl == 'drop_prelude_read' else None
extra = None
if ctl == 'wrong_coefficient':
    s0 = next(s for s in nonkept if row[s]); extra = (s0, next(t for t in range(v) if t not in row[s0]))
ABS = D.get('absorb', []); SINK = {x['s']: x for x in ABS}
# control 'stray_write': the first absorbed write lands on a NON-pivot member of T (outside the shear broadcast)
STRAY = (ABS[0]['place'][0][1], next(t for t in ABS[0]['T'] if t != ABS[0]['c'])) if ctl == 'stray_write' and ABS else None
EXTRA = []                                   # (key, kind, payload)
for x in ABS:
    for n, (key, i) in enumerate(x['place']):
        key = tuple(key)
        if ctl == 'late_write' and x is ABS[0] and n == 0:
            b_ = ops[i][1]; nxt = [k for k in range(i + 1, len(ops)) if ops[k][0] == b_]
            if nxt: key = (nxt[0], 1, 1)
        EXTRA.append((key, 'w', (x, i)))
    if len(x['T']) > 1:
        if not (ctl == 'drop_pre' and x is ABS[0]): EXTRA.append((tuple(x['pre']), 'pre', x))
        EXTRA.append((tuple(x['post']), 'post', x))
def run_extra(kind, x):
    if kind == 'w':
        x, i = x; b_ = ops[i][1]
        need(owner[phys[b_]] == b_, ('absorbed write after its control lost its register', b_, i)); remaining[b_] -= 1
        tg_ = STRAY[1] if STRAY and x is ABS[0] and i == STRAY[0] else x['c']
        f = cur[RG(b_)]; move(v + tg_, f); xor(v + tg_, RG(b_))
    elif kind == 'pre':                     # all of T to the join of their current frames, then Y_t ^= Y_c
        phi = canon([r_ for t in x['T'] for r_ in fb[cur[v + t]]])
        for t in x['T']: move(v + t, phi)
        for t in x['T']:
            if t != x['c']: xor(v + t, v + x['c'])
    else:
        U = EX[roots[x['j']]['frame']]
        for t in x['T']: move(v + t, U)
        if not (ctl == 'drop_post' and x is ABS[0]):
            for t in x['T']:
                if t != x['c']: xor(v + t, v + x['c'])
EXTRA.sort(key=lambda e: e[0]); xi = 0
def flush(key):
    global xi
    while xi < len(EXTRA) and EXTRA[xi][0] < key: run_extra(EXTRA[xi][1], EXTRA[xi][2]); xi += 1
try:
    for s in range(R):                                  # prelude at frame 0
        if s in kept or s == drop or s in SINK: continue
        tg = row[s] + ([extra[1]] if extra and extra[0] == s else [])
        for t in tg: move(v + t, kZ); xor(v + t, AUX + s)
    for s, x, f in leaf: move(AUX + s, f); xor(AUX + s, x)
    for mx in mixes:
        c, d = mx['carrier'], mx['passive']; move(c, EX[mx['M']]); move(d, EX[mx['M']]); xor(c, d)
    if ctl == 'desc_shrink':                             # one op frame loses a basis row
        i = next(i for i, (a, b, f) in enumerate(ops) if dim(f) > 2 and i >= n1 and a not in SINK)  # a sink's gates are skipped
        a, b, f = ops[i]; ops[i] = (a, b, canon(fb[f][:-1])); log('ctl: op %d frame dim %d -> %d' % (i, dim(f), dim(ops[i][2])))
    if ctl == 'early_birth':              # control: a recipient's readout moved before its donor's last gate
        b0 = sorted(recips)[len(recips) // 2]; a0 = donor_of[b0]
        lastA = max(i for i, (a, b, f) in enumerate(ops) if a0 in (a, b))
        reads = sorted([(lastA, d, s) if s == b0 else (t_, d, s) for t_, d, s in reads]); log('ctl: birth of %d moved before op %d' % (b0, lastA))
    if ctl == 'late_readout':
        tau0, d0, s0 = reads[0]; reads[0] = (tau0 + 1 + next(k for k, (a, b, f) in enumerate(ops[tau0:]) if s0 in (a, b)), d0, s0)
        reads.sort()
    def gate(i):
        a, b, f = ops[i]
        if a in SINK: return
        for z in (a, b):
            need(owner[phys[z]] == z, ('register used by a role that does not own it', z, i)); remaining[z] -= 1
        move(RG(a), f); move(RG(b), f); xor(RG(a), RG(b))
    for i in range(n1): gate(i)
    for j, r in enumerate(roots):                        # centre copies
        if r['kind'] != 'center': continue
        s = r['role']; f = EX[r['frame']]; need(owner[phys[s]] == s, 'centre role'); move(RG(s), f)
        need(all(cur[v + t] == kZ for t in r['targets']), 'centre targets not at 0')
        if (bits[RG(s)] & MASK) & ~allow(f): bad['center'] += 1
        ev.append((2, RG(s), f, kZ, dim(f), tuple(r['targets']))); rh['center'][dim(f)] += 1
        for t in r['targets']: bits[v + t] ^= bits[RG(s)]
    flush((n1, 0))
    ri = 0
    def readout(s):
        r = RG(s); f = sig[s]
        if s in recips:
            a = donor_of[s]; need(owner[phys[s]] == a and remaining[a] == 0, ('birth before the donor died', s, a))
            move(r, f); owner[phys[s]] = s
            if (bits[r] & MASK) & ~allow(f): bad['birth'] += 1
        else:
            need(not (bits[r] & MASK), ('gauge role holds data at its readout', s))
        need(cur[r] == f, ('readout frame', s))
        for t in gt[s]: move(v + t, f)
        for t in row[s]: xor(v + t, r)
    L = len(ops)
    for i in range(n1, L):
        while ri < len(reads) and max(reads[ri][0], n1) <= i:
            flush((max(reads[ri][0], n1), 0, ri)); readout(reads[ri][2]); ri += 1
        flush((i, 0.75)); gate(i)
    while ri < len(reads): flush((L, 0, ri)); readout(reads[ri][2]); ri += 1
    bymix = defaultdict(list)
    for mx in mixes: bymix[mx['root']].append(mx)
    for j, r in enumerate(roots):
        if r['kind'] != 'side': continue
        flush((L + 1, j, 0))
        s = r['role']; f = EX[r['frame']]
        if s in SINK:                                    # the root read is gone; group members still reach U here
            if len(r['targets']) > 1:
                for t in r['targets']: move(v + t, f)
        else:
            need(owner[phys[s]] == s, 'root role'); move(RG(s), f)
            for t in r['targets']: move(v + t, f); xor(v + t, RG(s))
        flush((L + 1, j, 1))
        for mx in bymix.get(j, ()):
            c = mx['carrier']; C = EX[mx['C']]
            if ctl == 'below_level': C = EX[mx['M']]
            else: need(C == f, 'mixing delivery frame is not the receivers\' root frame')
            move(c, C)
            for t in mx['receivers']: xor(v + t, c)
    flush((L + 9,)); need(xi == len(EXTRA), 'unplaced absorption events')
    for t in range(v): move(v + t, kcap[t])
    for s in range(R):
        if s not in recips and s not in SINK: move(AUX + s, kF)
    for a, b, f in reversed(ops):
        if a not in SINK: xor(RG(a), RG(b), False)
    for mx in mixes:
        c, d = mx['carrier'], mx['passive']; move(c, kF); move(d, kF)
        if ctl != 'broken_mix': xor(c, d, False)
    for s, x, f in leaf: xor(AUX + s, x, False)
    for t in range(v): move(t, kF)
    fwd = all(bits[i] == ((1 << i) ^ (1 << (i - v)) if v <= i < 2 * v else 1 << i) for i in range(size))
    report('complete forward F2 shear on all %d formal symbols (y_T += x_T; x and dirty scratch restored)' % size, fwd)
    report('data inside the exact frame at every gate, read and centre copy %s' % dict(bad), not bad)
    report('X and aux end at F, every Y at cap(T) = ker(3 chi_T - 1)', all(cur[i] == kF for i in range(v)) and
           all(cur[AUX + s] == kF for s in range(R) if s not in recips and s not in SINK) and all(cur[v + t] == kcap[t] for t in range(v)))
except Fail as e:
    report('forward word executes (%s)' % (e,), False)
    log('ALL FAIL', A.ctl); json.dump(dict(OK=OK, all=False, ctl=A.ctl, error=str(e)), open(A.out, 'w')) if A.out else None; sys.exit(0)
usedf = set(cur) | {e[2] for e in ev if e[0] == 0} | {e[3] for e in ev if e[0] in (0, 1)} | {e[2] for e in ev if e[0] == 2}
nd_bad = sum(1 for k in usedf if not nondeg(k))
report('G = I - J/9 nondegenerate over Q on every frame used (%d frames)' % len(usedf), nd_bad == 0)

# ------------------------------------------------------------------ reflected word
bank = lambda r: r + v if r < v else r - v if r < 2 * v else r
rev = [None] * size
for r in range(size): rev[bank(r)] = ~cur[r]
bits = [1 << i for i in range(size)]; rrh = Counter(); okr = True
for e in reversed(ev):
    k = e[0]
    if k == 0:
        _, r, b_, c_, rk = e; rr = bank(r); okr &= rev[rr] == ~c_; rev[rr] = ~b_; rrh[rk] += 1
    elif k == 1:
        _, r, s, c_, _ = e; rr, sr = bank(r), bank(s); okr &= rev[rr] == rev[sr] == ~c_; bits[rr] ^= bits[sr]
    else:
        _, r, f, z, rk, tg = e; rr = bank(r); okr &= rev[rr] == ~f
        for t in tg: okr &= rev[bank(v + t)] == ~z; bits[bank(v + t)] ^= bits[rr]
        rrh[rk] += 1
okr &= all(rev[bank(r)] == ~init[r] for r in range(size))
okb = all(bits[i] == ((1 << i) ^ (1 << (i + v)) if i < v else 1 << i) for i in range(size))
report('reflected word: frame continuity (complement + bank swap), returns to the complemented start', okr)
report('reflected word: complete opposite F2 shear (x_T += y_T, scratch restored)', okb)
report('reflected rank multiset == forward rank multiset', rrh == sum(rh.values(), Counter()))

# ------------------------------------------------------------------ histograms, certificate
H = Counter(); Wv = Q(2 * v)
for part in ('aux', 'src', 'tgt', 'center'):
    for r_, n_ in rh[part].items(): H[r_] += 3 * n_
for s in range(R):
    if s in recips or s in SINK: continue
    if s in kept:
        t = tdeg[s]; u = h - dim(sig[s]); Wv += Q(3, t)
        if m - t * u: H[m - t * u] += Q(3, t)
    else: Wv += 1
H[2] += 2 * v; H = +H
acc = {int(r): Q(c) for r, c in D['accounting']['hist'].items()}
report('child histogram from events == exported accounting (W %s)' % Wv, dict(H) == acc and Wv == Q(D['accounting']['W']))
mass = sum(r * n for r, n in H.items()); Dv = Wv * m - mass
report('deficit per vertex == 2v - 3h(h-2) = %d' % (2 * v - 3 * h * (h - 2)), Dv == 2 * v - 3 * h * (h - 2))
sys.path.insert(0, HERE); import cert   # local copy of w2-cover-gate/cert.py
K = math.lcm(Wv.denominator, *[n.denominator for n in H.values()])
c = cert.certify({r: int(n * K) for r, n in H.items()}, int(Wv * K), m, 10 ** 12)
th = Q(A.theta); ab = (1 - th) * c['a'] + th * cert.AOLD
report('certificate: a* = %s (%.10e) certified, next grid point rejected; stopped a_bit = %.10e' % (c['a'], float(c['a']), float(ab)),
       c['certified'] and c['next_rejected'])
if A.expect: report('certified a* equals the expected %s' % A.expect, c['a'] == Q(A.expect))
res = dict(OK=OK, all=all(OK.values()), ctl=A.ctl, a=str(c['a']), a_bit=str(ab), a_bit_f=float(ab), W=str(Wv),
           hist={str(r): str(n) for r, n in sorted(H.items())}, events=len(ev), frames=len(fb))
log('ALL', 'PASS' if res['all'] else 'FAIL ' + str([k for k, x in OK.items() if not x]), A.ctl)
if A.out: json.dump(res, open(A.out, 'w'), indent=1)
