"""Physical bit ledger of witness 2 (jfstack, h=23) with #144-style gauge omission. Own code; reads only our frozen
round-7 schedules (certificates/round7) through our jfdata reader.

Registers: X_t (data, leaf t+1), Y_t (targets), aux s (R roles); formal basis 2v+R, every bit a separate symbol
(arbitrary dirty scratch AND arbitrary x, y). Every register has an explicit path of exact canonical frames; a move
walks the path forward (each step charged with its rank), a gate requires both registers at the same exact frame.

Stage-1 word in time order (the modified word):
  prelude   every non-deferred role AND every omitted deferred role reads its garbage out at frame 0, with the exact
            adjoint coefficients (the F2 rows, recomputed symbolically and compared), before any positive-rank op;
  early V, phase 1, centre copies (macro U_c -> 0, rank 22), retained deferred readouts at sigma_u (dim order),
  late V, phase 2, output reads at t_T^perp; then L^-1 and V^-1 at the full frame.
Checks: complete forward F2 shear (y += x, everything else restored), exact nesting of every consecutive pair of every
path, data inside the frame at every aux gate, target chains = deleted subsequences of the certified chains,
reflected word (reverse time, complement frames, swap banks) continuity and the opposite shear, event rank
histograms == the accounting of omit.py, #144 child histogram rebuilt from events.

Usage: python3 ledger.py [--omit PKL] [--support F2|Z] [--ctl NAME] [--out JSON]"""
import sys, os, time, json, pickle, argparse
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.join(HERE, '..', '..')
sys.path.insert(0, REPO + '/independent/joint-frame-stack')
import jfdata
from xq import rref, ann, inside_ann

ap = argparse.ArgumentParser()
ap.add_argument('--omit', default=os.path.join(REPO, 'certificates', 'round9', 'w2_omission.json'))
ap.add_argument('--support', default='F2')
ap.add_argument('--ctl', default='')
ap.add_argument('--out', default='')
ap.add_argument('--nest', type=int, default=1)
ap.add_argument('--data', type=int, default=1)
A = ap.parse_args()
T0 = time.time()
def log(*a): print('[%5.0fs]' % (time.time() - T0), *a, flush=True)
OK = {}
def report(name, ok):
    OK[name] = bool(ok); log('%-100s %s' % (name, 'PASS' if ok else 'FAIL'))

class Fail(Exception): pass
def need(c, msg):
    if not c: raise Fail(msg)

K = jfdata.load(); EX, eqid = jfdata.load_frames()
X = K['X']; W = K['W']; G = X['G']; h = 23
trip = [tuple(t) for t in G['trip']]; v = len(trip); R = W['R']; tid = {T: i for i, T in enumerate(trip)}
ops = W['ops']; srcop = W['srcop']; f = W['f']; ro = jfdata.release_orders(K)
sel = ro['sel']; early_V, late_V, defer_order = ro['early_V'], ro['late_V'], ro['defer_order']
omitted = set(json.load(open(A.omit))['omitted_slots']) if A.omit != 'none' else set()
need(omitted <= sel, 'omitted roles must be deferred roles')
kept_def = sel - omitted
log('h %d v %d R %d deferred %d omitted %d retained %d' % (h, v, R, len(sel), len(omitted), len(kept_def)))

# ---------------------------------------------------------------- exact canonical frames
FULL = tuple(tuple(int(i == j) for j in range(h)) for i in range(h)); ZERO = ()
fid = {}; fb = []
def key(B):
    B = rref(B, h)
    if B not in fid: fid[B] = len(fb); fb.append(B)
    return fid[B]
kZ, kF = key(ZERO), key(FULL)
exk = {}
def kex(i):
    if i not in exk: exk[i] = key(EX[i])
    return exk[i]
ksig = {s: key(K['SIG'][s]) for s in sel}; kvs = {s: key(B) for s, B in K['VS'].items()}
tv = lambda T: [int(q in T) for q in range(h)]
kline = [key([tv(T)]) for T in trip]
ktperp = [key(ann(rref([[9 * x - 3 for x in tv(T)]], h), h)) for T in trip]
dim = lambda k: len(fb[k])
log('canonical frames', len(fb))

# ---------------------------------------------------------------- adjoint coefficients (exact, symbolic F2)
a = [1 << (v + s) for s in range(R)]
for op in ops:
    if op[0] == 'src': a[op[1]] ^= 1 << (op[2] - 1)
    else: a[op[1]] ^= a[op[2]]
y = [0] * v
for s, (c, T) in X['out'].items(): y[tid[T]] ^= a[s]
for s, c in X['ret'].items():
    for T in trip:
        if c in T: y[tid[T]] ^= a[s]
rowsF2 = [0] * R
for t, yy in enumerate(y):
    g = yy >> v
    while g:
        lb = g & -g; rowsF2[lb.bit_length() - 1] |= 1 << t; g ^= lb
del a, y
report('adjoint: own symbolic F2 garbage coefficients == frozen readout rows (%d roles)' % R, rowsF2 == list(W['F2']))
rowsread = list(rowsF2)                       # what the word actually reads (the F2 circuit)
if A.support == 'Z':
    sys.exit('the Z-support variant is not part of the release (the F2 rows are the circuit)')
    rowschain = rowsZ                          # #144's conservative convention: chains keep every positive Z support
else:
    rowschain = rowsF2
bits_of = lambda g: [t for t in range(v) if g >> t & 1]

# ---------------------------------------------------------------- paths
walk = jfdata.walks(K, eqid, ro)
def wkey(k):
    if isinstance(k, tuple) and k[0] == 'sig': return ksig[k[1]]
    if isinstance(k, tuple) and k[0] == 'v': return kvs[k[1]]
    return kex(k)
def dedup(p): return [k for i, k in enumerate(p) if i == 0 or k != p[i - 1]]
paths = [None] * (2 * v + R)
# X_t: <t_S>, V starts in time order (frozen X_S order = early then deferred by dim, checked by check_schedule), F
for x, order in W['xorder'].items():
    paths[x - 1] = dedup([kline[x - 1]] + [kvs[s] for s in order] + [kF])
# Y_t: 0, the retained deferred readouts touching t in time order (= deleted subsequence of the certified chain), t_T^perp
ych = [[kZ] for _ in range(v)]
for s in defer_order:
    if s in kept_def:
        for t in bits_of(rowschain[s]): ych[t].append(ksig[s])
for t in range(v): paths[v + t] = dedup(ych[t] + [ktperp[t]])
for s in range(R):
    wk = [wkey(k) for k in walk[s]]
    if s in sel and s in omitted: wk = wk[1:]          # omitted: sigma entry removed, chain starts at 0
    p = wk if s in kept_def else [kZ] + wk
    paths[2 * v + s] = dedup(p + [kF])
# control: a target chain that is NOT a subsequence of the certified chain (two deferred levels swapped)
if A.ctl == 'chain_not_subsequence':
    t = next(t for t in range(v) if len(paths[v + t]) >= 4 and dim(paths[v + t][1]) < dim(paths[v + t][2]))
    p = paths[v + t]; p[1], p[2] = p[2], p[1]; log('ctl: swapped levels %d,%d on target %d' % (dim(p[1]), dim(p[2]), t))

# exact nesting of every consecutive path pair
if A.nest:
    pairs = set()
    for p in paths:
        for a_, b_ in zip(p, p[1:]): pairs.add((a_, b_))
    _Z = {}
    def inside(ka, kb):
        if dim(kb) == h or dim(ka) == 0: return True
        if kb not in _Z: _Z[kb] = ann(fb[kb], h)
        return inside_ann(fb[ka], _Z[kb])
    badn = sum(1 for a_, b_ in pairs if not (dim(a_) < dim(b_) and inside(a_, b_)))
    report('exact strict nesting of every consecutive frame pair of every path (%d distinct pairs, %d paths)' % (len(pairs), len(paths)), badn == 0)

# deleted-subsequence check of the target chains against the certified (unmodified) chains
full_y = [[kZ] for _ in range(v)]
for s in defer_order:
    for t in bits_of(rowsF2[s] if A.support == 'F2' else rowschain[s]): full_y[t].append(ksig[s])
def is_subseq(p, q):
    it = iter(q); return all(any(x == z for z in it) for x in p)
report('target chains are deleted subsequences of the full certified chains (%d targets)' % v,
       all(is_subseq(paths[v + t], dedup(full_y[t] + [ktperp[t]])) for t in range(v)))

# ---------------------------------------------------------------- data-in-frame masks (aux gates)
MASK = (1 << v) - 1
_allow = {}
def allow(k):
    if k not in _allow:
        if dim(k) == h: _allow[k] = MASK
        else:
            Zk = ann(fb[k], h) if dim(k) else FULL
            _allow[k] = sum(1 << t for t, T in enumerate(trip) if all(z[T[0]] + z[T[1]] + z[T[2]] == 0 for z in Zk))
    return _allow[k]

# ---------------------------------------------------------------- forward word
size = 2 * v + R; AUX = 2 * v
pos = [0] * size; cur = [p[0] for p in paths]; init = list(cur)
bits = [1 << i for i in range(size)]
ev = []                                   # (kind, reg, from/frame, to, rank)  kind 0 move, 1 xor, 2 centre macro
rh = {'aux': Counter(), 'src': Counter(), 'tgt': Counter(), 'center': Counter()}
baddata = Counter()
def move(r, dst):
    if cur[r] == dst: return
    p = paths[r]
    need(dst in p[pos[r] + 1:], ('move off path', r, dim(cur[r]), dim(dst)))
    j = p.index(dst, pos[r] + 1)
    while pos[r] < j:
        a_, b_ = p[pos[r]], p[pos[r] + 1]; rk = dim(b_) - dim(a_); need(rk > 0, 'non-increasing step')
        ev.append((0, r, a_, b_, rk))
        rh['aux' if r >= AUX else 'src' if r < v else 'tgt'][rk] += 1
        pos[r] += 1; cur[r] = b_
def xor(t, s, chk_data=False):
    need(cur[t] == cur[s], ('gate frames differ', t, s))
    if chk_data and A.data:                     # data of every non-target register at the gate inside the frame
        m_ = allow(cur[t])
        for r in (t, s):
            if not (v <= r < 2 * v) and (bits[r] & MASK) & ~m_: baddata[(chk_data, 'X' if r < v else 'aux')] += 1
    ev.append((1, t, s, cur[t], 0)); bits[t] ^= bits[s]
rows_for = (lambda s: bits_of(rowsread[s]))
if A.ctl == 'wrong_coefficient':
    s_bad = sorted(omitted)[len(omitted) // 2]; t_add = next(t for t in range(v) if not rowsread[s_bad] >> t & 1)
    rows_for = (lambda s, _r=rows_for: _r(s) + ([t_add] if s == s_bad else []))
    log('ctl: role %d reads an extra target %d (coefficient 1 instead of 0)' % (s_bad, t_add))
drop = sorted(omitted)[len(omitted) // 3] if A.ctl == 'drop_prelude_read' else sorted(omitted)[0] if A.ctl == 'late_prelude' else None
def readout(s, frame):
    r = AUX + s
    need(cur[r] == frame, ('readout frame', s))
    need(not (bits[r] & MASK), ('deferred/prelude role holds data at its readout', s))
    for t in rows_for(s):
        move(v + t, frame); xor(v + t, r)
# 1. prelude at frame 0
for s in range(R):
    if s in kept_def or s == drop: continue
    readout(s, kZ)
n_prelude_reads = sum(len(bits_of(rowsread[s])) for s in omitted)
def source(s, active):
    r = AUX + s; xr = srcop[s] - 1; fr = kvs[s] if active else kF
    move(r, fr); move(xr, fr); xor(r, xr, 'V' if active else False)
for s in early_V: source(s, True)
A1 = set(W['A']); ph1 = [i for i in range(len(ops)) if i in A1 and ops[i][0] == 'add']
ph2 = [i for i in range(len(ops)) if i not in A1 and ops[i][0] == 'add']
kept = {c: x for x, cs in W['copies'].items() for c in cs}; s0 = W['s0']; opfr = X['opfr']
KCOP = jfdata.kept_copy_ops(ops, W['copies'], s0)   # kept copies by op index, never by role pair (see jfdata)
def gate(i, active):
    _, a_, b_ = ops[i]
    if not active: fr = kF
    elif i in KCOP: fr = kvs[a_]; need(kvs[a_] == kvs[b_], 'kept copy V starts differ')   # by op index, not pair
    else: fr = kex(opfr[i])
    lab = ('kept' if i in KCOP else 'gate') if active else False
    move(AUX + a_, fr); move(AUX + b_, fr); xor(AUX + a_, AUX + b_, lab)
for i in ph1: gate(i, True)
centres = {}
for s, c in X['ret'].items():
    r = AUX + s; fr = kex(X['rootid']['ret'][s]); move(r, fr)
    tg = [v + tid[T] for T in trip if c in T]
    need(all(cur[t] == kZ for t in tg), 'centre targets not at 0')
    if A.data: need(not ((bits[r] & MASK) & ~allow(fr)), 'centre data outside U_c')
    centres[r] = tg; ev.append((2, r, fr, kZ, dim(fr))); rh['center'][dim(fr)] += 1
    for t in tg: bits[t] ^= bits[r]
for s in defer_order:
    if s in kept_def: readout(s, ksig[s])
for s in late_V: source(s, True)
for i in ph2: gate(i, True)
if A.ctl == 'late_prelude':
    readout(drop, kZ)
for s, (c, T) in X['out'].items():
    r = AUX + s; fr = kex(X['rootid']['out'][s]); t = v + tid[T]
    need(fr == ktperp[tid[T]], 'output frame != t_T^perp'); move(r, fr); move(t, fr); xor(t, r, 'out')
for i in reversed(ph1 + ph2): gate(i, False)
for s in srcop: source(s, False)
for s in range(R): move(AUX + s, kF)
for t in range(v): move(t, kF)
fwd = all(bits[i] == ((1 << i) ^ (1 << (i - v)) if v <= i < 2 * v else 1 << i) for i in range(size))
report('complete forward F2 shear on all %d formal basis vectors (y_T += x_T, x and dirty scratch restored)' % size, fwd)
report('every aux gate and read: data inside the exact frame (%d frames) %s' % (len(_allow), dict(baddata)), not baddata)
report('all X and aux registers end at F; every Y at t_T^perp', all(cur[i] == kF for i in range(v)) and
       all(cur[AUX + s] == kF for s in range(R)) and all(cur[v + t] == ktperp[t] for t in range(v)))
report('every path fully walked (no skipped or unused frame)', all(pos[i] == len(paths[i]) - 1 for i in range(size)))

# ---------------------------------------------------------------- reflected word
bank = lambda r: r + v if r < v else r - v if r < 2 * v else r
rev = [None] * size
for r in range(size): rev[bank(r)] = ~cur[r]
bits = [1 << i for i in range(size)]; rrh = Counter(); okr = True
for k, r, b_, c_, rk in reversed(ev):
    rr = bank(r)
    if k == 0:
        okr &= rev[rr] == ~c_; rev[rr] = ~b_; rrh[rk] += 1
    elif k == 1:
        sr = bank(b_); okr &= rev[rr] == rev[sr] == ~c_; bits[rr] ^= bits[sr]
    else:
        okr &= rev[rr] == ~b_
        for t in centres[r]:
            okr &= rev[bank(t)] == ~c_; bits[bank(t)] ^= bits[rr]
        rrh[rk] += 1
okr &= all(rev[bank(r)] == ~init[r] for r in range(size))
okb = all(bits[i] == ((1 << i) ^ (1 << (i + v)) if i < v else 1 << i) for i in range(size))
report('reflected word: frame continuity under complement + bank swap, returns to the complemented start', okr)
report('reflected word: complete opposite F2 shear (x_T += y_T, scratch restored)', okb)
report('reflected rank multiset == forward rank multiset', rrh == sum(rh.values(), Counter()))

# ---------------------------------------------------------------- histograms vs accounting; #144 child histogram
sys.path.insert(0, HERE)
from omit import build, D as SLOTS
w2 = SLOTS['w2']
if A.support == 'Z':
    w2 = dict(w2); w2['tg'] = {s: bits_of(rowschain[s]) for s in w2['sel']}
cb = build(w2, omitted, False)
report('event aux histogram == omit.py accounting', +rh['aux'] == +cb['aux'])
report('event target-chain histogram == omit.py accounting', +rh['tgt'] == +cb['tgt'])
report('event source-chain histogram == omit.py accounting', +rh['src'] == +cb['src'])
report('centre copies: 23 transitions of rank 22', rh['center'] == Counter({22: 23}))
m = 3 * h; Wv = 2 * v + R; H = Counter()
for part in ('aux', 'src', 'tgt', 'center'):
    for r_, n_ in rh[part].items(): H[r_] += 3 * n_
for s in kept_def: H[3 * f[s]] += 1                  # exterior of a retained gauge: m - 3(h - f) = 3f
H[2] += 2 * v
rank = sum(r_ * n_ for r_, n_ in H.items())
report('child histogram from events == omit.py build (W %d, rank %d, D %d)' % (Wv, rank, Wv * m - rank),
       dict(H) == cb['hist'] and cb['k'] == 1 and Wv == cb['W'])
report('deficit per vertex == 2v - 3h(h-1) = 2024', Wv * m - rank == 2 * v - 3 * h * (h - 1))
res = dict(OK=OK, all=all(OK.values()), omitted=len(omitted), retained=len(kept_def), support=A.support, ctl=A.ctl,
           events=len(ev), frames=len(fb), prelude_reads_moved=n_prelude_reads, W=Wv, m=m, rank=rank,
           hist={int(k): v_ for k, v_ in sorted(H.items())},
           rank_hists={k: {int(r_): n_ for r_, n_ in sorted(c.items())} for k, c in rh.items()},
           changed=sorted(Counter((dim(wkey(walk[s][1])) - f[s], dim(wkey(walk[s][1]))) for s in omitted).items()))
log('ALL', 'PASS' if res['all'] else 'FAIL ' + str([k for k, x in OK.items() if not x]))
if A.out: json.dump(res, open(os.path.join(HERE, A.out), 'w'), indent=1, default=str)
