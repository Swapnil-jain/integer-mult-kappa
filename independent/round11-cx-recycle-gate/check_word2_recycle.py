"""Independent checker of a FROZEN paired-cube complex word with birth-read reuse AND terminal-output deletions
(#166's dirty-target terminal lemma as used by #176; audited SOUND in audit-176). Reads only the word JSON written by
export_word.py (+ its 'sinks' list); imports no construction code. Built from audit-176/t176.py (selection replaced
by the frozen list) and check_word.py (decoder identity, birth checks, certificate). Checks:
  1. decoder identity (H + K) x = x mod P on random x (every root, sinks included, is part of H);
  2. every frozen sink is eligible (side root role, not source/gauge/donor/recipient, destination-only additive
     writes after the centre cut, uniform root coefficient, no target event before the cut, pivot clean in
     [cut, last write]) and the sinks' target groups are disjoint;
  3. legality on F2 subspaces of the NEW word: spans, spliced role chains (sinks dropped), birth pairs (chronology,
     read before first touch, donor last frame inside the gauge), and every target chain in time order including the
     pre-shears, the pivot writes at their op frames and the post-shears to the root frame;
  4. ledger: recount == claim; children delta == sum of -3[r] - 3[h-r] per sink; W delta == -1 per sink; deficit;
  5. the 1e-12 certificate with the next point rejected;
  6. exact replay mod P, 2 seeds, arbitrary dirty slots AND dirty targets (shears, pivot writes, ORIGINAL
     compensation columns), and controls that must fail: birth (omitted recipient read, read after first touch, read
     before donor's last write) and terminal (reduced-decoder compensation, post-shear before the last write, a wrong
     write sign, a dropped pre-shear, a pivot with a correction inside its window).
Usage: python3 check_word2.py WORD.json.gz"""

import sys, gzip, json, random, time
from fractions import Fraction as Q
from collections import Counter, defaultdict

P = (1 << 61) - 1
t0 = time.time()
def log(*a): print('[%4.0fs]' % (time.time() - t0), *a, flush=True)
J = json.load(gzip.open(sys.argv[1], 'rt'))
h, v, m, R = J['h'], J['v'], J['m'], J['R']; args, signs, inputs = J['args'], J['signs'], J['inputs']
roots = J['roots']; ops = [tuple(o) for o in J['ops']]; sched = J['sched']; cut = J['cut']
sources = {int(x): s for x, s in J['sources'].items()}; rootroles = J['rootroles']
gauges = J['gauges']; tau = {int(b): t for b, t in J['tau'].items()}; pairs = [tuple(q) for q in J['pairs']]
frames = [tuple(f) for f in J['frames']]
# ---- collective hosting (acct): a fresh birth B is born on a dead copy pair's zero-signal difference.
# Mix op (pivot, control, -1) runs pivot -= control at B's first-op frame, just before B's first op; B is renamed to
# the pivot; no new read (the dirty response is recomputed from the new transcript and read at frame 0).
HOSTS = J.get('hosts', []); MIX = set(); UNUSED = set()
if HOSTS:
    ren = {z['birth']: z['pivot'] for z in HOSTS}; assert len(ren) == len(HOSTS)
    assert not {z[k] for z in HOSTS for k in ('pivot', 'control')} & set(J['rootroles']), 'a recycled copy is a root role'
    assert len({z[k] for z in HOSTS for k in ('pivot', 'control', 'birth')}) == 3 * len(HOSTS), 'hosts overlap'
    assert not set(ren) & {z['role'] for z in J.get('sinks', [])} and not set(ren) & set(int(b) for b in J['tau'])
    firstk = {}
    for k, i in enumerate(sched):
        for r_ in ops[i][:2]:
            if r_ in ren and r_ not in firstk: firstk[r_] = k
    ins = defaultdict(list)
    for z in HOSTS:
        k = firstk[z['birth']]; j = len(ops); ops.append((z['pivot'], z['control'], -1)); frames.append(frames[sched[k]])
        ins[k].append(j); MIX.add(j)
    ops = [(ren.get(d, d), ren.get(c, c), x) for d, c, x in ops]
    newsched = []; fmap = {}
    for k, i in enumerate(sched):
        fmap[k] = len(newsched); newsched += ins.get(k, []); newsched.append(i)
    fmap[len(sched)] = len(newsched)
    sched = newsched; cut = fmap[cut]
    J['tau'] = {b: fmap[t] for b, t in J['tau'].items()}
    J['rootroles'] = [ren.get(s_, s_) for s_ in J['rootroles']]
    J['pairs'] = [[ren.get(a_, a_), ren.get(b_, b_)] for a_, b_ in J['pairs']]
    UNUSED = set(ren)
rootroles = J['rootroles']; tau = {int(b): t for b, t in J['tau'].items()}; pairs = [tuple(q) for q in J['pairs']]
md = lambda q: Q(q).numerator % P * pow(Q(q).denominator % P, P - 2, P) % P
OK = {}
def check(name, ok): OK[name] = bool(ok); log(('PASS ' if ok else 'FAIL ') + name)

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
    piv = {y.bit_length() - 1: y for y in basis(B)}
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
n = len(args)
def dag(x):
    val = [0] * n
    for k in range(1, n): val[k] = x[k - 1] if args[k] is None else (val[args[k][0]] + signs[k] * val[args[k][1]]) % P
    return val
spans = [()] * n
for x in range(1, n): spans[x] = (inputs[x - 1],) if args[x] is None else basis(spans[args[x][0]] + spans[args[x][1]])

pos = {i: k for k, i in enumerate(sched)}
rops = [[] for _ in range(R)]
for i in sched: a, b, _ = ops[i]; rops[a].append(i); rops[b].append(i)
hold = {sr: xx for xx, sr in sources.items()}; coef = {}
for i in sched:
    d, c, xx = ops[i]
    if i in MIX:
        assert hold.get(d) is not None and hold[d] == hold[c], 'mix of unequal values'; hold[d] = None; coef[i] = (1, -1); continue
    if hold.get(d) is None: assert hold[c] == xx; coef[i] = (1, 1)
    else:
        aa, bb = args[xx]; sg = signs[xx]
        if hold[d] == aa: assert hold[c] == bb; coef[i] = (1, sg)
        else: assert hold[d] == bb and hold[c] == aa; coef[i] = (sg, 1)
    hold[d] = xx
rootframe, rootkind, rootann, rootidx = {}, {}, {}, {}
for j, (r, sr) in enumerate(zip(roots, rootroles)):
    rootidx[sr] = j; rootkind[sr] = r['kind']
    if r['kind'] == 'center': rootframe[sr] = spans[r['node']]
    else: rootann[sr] = basis([inputs[t] for t in r['targets']]); rootframe[sr] = perp(rootann[sr])
start = [()] * R
for xx, sr in sources.items(): start[sr] = (inputs[xx - 1],)
gsig = {}
for z in gauges: gsig[z['role']] = perp(z['A']); start[z['role']] = gsig[z['role']]
gauged = set(gsig); donor = {a: b for a, b in pairs}; recip = {b: a for a, b in pairs}

def compensation(drop=()):
    """old-response coefficients cf[role] (future coefficient on the targets), backward over the word.
    drop: sinks treated as absent (the WRONG reduced decoder, used only as a control)."""
    seeds = [(sr, r['kind'], {t: md(c) for t, c in zip(r['targets'], r['coefficients'])}) for r, sr in zip(roots, rootroles)
             if sr not in drop]
    cf = [dict() for _ in range(R)]
    def addto(dst, src, k):
        for t, u in src.items(): dst[t] = (dst.get(t, 0) + k * u) % P
    for sr, kind, sd in seeds:
        if kind == 'side': addto(cf[sr], sd, 1)
    for k in range(len(sched) - 1, -1, -1):
        i = sched[k]; d, c, _ = ops[i]; ca, cb = coef[i]
        if k == cut - 1:
            for sr, kind, sd in seeds:
                if kind == 'center': addto(cf[sr], sd, 1)
        if d in drop: continue
        addto(cf[c], cf[d], cb)
        if ca != 1: cf[d] = {t: ca * u % P for t, u in cf[d].items()}
    if cut == 0:
        for sr, kind, sd in seeds:
            if kind == 'center': addto(cf[sr], sd, 1)
    return cf
CF = compensation()

# ---------------- sink eligibility and pivots ----------------
reads_on = defaultdict(list)            # target -> [(tau, gauge role)]
for z in gauges:
    for t in z['targets']: reads_on[t].append((tau[z['role']], z['role']))
srcroles = set(sources.values()); why = Counter(); cand = []
for j, (r, sr) in enumerate(zip(roots, rootroles)):
    if r['kind'] != 'side': continue
    T = list(r['targets'])
    if sr in srcroles: why['source'] += 1; continue
    if sr in gauged or sr in donor or sr in recip: why['gauge/alias'] += 1; continue
    if len(set(r['coefficients'])) != 1: why['non-uniform root'] += 1; continue
    if not rops[sr]: why['no ops'] += 1; continue
    if any(ops[i][0] != sr for i in rops[sr]): why['controls another op'] += 1; continue
    if any(coef[i][0] != 1 for i in rops[sr]): why['non-additive write'] += 1; continue
    if any(pos[i] < cut for i in rops[sr]): why['centre-phase write'] += 1; continue
    L = pos[rops[sr][-1]]
    if any(tt < cut for t in T for tt, _ in reads_on[t]): why['target event before cut'] += 1; continue
    clean = [t for t in T if not any(cut <= tt <= L for tt, _ in reads_on[t])]
    dirty = [t for t in T if t not in clean]
    if not clean: why['no clean pivot'] += 1; continue
    cand.append(dict(s=sr, j=j, T=T, c=clean[0], dirty=dirty, L=L, u=md(r['coefficients'][0]), r=len(rootframe[sr])))
log('side roots %d; eligible %d by |T| %s; rejected %s' % (sum(r['kind'] == 'side' for r in roots), len(cand),
    dict(Counter(len(x['T']) for x in cand)), dict(why)))

# ---------------- target-chain legality for one sink (all others unchanged) ----------------
def target_events(sel):
    """per target, the ordered list of (key, annihilator, tag) of the NEW word."""
    S = {x['s']: x for x in sel}; ev = defaultdict(list)
    for k, z in enumerate(gauges):   # same-clock gauge reads: reverse list order (the exporter's convention)
        A = basis(z['A'])
        for t in z['targets']: ev[t].append(((tau[z['role']], 0, len(gauges) - 1 - k), A, 'g'))
    for x in sel:
        if len(x['T']) > 1:
            for t in x['T']: ev[t].append(((cut, -1, 0), FULL, 'pre'))
        for i in rops[x['s']]: ev[x['c']].append(((pos[i], 1, 0), perp(frames[i]), 'w'))
        for t in x['T']: ev[t].append(((x['L'], 2, 0), rootann[x['s']], 'post'))
    for sr, j in rootidx.items():
        if rootkind[sr] == 'side' and sr not in S:
            for t in roots[j]['targets']: ev[t].append(((float('inf'), j, 0), rootann[sr], 'root'))
    for t in ev: ev[t].sort(key=lambda e: e[0])
    return ev

def ledger(sel):
    """full ledger of the word with sinks `sel` deleted (our check_word normalisation); returns (C, W, bad)."""
    S = {x['s'] for x in sel}; bad = Counter(); local = Counter()
    for i, (_, _, xx) in enumerate(ops):
        if i not in MIX and not inside(spans[xx], frames[i]): bad['span'] += 1
    def chain(sr): return [start[sr]] + [frames[i] for i in rops[sr]] + ([rootframe[sr]] if sr in rootframe else [])
    for sr in range(R):
        if sr in recip or sr in S or sr in UNUSED: continue
        seq = []; cur = sr
        while cur is not None: seq += chain(cur); cur = donor.get(cur)
        seq.append(FULL)
        for A_, B_ in zip(seq, seq[1:]):
            if not inside(A_, B_): bad['role nesting'] += 1
        if len(start[sr]) == 1 and sr not in gsig: local[1] += 1
        d = [len(X) for X in seq]
        for d0, d1 in zip(d, d[1:]):
            if d1 > d0: local[d1 - d0] += 1
        cur = sr
        while cur is not None:
            if rootkind.get(cur) == 'center': local[len(rootframe[cur])] += 1
            cur = donor.get(cur)
    ev = target_events(sel); target = Counter()
    for t in range(v):
        cur = FULL
        for key, A, tag in ev.get(t, []):
            if not inside(A, cur): bad['target nesting (%s)' % tag] += 1
            target[len(cur) - len(A)] += 1; cur = A
        if not inside((inputs[t],), cur): bad['target input'] += 1
        target[len(cur) - 1] += 1
    C = Counter()
    for hist in (local, {1: v, 2: v, h - 4: v}, target):
        for r_, k in hist.items():
            if r_: C[r_] += 3 * k
    C[2] += 2 * v
    for z in gauges:
        if z['role'] not in recip: C[3 * (h - len(basis(z['A'])))] += 1
    Wv = 2 * v + R - len(pairs) - len(S) - len(UNUSED)
    return +C, Wv, bad, local, target

# ---------------- selection: disjoint groups, each individually frame-legal ----------------

# ---------------- exact replay ----------------
def replay(sel, seed, ctl=None, omit=None, after_first=None, before_death=None):
    rnd = random.Random(seed); S = {x['s']: x for x in sel}
    cf = compensation(drop=set(S)) if ctl == 'reduced_decoder' else CF
    multi = [x for x in sel if len(x['T']) > 1 and len(rops[x['s']]) > 1]
    tx = multi[0]['s'] if multi else (sel[0]['s'] if sel else None)
    def rs(sr):
        while sr in recip: sr = recip[sr]
        return sr
    slot = {sr: rs(sr) for sr in range(R) if sr not in S}; phys = sorted(set(slot.values()))
    z = {q: rnd.randrange(P) for q in phys}; a = dict(z)
    x = [rnd.randrange(P) for _ in range(v)]; y0 = [rnd.randrange(P) for _ in range(v)]; y = list(y0)
    def read(sr, sg, vec):
        val = a[slot[sr]]
        for t, u in vec.items(): y[t] = (y[t] + sg * u * val) % P
    rt = defaultdict(list)
    for b in gauged:
        t_ = tau[b]
        if b == after_first: t_ = pos[rops[b][0]] + 1
        if b == before_death: t_ = max(pos[i] for i in rops[recip[b]] if ops[i][0] == recip[b])
        rt[t_].append(b)
    lastw = {}
    for s, xx in S.items():
        L = rops[s][-1]
        if ctl == 'post_early' and s == tx: lastw[rops[s][-2] if len(rops[s]) > 1 else None] = s
        else: lastw[L] = s
    pre_skip = tx if ctl == 'drop_pre' else None
    flip = rops[tx][0] if ctl == 'wrong_sign' else None
    for sr in range(R):
        if sr not in gauged and sr not in S: read(sr, -1, cf[sr])
    for xx, sr in sources.items(): a[slot[sr]] = (a[slot[sr]] + x[xx - 1]) % P
    def shear(xx, sg):
        for t in xx['T']:
            if t != xx['c']: y[t] = (y[t] + sg * y[xx['c']]) % P
    chron = []
    def centre_and_pre():
        for j, (r, sr) in enumerate(zip(roots, rootroles)):
            if r['kind'] == 'center': read(sr, +1, {t: md(c) for t, c in zip(r['targets'], r['coefficients'])})
        for s, xx in S.items():
            if s != pre_skip and len(xx['T']) > 1: shear(xx, -1)
    for k, i in enumerate(sched):
        if k == cut: centre_and_pre()
        for b in rt.get(k, ()):
            if b != omit: read(b, -1, cf[b])
        d, c, _ = ops[i]; ca, cb = coef[i]
        if i in MIX and ctl in ('host_skip_mix', 'host_plus') and i == min(MIX):
            if ctl == 'host_skip_mix': continue
            cb = 1
        if d in S:
            xx = S[d]; sg = -1 if i == flip else 1
            y[xx['c']] = (y[xx['c']] + sg * xx['u'] * cb * a[slot[c]]) % P
        else:
            dd, cc = slot[d], slot[c]; assert dd != cc
            a[dd] = (ca * a[dd] + cb * a[cc]) % P; chron.append((dd, cc, ca, cb))
        if i in lastw:
            xx = S[lastw[i]]
            if len(xx['T']) > 1: shear(xx, +1)
    if cut == len(sched): centre_and_pre()
    for b in rt.get(len(sched), ()):
        if b != omit: read(b, -1, cf[b])
    for j, (r, sr) in enumerate(zip(roots, rootroles)):
        if r['kind'] == 'side' and sr not in S: read(sr, +1, {t: md(c) for t, c in zip(r['targets'], r['coefficients'])})
    for dd, cc, ca, cb in reversed(chron): a[dd] = (a[dd] - cb * a[cc]) * pow(ca % P, P - 2, P) % P
    for xx, sr in sources.items(): a[slot[sr]] = (a[slot[sr]] - x[xx - 1]) % P
    val = dag(x); want = list(y0)
    for r in roots:
        for t, c in zip(r['targets'], r['coefficients']): want[t] = (want[t] + md(c) * val[r['node']]) % P
    return a == z, y == want


# ---------------- 1. decoder identity ----------------
rng = random.Random(7); x = [rng.randrange(P) for _ in range(v)]; val = dag(x); y = [0] * v
for r in roots:
    for t, c in zip(r['targets'], r['coefficients']): y[t] = (y[t] + md(c) * val[r['node']]) % P
half = md(Q(1, 2)); nbad = 0
for t in range(v):
    k0 = t - t % 8; acc = 0
    for s_ in range(k0, k0 + 8):
        w_ = bin(inputs[t] ^ inputs[s_]).count('1'); acc += x[s_] if w_ == 6 else -x[s_] if w_ == 2 else 0
    nbad += (y[t] + acc * half - x[t]) % P != 0
check('decoder identity H x + K x = x on random x mod P (%d wrong of %d)' % (nbad, v), nbad == 0)

# ---------------- 2. frozen sinks ----------------
byrole = {x_['s']: x_ for x_ in cand}; sel = []; okS = True; usedT = set()
for z in J.get('sinks', []):
    c_ = byrole.get(z['role'])
    if c_ is None or z['pivot'] not in c_['T'] or z['pivot'] in c_['dirty'] or sorted(z['targets']) != sorted(c_['T']):
        okS = False; continue
    if usedT & set(c_['T']): okS = False
    usedT |= set(c_['T']); sel.append(dict(c_, c=z['pivot']))
check('%d frozen sinks: each eligible with a clean pivot, target groups disjoint' % len(J.get('sinks', [])), okS and len(sel) == len(J.get('sinks', [])))

# ---------------- 3. legality (new word) and birth pairs ----------------
bad = Counter()
check('pairs are a matching (%d)' % len(pairs), len(donor) == len(recip) == len(pairs))
if HOSTS:   # recycling chronology: both copies retired before the mix; the control is never touched again
    badh = Counter(); birthops = {}
    for z in HOSTS:
        j = next(i for i in MIX if ops[i][0] == z['pivot'] and ops[i][1] == z['control'])
        k = pos[j]
        if any(pos[i] > k for i in rops[z['control']] if i != j): badh['control used after its mix'] += 1
        pre = [i for i in rops[z['pivot']] if pos[i] < k]
        if not pre: badh['pivot never written'] += 1
        if z['control'] in rootframe or z['pivot'] in gauged or z['control'] in gauged: badh['root/gauge copy'] += 1
        if J['ops'][0] is None: pass
    check('recycling chronology and kinds (%s)' % dict(badh), not badh)
for a_, b_ in pairs:
    if a_ in rootframe or not rops[a_]: bad['donor kind'] += 1
    if b_ not in gsig: bad['recipient kind'] += 1
    if not pos[rops[a_][-1]] < tau[b_]: bad['chronology'] += 1
    if rops[b_] and not tau[b_] <= pos[rops[b_][0]]: bad['read after first touch'] += 1
    if not inside(frames[rops[a_][-1]], gsig[b_]): bad['donor frame outside gauge'] += 1
for z in gauges:
    b_ = z['role']
    if rops[b_] and not tau[b_] <= pos[rops[b_][0]]: bad['gauge read after first touch'] += 1
check('birth pairs and gauge clocks legal (%s)' % dict(bad), not bad)
C0, W0, bad0, loc0, tgt0 = ledger([])
C1, W1, bad1, loc1, tgt1 = ledger(sel)
check('new word: spans, role chains, all target chains nested incl. shears and pivot writes (%s)' % dict(bad1), not bad1)

# ---------------- 4. ledger ----------------
pred = Counter(C0)
for x_ in sel: pred[x_['r']] -= 3; pred[h - x_['r']] -= 3
check('children delta = sum over sinks of -3[r] - 3[h-r]; W delta = -%d' % len(sel), +pred == C1 and W1 == W0 - len(sel))
s1 = sum(r_ * k for r_, k in C1.items())
check('telescoping deficit W m - s = 2v - 3 loss (%d)' % (W1 * m - s1), W1 * m - s1 == 2 * v - 3 * J['loss'])
check('every child narrower than m (max %d)' % max(C1), max(C1) < m)
CL = J['claim_hosts'] if HOSTS else J['claim']
claimC = {int(r_): Q(k) for r_, k in CL['C'].items()}
check('child histogram and W (recounted from chains) equal the exporter claim', {r_: Q(k) for r_, k in C1.items()} == claimC and Q(W1) == Q(CL['W']))
if HOSTS:
    import math
    def Ff(a_): return math.fsum(k * (r_ / m) ** (1 - a_) for r_, k in C1.items()) / W1
    lo_, hi_ = 0.0, 0.05
    for _ in range(200):
        mid = (lo_ + hi_) / 2; lo_, hi_ = (mid, hi_) if Ff(mid) < 1 else (lo_, mid)
    J['claim']['a_c'] = '%d/%d' % (int(lo_ * 10 ** 12), 10 ** 12)

# ---------------- 5. certificate ----------------
def ln_upper(z, N=28):
    z = Q(z); k = 0
    while z > 2: z /= 2; k += 1
    def ta(u): return 2 * (sum(u ** (2 * j + 1) / (2 * j + 1) for j in range(N)) + u ** (2 * N + 1) / ((2 * N + 1) * (1 - u * u)))
    return k * ta(Q(1, 3)) + (ta((z - 1) / (z + 1)) if z > 1 else 0)
lu = {w: ln_upper(Q(m, w)) for w in C1}
def Fup(a):
    tot = Q(0)
    for w, k in C1.items():
        u = a * lu[w]; assert 0 <= u < 1
        tot += k * Q(w, m) * (1 + u + u * u / (2 * (1 - u / 3)))
    return tot / W1
a = Q(J['claim']['a_c']); den = 10 ** 12
if HOSTS:   # certify the new word's own 1e-12 grid floor (rigorous Fup), starting from the float root
    k_ = int(Q(J['claim']['a_c']) * den) + 3
    while not Fup(Q(k_, den)) < 1: k_ -= 1
    a = Q(k_, den)
check('a_c = %s: grid point accepted, denominator divides 1e12' % a, (a * den).denominator == 1 and Fup(a) < 1)
check('next grid point rejected', not Fup(a + Q(1, den)) < 1)

# ---------------- 6. replay and controls ----------------
for sd in (1, 2):
    check('replay seed %d (%d sinks, dirty slots and targets): slots restored, every target exact' % (sd, len(sel)), replay(sel, sd) == (True, True))
if pairs:
    check('control omitted recipient read rejected', replay(sel, 3, omit=pairs[0][1]) != (True, True))
    b_ = next(b for _, b in pairs if rops[b]); check('control read after first touch rejected', replay(sel, 4, after_first=b_) != (True, True))
    b_ = next(b for a_, b in reversed(pairs) if any(ops[i][0] == a_ for i in rops[a_]))
    check('control read before donor last write rejected', replay(sel, 5, before_death=b_) != (True, True))
if sel:
    for c_ in ('reduced_decoder', 'post_early', 'wrong_sign', 'drop_pre'):
        check('control %s rejected' % c_, replay(sel, 11, c_) != (True, True))
    bp = next((x_ for x_ in cand if x_['dirty']), None)
    if bp is not None:
        alt = [x_ for x_ in sel if not set(x_['T']) & set(bp['T'])] + [dict(bp, c=bp['dirty'][0])]
        check('control pivot with a correction in its window rejected', replay(alt, 12) != (True, True))
if HOSTS:
    for c_ in ('host_skip_mix', 'host_plus'):
        check('control %s rejected' % c_, replay(sel, 13, c_) != (True, True))
log('RESULT check_word2 %s p=%d a_c=%s sinks=%d: %s' % (J['name'], J['p'], a, len(sel), 'ALL PASS' if all(OK.values()) else 'FAIL'))
sys.exit(0 if all(OK.values()) else 1)
