"""Our own exact implementation of the terminal-output deletion (#166 as used by #176) on a frozen paired-cube complex
word in our export_word schema (our p=7/p=8 words, or #176's p=11 word through conv176.py). Imports nothing.

Mechanism (per sink s, target group T, uniform root coefficient u, pivot c in T):
  * s is a side root role, not a source, gauge, donor or recipient; it is only ever a gate DESTINATION, every write
    has ca = 1 (z_s <- z_s + cb * ctl), and all writes lie after the centre cut;
  * the register s is deleted: no slot, no dirty column, no old-value read, no forward/inverse ops, no root read;
  * pre-shear at the cut (after the centre reads, before every deferred correction):  y_t -= y_c  (t in T \\ c);
  * each write z_s += cb * ctl becomes  y_c += u * cb * ctl  at the same time and the same op frame;
  * post-shear right after the last write:  y_t += y_c  (t in T \\ c), all of T raised to the sink's root frame U;
  * every retained role keeps its ORIGINAL old-response coefficients (computed with the sink present).
Pivot rule: c receives no correction in [cut, last write]; T receives no target event before the cut.
Checks:
  value: exact replay mod 2^61-1 with arbitrary dirty start on EVERY slot and arbitrary target start y0, 2 seeds:
         every slot restored, y - y0 = sum over all roots (sinks included) of coefficient x DAG value;
  frames: every new target chain nested in time (pivot: 0 -> write frames -> U -> later reads; others: old events,
         then U at the post-shear, then later reads), and recount of the full ledger before/after;
  cost:  delta children == -3[r] - 3[h-r] per sink, delta W == -1 per sink, deficit unchanged;
  controls (must fail): reduced-decoder compensation, post-shear before the last write, a pivot with a correction
         in its window, one write with the wrong sign, one pre-shear dropped.
Usage: python3 -I t176.py WORD.json.gz [--sinks THEIR_SINKS.json] [--out OUT.json]"""
import sys, gzip, json, random, time
from fractions import Fraction as Q
from collections import Counter, defaultdict

P = (1 << 61) - 1
t0 = time.time()
def log(*a): print('[%4.0fs]' % (time.time() - t0), *a, flush=True)
argv = sys.argv[1:]; J = json.load(gzip.open(argv[0], 'rt'))
opt = {argv[k]: argv[k + 1] for k in range(1, len(argv) - 1, 2)}
h, v, m, R = J['h'], J['v'], J['m'], J['R']; args, signs, inputs = J['args'], J['signs'], J['inputs']
roots = J['roots']; ops = [tuple(o) for o in J['ops']]; sched = J['sched']; cut = J['cut']
sources = {int(x): s for x, s in J['sources'].items()}; rootroles = J['rootroles']
gauges = J['gauges']; tau = {int(b): t for b, t in J['tau'].items()}; pairs = [tuple(q) for q in J['pairs']]
frames = [tuple(f) for f in J['frames']]
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
        if not inside(spans[xx], frames[i]): bad['span'] += 1
    def chain(sr): return [start[sr]] + [frames[i] for i in rops[sr]] + ([rootframe[sr]] if sr in rootframe else [])
    for sr in range(R):
        if sr in recip or sr in S: continue
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
    Wv = 2 * v + R - len(pairs) - len(S)
    return +C, Wv, bad, local, target

# ---------------- selection: disjoint groups, each individually frame-legal ----------------
C0, W0, bad0, loc0, tgt0 = ledger([])
check('baseline ledger legal (%s), W %d, deficit %d' % (dict(bad0), W0, W0 * m - sum(r * k for r, k in C0.items())),
      not bad0 and W0 * m - sum(r * k for r, k in C0.items()) == 2 * v - 3 * J['loss'])
used = set(); sel = []; fr_bad = Counter()
for x in cand:
    if used & set(x['T']): continue
    b = None   # only T's chains change; the joint ledger below re-checks everything
    if b is None:   # large word: test only the affected targets
        ev = target_events([x]); b = Counter()
        for t in x['T']:
            cur = FULL
            for key, A, tag in ev[t]:
                if not inside(A, cur): b[tag] += 1
                cur = A
            if not inside((inputs[t],), cur): b['input'] += 1
    if b: fr_bad.update(b); continue
    sel.append(x); used |= set(x['T'])
log('selected %d disjoint, frame-legal sinks (by |T| %s); frame rejections %s' % (len(sel),
    dict(Counter(len(x['T']) for x in sel)), dict(fr_bad)))
C1, W1, bad1, loc1, tgt1 = ledger(sel)
check('joint new word: spans, role chains, all target chains nested (%s)' % dict(bad1), not bad1)
pred = Counter(C0)
for x in sel: pred[x['r']] -= 3; pred[h - x['r']] -= 3
check('children delta = sum over sinks of -3[r] - 3[h-r]; W delta = -%d' % len(sel), +pred == C1 and W1 == W0 - len(sel))
check('deficit unchanged (%d)' % (W1 * m - sum(r * k for r, k in C1.items())),
      W1 * m - sum(r * k for r, k in C1.items()) == W0 * m - sum(r * k for r, k in C0.items()))

# ---------------- exact replay ----------------
def replay(sel, seed, ctl=None):
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
    for b in gauged: rt[tau[b]].append(b)
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
        for b in rt.get(k, ()): read(b, -1, cf[b])
        d, c, _ = ops[i]; ca, cb = coef[i]
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
    for b in rt.get(len(sched), ()): read(b, -1, cf[b])
    for j, (r, sr) in enumerate(zip(roots, rootroles)):
        if r['kind'] == 'side' and sr not in S: read(sr, +1, {t: md(c) for t, c in zip(r['targets'], r['coefficients'])})
    for dd, cc, ca, cb in reversed(chron): a[dd] = (a[dd] - cb * a[cc]) * pow(ca % P, P - 2, P) % P
    for xx, sr in sources.items(): a[slot[sr]] = (a[slot[sr]] - x[xx - 1]) % P
    val = dag(x); want = list(y0)
    for r in roots:
        for t, c in zip(r['targets'], r['coefficients']): want[t] = (want[t] + md(c) * val[r['node']]) % P
    return a == z, y == want

for sd in (1, 2):
    check('replay seed %d, %d sinks deleted, %d slots, dirty targets: slots restored and every target exact'
          % (sd, len(sel), len(set(sr for sr in range(R) if sr not in recip)) - len(sel)), replay(sel, sd) == (True, True))
if sel:
    for c_ in ('reduced_decoder', 'post_early', 'wrong_sign', 'drop_pre'):
        r_ = replay(sel, 11, c_); check('control %s rejected %s' % (c_, r_), r_ != (True, True))
    # bad pivot: a selected sink re-pivoted on a target that gets a correction inside the window
    bp = next((x for x in sel if x['dirty']), None)
    if bp is not None: alt = [dict(x) for x in sel if x['s'] != bp['s']] + [dict(bp, c=bp['dirty'][0])]
    else:
        bp = next((x for x in cand if x['dirty']), None)
        if bp is not None: alt = [dict(bp, c=bp['dirty'][0])]
    if bp is not None:
        r_ = replay(alt, 12); check('control pivot with a correction in its window rejected %s' % (r_,), r_ != (True, True))
    else: log('no sink with a corrected target: bad-pivot control not applicable')

res = dict(word=J['name'], p=J['p'], eligible=len(cand), selected=len(sel), rejected=dict(why), frame_rej=dict(fr_bad),
           sinks=[dict(role=x['s'], root=x['j'], targets=x['T'], pivot=x['c'], r=x['r'], writes=rops[x['s']]) for x in sel],
           C0={str(k): c for k, c in sorted(C0.items())}, W0=W0, C1={str(k): c for k, c in sorted(C1.items())}, W1=W1,
           local1={str(k): c for k, c in sorted(loc1.items())}, target1={str(k): c for k, c in sorted(tgt1.items())},
           ok=OK)
if '--sinks' in opt:   # compare with #176's frozen selection (data only)
    TS = json.load(open(opt['--sinks']))
    theirs = {z['role'] for z in TS['sinks']}; ours = {x['s'] for x in cand}
    log('their %d sinks; in our eligible set: %d; our eligible not theirs: %d' % (len(theirs), len(theirs & ours), len(ours - theirs)))
    sel2 = [x for x in cand if x['s'] in theirs]
    bad = [z['role'] for z in TS['sinks'] if z['role'] not in ours]
    check('every #176 sink is eligible under our rule (missing %s)' % bad, not bad)
    C2, W2, bad2, loc2, tgt2 = ledger(sel2)
    check('#176 selection: all chains nested (%s)' % dict(bad2), not bad2)
    norm = lambda d: {int(k): int(c) for k, c in d.items() if int(c) and int(k)}   # rank-0 steps cost nothing
    check('#176 selection: recount == their local histogram', norm(loc2) == norm(TS['local_histogram']))
    check('#176 selection: recount == their target histogram', norm(tgt2) == norm(TS['target_histogram']))
    check('#176 selection: recount == their child histogram, W %d == %d' % (W2, TS['new_W']),
          norm(C2) == norm(TS['child_histogram']) and W2 == TS['new_W'])
    for sd in (4, 5):
        check('#176 selection replay seed %d' % sd, replay(sel2, sd) == (True, True))
    res['theirs'] = dict(n=len(theirs), C={str(k): c for k, c in sorted(C2.items())}, W=W2)
if '--out' in opt: json.dump(res, open(opt['--out'], 'w'), indent=1)
log('RESULT t176 %s p=%d: %d sinks, %s' % (J['name'], J['p'], len(sel), 'ALL PASS' if all(OK.values()) else 'FAIL'))
sys.exit(0 if all(OK.values()) else 1)
