"""Round-seven witness 2, part B (stdlib): the compiled programs and the explicit two-phase stage-1 schedule.

 A. side DAG semantics from the definitions (disjoint additions; outputs = leave-two-out sums; retained = star(c)).
 B. jf program (COARSE=col compile, frozen): EXACT symbolic F2 replay of the plain word with dirty
    scratch (every scratch bit a separate symbol): pass-1 readouts, V, L, J, L^-1, V^-1; and the same for the lifted +
    late-copy program. Targets must receive exactly x_T with zero garbage; every role restored.
 C. phase 1: own read/write closure of the retained totals' last writes; equal to the frozen phase-1 set; the
    two-phase order preserves the version of every role seen by every op (so L is unchanged as a map); deferred
    roles untouched in phase 1; every retained total completes in phase 1; ops after phase 1 on a retained role
    only read it and run at a frame EQUAL to its root U_c (exact frames, part A).
 D. explicit time-ordered stage-1 word (pass-1 readouts at 0, early V, phase 1, centre reads, deferred readouts
    sorted by dim sigma, deferred V gates sorted by dim F0, phase 2, output reads), EXACT symbolic F2 replay:
    y_T += x_T exactly, no garbage, scratch restored; data at every read = its definition.
 E. frames in TIME order on every role: every event (gate, V gate, readout, centre read, output read) at a frame of
    the role's chain, pointer monotone; kept V copies at the shared V start; data of both roles of every gate inside
    the gate frame (t_S in frame, exact integer dot products); deferred roles hold zero data at their readout.
 F. negative controls of the word and of the walk.
The walk is the ACTUAL reordered op sequence (phase 1, centre reads, deferred readouts, deferred V gates, phase 2,
output reads), not the original op order.
Usage: python3 check_schedule.py"""
import time, sys, os, random
from collections import defaultdict, Counter
T0 = time.time()
def log(*a): print('[%5.0fs]' % (time.time() - T0), *a, flush=True)
OK = {}
def report(name, ok):
    OK[name] = bool(ok); log('%-104s %s' % (name, 'PASS' if ok else 'FAIL'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jfdata
from xq import ann
K = jfdata.load(); EXa, eqida = jfdata.load_frames()
D = K['D']; X = K['X']; W = K['W']; E = dict(EX=EXa, eqid=eqida)
h = 23; G = D['G']; trip = [tuple(t) for t in G['trip']]; v = len(trip); args = G['args']; act = sorted(G['active'])
tid = {T: i for i, T in enumerate(trip)}; R = X['R']

# ------------------------------------------------------------------ A. DAG
star = [sum(1 << i for i, T in enumerate(trip) if q in T) for q in range(h)]
sup = {}; okA = True
for n in act:
    if args[n] is None: okA &= 1 <= n <= v; sup[n] = 1 << (n - 1)
    else:
        a, b = args[n]; okA &= not (sup[a] & sup[b]); sup[n] = sup[a] | sup[b]
for (c, T), n in G['outputs'].items():
    A_, B_ = [q for q in T if q != c]; okA &= c in T and sup[n] == star[c] & ~star[A_] & ~star[B_]
okA &= len(G['outputs']) == 3 * v and sorted(G['retained']) == list(range(h))
okA &= all(sup[n] == star[c] for c, n in G['retained'].items())
okA &= all(sup[n] == G['sup'][n] for n in act)
report('A. COARSE=col DAG: disjoint additions, out(c,T) = {S: S cap T = {c}}, retained = star(c) (%d outputs, %d nodes)' % (len(G['outputs']), len(act)), okA)

# ------------------------------------------------------------------ B. plain words, exact symbolic F2
def plain_word(ops, out, ret, Rr, srcs):
    """symbolic plain word. role value = data bits [0, v) + scratch-symbol bits [v, v+R). Returns ok, details."""
    a = [1 << (v + s) for s in range(Rr)]; a0 = list(a)
    for i, op in enumerate(ops):                       # forward L with V folded in (src op = the V gate x_S)
        if op[0] == 'src': a[op[1]] ^= 1 << (op[2] - 1)
        else: a[op[1]] ^= a[op[2]]
    y = [0] * v; okd = True
    for s, (c, T) in out.items():
        n = G['outputs'][c, T]; okd &= (a[s] & ((1 << v) - 1)) == sup[n]; y[tid[T]] ^= a[s]
    for s, c in ret.items():
        okd &= (a[s] & ((1 << v) - 1)) == star[c]
        for T in trip:
            if c in T: y[tid[T]] ^= a[s]
    # garbage coefficients = scratch symbols reaching y_T; the pass-1 readouts subtract exactly these
    gar = [yy >> v for yy in y]
    for t in range(v): y[t] ^= gar[t] << v
    okx = all(y[t] == 1 << t for t in range(v))
    for op in reversed(ops):
        if op[0] == 'src': a[op[1]] ^= 1 << (op[2] - 1)
        else: a[op[1]] ^= a[op[2]]
    return okd and okx and a == a0, gar
C = D['C']
ok1, garJ = plain_word(C['ops'], C['out'], C['ret'], C['R'], None)
report('B. jf program: symbolic F2 plain word, dirty scratch: data at reads = definitions, y_T += x_T, restored', ok1)
ok2, garL = plain_word(X['ops'], X['out'], X['ret'], R, None)
report('B. lifted + late-copy program: same symbolic F2 plain word', ok2)
ok3, garV = plain_word(W['ops'], X['out'], X['ret'], R, None)
report('B. final L with the V-leaf gates (frozen sched ops): same symbolic F2 plain word', ok3)
F2 = W['F2']
coefL = [0] * R
for t, g in enumerate(garV):
    while g:
        lb = g & -g; s = lb.bit_length() - 1; coefL[s] |= 1 << t; g ^= lb
report('B. own garbage coefficients (symbolic) == frozen F2 readout rows (all %d roles)' % R, coefL == list(F2))
del garJ, garL, garV

# ------------------------------------------------------------------ C. phase 1
ops = W['ops']; srcop = W['srcop']; f = W['f']; sel = {s for s in range(R) if f[s] > 0}
ret = X['ret']; out = X['out']
lastw = {}; readers = defaultdict(list); dep = defaultdict(list)
for i, op in enumerate(ops):
    if op[0] != 'add': continue
    a_, b_ = op[1], op[2]
    if b_ in lastw: dep[i].append(lastw[b_])
    if a_ in lastw: dep[i].append(lastw[a_])
    dep[i].extend(readers[a_]); readers[a_] = []; readers[b_].append(i); lastw[a_] = i
def closure(extra):
    A1 = set(); st = [lastw[s] for s in ret]
    while st:
        i = st.pop()
        if i not in A1: A1.add(i); st.extend(dep[i]); st.extend(extra.get(i, ()))
    return A1
# frame-order edges: an op touching role b at chain position k must follow every earlier op touching b at a lower
# position (positions from the chain walk in ORIGINAL op order, frames compared as exact subspaces)
eqid0 = E['eqid']; pos0 = {}; pp = [0] * R
for i, op in enumerate(ops):
    if op[0] != 'add': continue
    for s_ in (op[1], op[2]):
        ch = [eqid0[c] for c in X['chains'][s_]]; k = pp[s_]
        while k < len(ch) and ch[k] != eqid0[X['opfr'][i]]: k += 1
        assert k < len(ch); pp[s_] = k; pos0[i, s_] = k
bylev = defaultdict(lambda: defaultdict(list))
for (i, s_), k in pos0.items(): bylev[s_][k].append(i)
fo = defaultdict(list)
for s_, d_ in bylev.items():
    ks = sorted(d_)
    for a_, b_ in zip(ks, ks[1:]):
        for i in d_[b_]: fo[i].extend(d_[a_])
Arw = closure({}); Afo = closure(fo)
log('   read/write closure %d ops; with frame-order edges %d ops; frozen phase 1 %d ops' % (len(Arw), len(Afo), len(W['A'])))
A1 = set(W['A'])
report('C. frozen phase 1 == own closure (%s)' % ('read/write + frame order' if A1 == Afo else 'read/write only' if A1 == Arw else 'neither'), A1 in (Afo, Arw))
ph1 = [i for i in range(len(ops)) if i in A1]; ph2 = [i for i in range(len(ops)) if i not in A1]
def versions(order):
    cnt = Counter(); seen = {}
    for i in order:
        op = ops[i]
        if op[0] != 'add': continue
        seen[i] = (cnt[op[1]], cnt[op[2]]); cnt[op[1]] += 1
    return seen
report('C. two-phase order: every op sees the same version of both its roles (L unchanged as a map)', versions(range(len(ops))) == versions(ph1 + ph2))
touch1 = {s for i in ph1 for s in ops[i][1:3]}
report('C. no phase-1 op touches a deferred role (%d deferred, f > 0)' % len(sel), not (touch1 & sel))
deferrable = set(W['deferrable']); late_v = {s for s in srcop if s in deferrable}
report('C. deferred V gates (%d) on roles untouched in phase 1; deferred set inside the deferrable set' % len(late_v), not (touch1 & late_v) and sel <= deferrable)
report('C. every retained total completes (last write) in phase 1', all(lastw[s] in A1 for s in ret))
eqid = E['eqid']; rootret = X['rootid']['ret']; rootout = X['rootid']['out']; opfr = X['opfr']
post = [(i, s) for i in ph2 for s in ops[i][1:3] if ops[i][0] == 'add' and s in ret]
okpost = all(ops[i][2] == s and eqid[opfr[i]] == eqid[rootret[s]] for i, s in post)
report('C. %d ops after phase 1 on retained roles: all read-only, at a frame equal (exactly) to U_c' % len(post), okpost)

# ------------------------------------------------------------------ D. explicit time-ordered word
xorder = W['xorder']; copies = W['copies']; s0 = W['s0']
dimF0 = {}
FRd = None
chains = X['chains']; fnode = X['fnode']
EX = E['EX']
leafroles = {s for s in range(R) if fnode[s] is not None and args[fnode[s]] is None}
for s in range(R):
    dimF0[s] = len(K['VS'][s]) if s in leafroles else len(EX[chains[s][0]])
early_V = [s for x in sorted(xorder) for s in xorder[x] if s not in late_v]
late_V = sorted(late_v, key=lambda s: (dimF0[s], s))
defer_order = sorted(sel, key=lambda s: (f[s], s))
# the frozen X_S order must agree with the time order (early ones in xorder order, then deferred by dim F0)
okxo = True
for x, o in xorder.items():
    e_ = [s for s in o if s not in late_v]; l_ = [s for s in o if s in late_v]
    okxo &= o == e_ + l_ and [dimF0[s] for s in l_] == sorted(dimF0[s] for s in l_) and [dimF0[s] for s in e_] == sorted(dimF0[s] for s in e_)
report('D. frozen X_S order == time order (early V first, deferred V sorted by dim F0)', okxo)
ro = jfdata.release_orders(K)
report('D. release orders of the certificate (jfdata) == these', (ro['sel'], ro['early_V'], ro['late_V'], ro['defer_order']) == (sel, early_V, late_V, defer_order))
MASK = (1 << v) - 1
rows = {s: [t for t in range(v) if F2[s] >> t & 1] for s in range(R)}
def word(ctl=None):
    a = [1 << (v + s) for s in range(R)]; a0 = list(a); y = [0] * v; okd = True
    drop = None; dsel = set(sel)
    if ctl == 'drop': drop = defer_order[len(defer_order) // 2]
    if ctl == 'phase1_role':
        dsel.add(next(s for s in range(R) if s not in deferrable and F2[s] and any(ops[i][0] == 'add' and ops[i][1] == s for i in ph1)))
    def readout(s):
        nonlocal okd
        if s == drop: return
        if s in sel and a[s] & MASK: okd = False                     # a deferred role holds zero data at its readout
        for t in rows[s]: y[t] ^= a[s]
    for s in range(R):
        if s not in dsel: readout(s)
    for s in early_V: a[s] ^= 1 << (srcop[s] - 1)
    if ctl == 'V_before_readout':
        for s in late_V: a[s] ^= 1 << (srcop[s] - 1)
    def run(idx):
        for i in idx:
            op = ops[i]
            if op[0] == 'add': a[op[1]] ^= a[op[2]]
    run(ph1)
    for s, c in ret.items():
        okd &= (a[s] & MASK) == star[c]
        for T in trip:
            if c in T: y[tid[T]] ^= a[s]
    if ctl != 'late_readout':
        for s in sorted(dsel, key=lambda s: (f[s], s)): readout(s)
    if ctl != 'V_before_readout':
        for s in late_V: a[s] ^= 1 << (srcop[s] - 1)
    run(ph2)
    if ctl == 'late_readout':
        for s in defer_order: readout(s)
    for s, (c, T) in out.items():
        okd &= (a[s] & MASK) == sup[G['outputs'][c, T]]; y[tid[T]] ^= a[s]
    for i in reversed(ph1 + ph2):
        op = ops[i]
        if op[0] == 'add': a[op[1]] ^= a[op[2]]
    for s, x in srcop.items(): a[s] ^= 1 << (x - 1)
    return okd and a == a0 and all(y[t] == 1 << t for t in range(v))
report('D. EXACT symbolic F2 stage-1 word on the two-phase schedule: y_T += x_T, zero garbage, scratch restored', word())
for ctl in ('drop', 'late_readout', 'V_before_readout', 'phase1_role'):
    report('F. NEG word control %-17s breaks the identity' % ctl, not word(ctl))

# ------------------------------------------------------------------ E. frames in time order
startkey = {}
for s in range(R): startkey[s] = ('v', s) if s in leafroles else eqid[chains[s][0]]
walk = {}
for s in range(R):
    seq = ([('sig', s)] if s in sel else []) + [startkey[s]] + [eqid[i] for i in chains[s][1:]]
    walk[s] = [k for n_, k in enumerate(seq) if n_ == 0 or k != seq[n_ - 1]]
ptr = [-1] * R; bad = Counter()
def at(s, key, what):
    ch = walk[s]; k = max(ptr[s], 0)
    while k < len(ch) and ch[k] != key: k += 1
    if k == len(ch): bad[what] += 1
    else: ptr[s] = k
kept = {c: x for x, cs in copies.items() for c in cs}
for s in range(R):
    if s not in sel and ptr[s] != -1: bad['pass1'] += 1
for s in early_V: at(s, ('v', s), 'earlyV')
def gate(i):
    op = ops[i]
    if op[0] != 'add': return
    a_, b_ = op[1], op[2]
    if a_ in kept and b_ == s0[kept[a_]]:                        # kept V copy: at the shared V start
        at(a_, ('v', a_), 'keptcopy'); at(b_, ('v', b_), 'keptcopy'); return
    k = eqid[opfr[i]]; at(a_, k, 'gate'); at(b_, k, 'gate')
for i in ph1: gate(i)
for s in ret: at(s, eqid[rootret[s]], 'centre')
for s in defer_order:
    if ptr[s] != -1: bad['deferred touched before readout'] += 1
    at(s, ('sig', s), 'readout')
for s in late_V: at(s, ('v', s), 'lateV')
for i in ph2: gate(i)
for s in out: at(s, eqid[rootout[s]], 'output')
report('E. every event in TIME order at a frame of both roles\' chains, pointers monotone %s' % (dict(bad) or ''), not bad)
okkc = all(K['VS'][c] == K['VS'][s0[x]] for c, x in kept.items())
report('E. kept V copies (%d) start at the V start of their source role (same frame)' % len(kept), okkc)
# negative control: the same walk with phase 2 before phase 1 must fail somewhere
ptr = [-1] * R; bad = Counter()
for i in ph2 + ph1: gate(i)
report('F. NEG walk with the phases swapped is detected (%d violations)' % sum(bad.values()), sum(bad.values()) > 0)

# data inside the gate frame (exact): t_S in frame <=> Z t_S = 0 for the integer annihilator Z; t_S has three ones,
# so each product is the exact integer z_a + z_b + z_c
allowed = {}
def allow(k):
    if k not in allowed:
        B = EX[k]
        if len(B) == h: allowed[k] = MASK; return MASK
        Z_ = ann(B, h)
        allowed[k] = sum(1 << t for t, T in enumerate(trip) if all(z[T[0]] + z[T[1]] + z[T[2]] == 0 for z in Z_))
    return allowed[k]
a = [0] * R; badd = 0
for s in early_V: a[s] ^= 1 << (srcop[s] - 1)
def dgate(i):
    global badd
    op = ops[i]
    if op[0] != 'add': return
    a_, b_ = op[1], op[2]
    if a_ in kept and b_ == s0[kept[a_]]: a[a_] ^= a[b_]; return            # V-start membership is checked in part C
    m = allow(eqid[opfr[i]]); badd += bool((a[a_] | a[b_]) & ~m); a[a_] ^= a[b_]; badd += bool(a[a_] & ~m)
for i in ph1: dgate(i)
for s, c in ret.items(): badd += bool(a[s] & ~allow(eqid[rootret[s]]))
for s in late_V: a[s] ^= 1 << (srcop[s] - 1)
for i in ph2: dgate(i)
for s in out: badd += bool(a[s] & ~allow(eqid[rootout[s]]))
report('E. data of both roles of every gate (and of every read) inside the frame, exactly (%d frames tested)' % len(allowed), badd == 0)
report('E. this walk == the walk of the certificate (jfdata.walks)', walk == jfdata.walks(K, eqid, ro))
log('ALL %s' % ('PASS' if all(OK.values()) else 'FAIL: ' + str([k for k, x in OK.items() if not x])))
sys.exit(0 if all(OK.values()) else 1)
