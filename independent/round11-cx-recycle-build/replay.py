"""Exact scalar replay (mod P = 2^61-1) of the physical cube word with birth-read reuse, arbitrary dirty scratch.

Independent of the ledger: it uses only the DAG (args, signs, roots), the op list (dest, ctl, node), the chronological
schedule, the gauge read clocks and the pairs. Every physical slot starts with a random value z. Each role's
old-value read subtracts (its future coefficient) x (slot value) at its read time: ordinary roles at time 0, gauged
roles at their read clock (merged recipients read their donor's leftover). Root reads add coefficient x value
(centres at the cut, sides at the end). Cleanup undoes every op in reverse. Checks: every slot returns to z, and
y - y0 equals sum over roots of coefficient x DAG value of the root node, computed straight from the DAG.
Controls: omit one recipient's read; read a recipient after its first op; read a recipient before its donor dies."""
import sys, os, random
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
P = (1 << 61) - 1
def md(q): q = Q(q); return q.numerator % P * pow(q.denominator % P, P - 2, P) % P

def prepare(W, F):
    G = W['wit']['g']; s = W['s']; args = G['args']; inputs = G['inputs']; v = len(inputs); R = W['pr']['R']
    sign = lambda x: G['signs'][x - 1]
    ops = s['ops']; sched = F['sched']; cut = F['cut']
    # op coefficients from the operand each dest register holds (tracked through the schedule)
    hold = {}
    for x, sr in s['sources'].items(): hold[sr] = x
    coef = {}
    for i in sched:
        d, c, x = ops[i]
        if hold.get(d) is None:            # fan-out copy: fresh role d += c
            assert hold[c] == x; coef[i] = (1, 1)
        else:
            aa, bb = args[x]; sg = sign(x)
            if hold[d] == aa: assert hold[c] == bb; coef[i] = (1, sg)
            else: assert hold[d] == bb and hold[c] == aa; coef[i] = (sg, 1)
        hold[d] = x
    roots = G['roots']; rr = s['rootroles']
    seeds = []
    for r, sr in zip(roots, rr):
        seeds.append((sr, r['kind'], {t: md(q) for t, q in zip(r['targets'], r['coefficients'])}, r['node']))
    # future coefficients: backward over the chronological events (ops; centre reads at the cut; side reads at end)
    cf = [dict() for _ in range(R)]; birth = [None] * R
    def addto(dst, src, k):
        for t, u in src.items(): dst[t] = (dst.get(t, 0) + k * u) % P
    for sr, kind, sd, _ in seeds:
        if kind == 'side': addto(cf[sr], sd, 1)
    for k in range(len(sched) - 1, -1, -1):
        i = sched[k]; d, c, x = ops[i]; ca, cb = coef[i]
        # centre reads happen just after the last phase-1 op, i.e. before sched[cut]
        if k == cut - 1:
            for sr, kind, sd, _ in seeds:
                if kind == 'center': addto(cf[sr], sd, 1)
        addto(cf[c], cf[d], cb)
        if ca != 1: cf[d] = {t: ca * u % P for t, u in cf[d].items()}
    if cut == 0:
        for sr, kind, sd, _ in seeds:
            if kind == 'center': addto(cf[sr], sd, 1)
    return dict(coef=coef, seeds=seeds, cf=cf, v=v, R=R, args=args, sign=sign, inputs=inputs)

def dag_values(Pd, x):
    args = Pd['args']; val = [0] * len(args)
    for n in range(1, len(args)):
        if args[n] is None: val[n] = x[n - 1]
        else: a, b = args[n]; val[n] = (val[a] + Pd['sign'](n) * val[b]) % P
    return val

def run(W, F, Pd, pairs, tau, seed, omit=None, after_first=None, before_death=None):
    """returns (slots restored, y correct)."""
    rng = random.Random(seed); R = Pd['R']; v = Pd['v']; ops = W['s']['ops']; sched = F['sched']; cut = F['cut']
    gauged = {z['role'] for z in W['s']['sel']}
    merge = {b: a for a, b in pairs}
    def root_slot(sr):
        while sr in merge: sr = merge[sr]
        return sr
    slot = {sr: root_slot(sr) for sr in range(R)}
    phys = sorted(set(slot.values()))
    z = {q: rng.randrange(P) for q in phys}; a = dict(z)
    x = [rng.randrange(P) for _ in range(v)]; y0 = [rng.randrange(P) for _ in range(v)]; y = list(y0)
    cf = Pd['cf']; coef = Pd['coef']
    def read(sr, sgn, cvec):
        val = a[slot[sr]]
        for t, u in cvec.items(): y[t] = (y[t] + sgn * u * val) % P
    rt = {}
    for b in gauged:
        t = tau[b]
        if b == after_first: t = F['pos'][F['rops'][b][0]] + 1
        if b == before_death:    # before the donor's last WRITE (its last op is as control, which leaves the slot unchanged)
            t = max(F['pos'][i] for i in F['rops'][merge[b]] if ops[i][0] == merge[b])
        rt.setdefault(t, []).append(b)
    for sr in range(R):
        if sr not in gauged: read(sr, -1, cf[sr])
    for xn, sr in W['s']['sources'].items(): a[slot[sr]] = (a[slot[sr]] + x[xn - 1]) % P
    chron = []
    for k, i in enumerate(sched):
        if k == cut:
            for sr, kind, sd, _ in Pd['seeds']:
                if kind == 'center': read(sr, +1, sd)
        for b in rt.get(k, ()):
            if b != omit: read(b, -1, cf[b])
        d, c, _ = ops[i]; ca, cb = coef[i]; dd, cc = slot[d], slot[c]
        assert dd != cc, 'aliased gate ports'
        a[dd] = (ca * a[dd] + cb * a[cc]) % P; chron.append((dd, cc, ca, cb))
    if cut == len(sched):
        for sr, kind, sd, _ in Pd['seeds']:
            if kind == 'center': read(sr, +1, sd)
    for b in rt.get(len(sched), ()):
        if b != omit: read(b, -1, cf[b])
    for sr, kind, sd, _ in Pd['seeds']:
        if kind == 'side': read(sr, +1, sd)
    for dd, cc, ca, cb in reversed(chron):
        a[dd] = (a[dd] - cb * a[cc]) * pow(ca % P, P - 2, P) % P
    for xn, sr in W['s']['sources'].items(): a[slot[sr]] = (a[slot[sr]] - x[xn - 1]) % P
    val = dag_values(Pd, x); want = list(y0)
    for sr, kind, sd, node in Pd['seeds']:
        for t, u in sd.items(): want[t] = (want[t] + u * val[node]) % P
    return a == z, y == want

def validate(W, F, pairs, tau, log=print):
    Pd = prepare(W, F); res = {}
    for sd in (1, 2): res['seed%d' % sd] = run(W, F, Pd, pairs, tau, sd)
    log('replay with %d pairs (%d slots): %s' % (len(pairs), Pd['R'] - len(pairs), res))
    ok = all(r == (True, True) for r in res.values())
    ctl = {}
    if pairs:
        b = pairs[0][1]; ctl['omit recipient read'] = run(W, F, Pd, pairs, tau, 3, omit=b)
        b = next((b for _, b in pairs if F['rops'][b]), None)
        if b is not None: ctl['read after first touch'] = run(W, F, Pd, pairs, tau, 4, after_first=b)
        ops = W['s']['ops']   # a donor whose last op writes its slot (as ctl the slot is unchanged, nothing to detect)
        b = next(b for a_, b in reversed(pairs) if any(ops[i][0] == a_ for i in F['rops'][a_])); ctl['read before donor death'] = run(W, F, Pd, pairs, tau, 5, before_death=b)
    log('controls (must fail): %s' % ctl)
    cok = all(r != (True, True) for r in ctl.values())
    return ok, cok, res, ctl
