"""Reader for the frozen round-seven witness 2 (stdlib only).

certificates/round7/jfstack_23.json.gz:
  dag    side DAG (leaves = triples of [h]; additions; outputs (c, T); retained totals; span dims);
  jf     the joint-frame-compiled program (PR #57's compiler on the COARSE=col DAG): ops, out/ret slots, role node lists;
  lift   the lifted + late-copy program: ops, the frame id of every gate (opfr), every role's chain of frame ids, the
         region frame id (Mid) and lift index j of every region, root and use frame ids, late-copy components (capdef),
         the compiled region of every op (opreg) and the destination of every use (usedest);
  flag   the integer flag matrix F (the lifted frames use its first j rows);
  sched  the final stage-1 schedule: ops with the V gates, V layer (srcop, X_S order, kept copies and their source
         role s0), the deferrable set, the frozen phase-1 op set, f = dim sigma per role, the F2 readout rows, and the
         x / y chain steps used by the accounting;
  sigma, vstart  exact deferred readout frames and V starts (canonical integer bases);
  frame_dims     dimension of every used frame id.
certificates/round7/jfstack_frames_23.json.gz: the exact canonical basis of every used frame id.
Ops: ('add', a, b) is a[a] += a[b]; ('src', s, n) is the V gate a[s] += x_S of leaf n."""
import gzip, json, os
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'certificates', 'round7')


def _load(name):
    with gzip.open(os.path.join(DATA, name), 'rt') as fh:
        return json.load(fh)


def _ops(ops): return [tuple(o) for o in ops]
def _fr(B): return tuple(tuple(r) for r in B)


def load(h=23):
    J = _load('jfstack_%d.json.gz' % h); assert J['h'] == h
    trip = list(combinations(range(h), 3)); v = len(trip)
    d = J['dag']; args = [None if a is None else tuple(a) for a in d['args']]
    sup = [0] * len(args)
    for n in sorted(d['active']):
        sup[n] = 1 << (n - 1) if args[n] is None else sup[args[n][0]] | sup[args[n][1]]
    G = dict(h=h, trip=trip, args=args, dim=d['dim'], active=set(d['active']), sup=sup,
             outputs={(c, tuple(T)): n for c, T, n in d['outputs']}, retained={c: n for c, n in d['retained']})
    jf = J['jf']
    C = dict(R=jf['R'], ops=_ops(jf['ops']), hold=jf['hold'], out={s: (c, tuple(T)) for s, c, T in jf['out']},
             ret={s: c for s, c in jf['ret']})
    L = J['lift']
    X = dict(G=G, R=L['R'], ops=_ops(L['ops']), opfr=L['opfr'], chains=L['chains'], Mid={n: i for n, i in L['Mid']},
             j={n: i for n, i in L['j']}, phys=L['phys'], fnode=L['fnode'],
             rootid={k: {p: i for p, i in mp} for k, mp in L['rootid'].items()}, dfr={u: i for u, i in L['dfr']},
             capdef={i: tuple(ids) for i, ids in L['capdef']}, out={s: (c, tuple(T)) for s, c, T in L['out']},
             ret={s: c for s, c in L['ret']})
    meta = dict(opreg=L['opreg'], usedest=[(k, (x[0], tuple(x[1])) if k == 'out' else x) for k, x in L['usedest']])
    S = J['sched']
    F2 = [sum(1 << t for t in row) for row in S['F2']]
    W = dict(R=S['R'], ops=_ops(S['ops']), srcop={s: n for s, n in S['srcop']}, copies={x: cs for x, cs in S['copies']},
             s0={x: s for x, s in S['s0']}, deferrable=S['deferrable'], A=S['phase1'], f=S['f'],
             xorder={x: o for x, o in S['xorder']}, xdata=S['xdata'], yd=S['yd'], F2=F2)
    SIG = {s: _fr(B) for s, B in J['sigma']}; VS = {s: _fr(B) for s, B in J['vstart']}
    dims = {i: dd for i, dd in J['frame_dims']}
    return dict(h=h, D=dict(G=G, C=C), X=X, meta=meta, W=W, SIG=SIG, VS=VS, dims=dims, Fl=J['flag'], claimed=J['claimed'])


def load_frames(h=23):
    """exact canonical bases EX[id] and the equality classes eqid[id] (first id, in increasing order, of the subspace)."""
    F = _load('jfstack_frames_%d.json.gz' % h); EX = {i: _fr(B) for i, B in F['EX']}
    canon = {}
    for i in sorted(EX): canon.setdefault(EX[i], i)
    return EX, {i: canon[EX[i]] for i in EX}


def kept_copy_ops(ops, copies, s0):
    """The op index of every kept V copy: the FIRST op touching the copy role, which must be ('add', c, s0[x]).

    Gotcha: do not identify a kept copy by its role pair (a in kept and b == s0[kept[a]]). Words compiled with
    skip-suffix strips and the gm point order gate the same (copy, source) pair a second time, later, at a higher
    frame; that later op is an ordinary gate and must get the walk and data-in-frame checks. Witness 2 never does
    this, so both rules agree on it (tests/test_round10.py has the truth table)."""
    kept = {c: x for x, cs in copies.items() for c in cs}
    first = {}
    for i, op in enumerate(ops):
        for r in (op[1:2] if op[0] == 'src' else op[1:3]):
            first.setdefault(r, i)
    idx = set()
    for c, x in kept.items():
        i = first.get(c)
        assert i is not None and tuple(ops[i]) == ('add', c, s0[x]), ('kept copy is not its role\'s first op', c)
        idx.add(i)
    return idx


def release_orders(K):
    """the release rule: deferred readouts by (dim sigma, id), deferred V gates by (dim F0, id), early V in X_S order."""
    W, X = K['W'], K['X']; R = W['R']; f = W['f']; srcop = W['srcop']
    sel = {s for s in range(R) if f[s] > 0}; deferrable = set(W['deferrable'])
    late_v = {s for s in srcop if s in deferrable}
    args = X['G']['args']; fnode = X['fnode']
    leafroles = {s for s in range(R) if fnode[s] is not None and args[fnode[s]] is None}
    dimF0 = {s: (len(K['VS'][s]) if s in leafroles else K['dims'][X['chains'][s][0]]) for s in range(R)}
    early_V = [s for x in sorted(W['xorder']) for s in W['xorder'][x] if s not in late_v]
    late_V = sorted(late_v, key=lambda s: (dimF0[s], s))
    defer_order = sorted(sel, key=lambda s: (f[s], s))
    return dict(sel=sel, late_v=late_v, leafroles=leafroles, dimF0=dimF0, early_V=early_V, late_V=late_V, defer_order=defer_order)


def walks(K, eqid, ro):
    """per-role frame sequence in time: [sigma_u] (deferred roles), the start (V start for leaf roles), the chain."""
    X = K['X']; R = X['R']; chains = X['chains']; out = {}
    for s in range(R):
        seq = ([('sig', s)] if s in ro['sel'] else []) + [('v', s) if s in ro['leafroles'] else eqid[chains[s][0]]] + \
              [eqid[i] for i in chains[s][1:]]
        out[s] = [k for n_, k in enumerate(seq) if n_ == 0 or k != seq[n_ - 1]]
    return out


def key_dim(K, k):
    if isinstance(k, tuple) and k[0] == 'sig': return len(K['SIG'][k[1]])
    if isinstance(k, tuple) and k[0] == 'v': return len(K['VS'][k[1]])
    return K['dims'][k]


def accounting(K, eqid):
    """chain-step ranks, deferred readout levels per target and X_S chains, rebuilt from the schedule and the frame
    dimensions (the exact frames are checked by check_frames.py and check_design.py)."""
    from collections import Counter
    h = K['h']; W = K['W']; R = W['R']; f = W['f']; ro = release_orders(K); wk = walks(K, eqid, ro)
    rk = Counter()
    for s in range(R):
        ds = ([] if s in ro['sel'] else [0]) + [key_dim(K, k) for k in wk[s]] + [h]
        assert all(b >= a for a, b in zip(ds, ds[1:])), ('decreasing chain', s)
        for a, b in zip(ds, ds[1:]):
            if b > a: rk[b - a] += 1
    v = len(K['X']['G']['trip']); lev = {t: set() for t in range(v)}
    for s in ro['defer_order']:
        g = W['F2'][s]
        while g:
            lb = g & -g; lev[lb.bit_length() - 1].add(f[s]); g ^= lb
    ylev = {t: sorted(x) for t, x in lev.items()}
    xd = []
    for x, order in sorted(W['xorder'].items()):
        early = [s for s in order if s not in ro['late_v']]
        late = sorted((s for s in order if s in ro['late_v']), key=lambda s: (len(K['VS'][s]), s))
        ds = [1] + [len(K['VS'][s]) for s in early + late] + [h]
        ds = [d for i, d in enumerate(ds) if i == 0 or d != ds[i - 1]]
        xd.append([b - a for a, b in zip(ds, ds[1:])])
    return rk, ylev, xd, ro


class Shape:
    """the fields of the accounting that deferred.histogram reads (h, v, R and f)."""
    def __init__(self, K):
        self.h = K['h']; self.v = len(K['X']['G']['trip']); self.R = K['W']['R']; self.f = K['W']['f']
