"""Frame descent on the cube word (the cover analogue of #126/#129 physical operation frames / our shrunk lift).

#144 places every node x at its CAP U_x = perp(ann[x]) (the largest frame inside everything it reaches). Any nested
assignment span(x) <= U_x with U_operand <= U_x <= U_successor (successors: consuming ops, carried-arc targets, side
root caps; copied-centre nodes stay exactly at their star span) is a legal word with the same roles, deficit and
data/target histograms; only the aux ledger changes. The ledger is written as signed terms f(linear form in ranks)
(f(r) = r ln(m/r), the cover's local Lambda), and coordinate descent moves each U_x between
lo_x = sum of operand frames and hi_x = intersection of successor frames, including intermediate dimensions
(lo + a prefix of a fixed complement inside hi). Output: new annihilators (ann' = perp U), used by select_gauges and
validated by recount (nesting on actual subspaces)."""
import math, sys
from collections import defaultdict, Counter
sys.path.insert(0, '.')
from cube_word import basis, perp, contained

def terms_of(prof, wit):
    """ledger terms: list of (weight, {node: coeff}, const) so that the internal histogram is
    sum_w w * [const + sum coeff*rank]; exactly #144's ledger (compile_frames)."""
    g = wit['g']; h = prof['h']; args = g['args']; roots = g['roots']; order = wit['order']; uses = wit['uses']
    uval = wit['uval']; utgt = wit['utgt']; rframe = wit['rframe']; n = wit['n']; arcs = wit['arcs']
    T = []
    for x in order:
        k = len(uses[x]) - 1
        if k: T.append((k, {x: 1}, 0))
        if args[x]:
            T.append((1, {x: -1}, h))
            for y in args[x]: T.append((1, {x: 1, y: -1} if x != y else {}, 0))
        else: T.append((1, {}, 1)); T.append((1, {x: 1}, -1))
    for j, root in enumerate(roots):
        x = root['node']; rt = len(rframe[j])
        if root['kind'] == 'center': T.append((1, {x: 1}, 0)); T.append((1, {x: -1}, h))
        else: T.append((1, {x: -1}, rt)); T.append((1, {}, h - rt))
    for x, u in arcs.items():
        val = uval[u]; t = utgt[u]
        tl = ({t: 1}, 0) if t < n else ({}, len(rframe[t - n]))
        T.append((-1, {x: -1}, h)); T.append((-1, {val: 1}, 0))
        c = dict(tl[0]); c[val] = c.get(val, 0) - 1; T.append((-1, c, tl[1]))
        c = dict(tl[0]); c[x] = c.get(x, 0) - 1; T.append((1, c, tl[1]))
    return T

def run(prof, wit, sweeps=6, verbose=True, mode='full'):
    g = wit['g']; h = prof['h']; m = 3 * h; args = g['args']; roots = g['roots']; order = wit['order']; n = wit['n']
    inputs = g['inputs']; uval = wit['uval']; utgt = wit['utgt']; rframe = wit['rframe']; arcs = wit['arcs']
    f = lambda r: r * math.log(m / r) if r > 0 else (0.0 if r == 0 else float('inf'))
    U = {x: perp(wit['ann'][x], h) for x in order}
    rank = {x: len(U[x]) for x in order}
    T = terms_of(prof, wit); idx = defaultdict(list)
    for i, (w, c, k) in enumerate(T):
        for x in c: idx[x].append(i)
    def val(i, over=None):
        w, c, k = T[i]; s = k
        for x, a in c.items(): s += a * (over[1] if over and x == over[0] else rank[x])
        return w * f(s)
    lam = sum(val(i) for i in range(len(T)))
    # structure: successors (consumers, arc targets, root frames), centre-fixed nodes
    succ = defaultdict(list); rootfr = defaultdict(list); fixed = set()
    for x in order:
        if args[x]:
            for y in args[x]: succ[y].append(x)
    for j, r in enumerate(roots):
        if r['kind'] == 'center': fixed.add(r['node'])
        else: rootfr[r['node']].append(rframe[j])
    carried_in = defaultdict(list)   # donors whose carried register enters op t: U_donor <= U_t (part of t's lower bound)
    for x, u in arcs.items():
        t = utgt[u]
        if t < n: succ[x].append(t); carried_in[t].append(x)
        else: rootfr[x].append(rframe[t - n])
    spans = {}
    for x in order:
        spans[x] = (inputs[x - 1],) if not args[x] else basis(spans[args[x][0]] + spans[args[x][1]])
    if verbose: print('descent start Lambda_int %.1f' % lam, flush=True)
    for sw in range(sweeps):
        moved = 0; gain = 0.0
        seq = list(reversed(order)) if sw % 2 == 0 else list(order)
        for x in seq:
            if x in fixed: continue
            lo = basis(sum((U[y] for y in args[x]), ())) if args[x] else (inputs[x - 1],)
            lo = basis(lo + spans[x] + tuple(z for d in carried_in[x] for z in U[d]))
            cons = [perp(U[s], h) for s in succ[x]] + [perp(F, h) for F in rootfr[x]]
            hi = perp(basis(tuple(z for c in cons for z in c)), h) if cons else tuple(1 << j for j in range(h - 1, -1, -1))
            if not contained(lo, hi): continue
            # complement of lo inside hi (fixed order): candidate frames lo + comp[:k]
            comp = []; cur = list(lo)
            for z in hi:
                if len(basis(cur + [z])) > len(basis(cur)): comp.append(z); cur.append(z)
            ks = range(len(comp) + 1) if mode == 'full' else sorted({0, len(comp)})
            base = sum(val(i) for i in idx[x]); best = (0.0, None)
            for k in ks:
                r = len(lo) + k
                if r == rank[x]: continue
                d = sum(val(i, (x, r)) for i in idx[x]) - base
                if d < best[0] - 1e-9: best = (d, k)
            if best[1] is not None:
                k = best[1]; U[x] = basis(lo + tuple(comp[:k])); rank[x] = len(U[x]); moved += 1; gain += best[0]
        lam += gain
        if verbose: print('sweep %d moved %d gain %.1f Lambda_int %.1f' % (sw, moved, gain, lam), flush=True)
        if not moved: break
    ann = list(wit['ann'])
    for x in order: ann[x] = perp(U[x], h)
    H = Counter()
    for w, c, k in T: H[k + sum(a * rank[x] for x, a in c.items())] += w
    assert min(H.values()) >= 0, 'negative ledger count'
    return ann, rank, +H

def apply(prof, wit, ann, rank, H):
    wit2 = dict(wit, ann=ann, rank=[rank.get(x, 0) if isinstance(rank, dict) else rank[x] for x in range(wit['n'])])
    prof2 = dict(prof, H=H)
    return prof2, wit2
