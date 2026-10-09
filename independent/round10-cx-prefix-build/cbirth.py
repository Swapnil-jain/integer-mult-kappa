"""Birth-read reuse on our paired-cube complex word (cube-stack's descent + NDS words), own code.

Chronological schedule: the word's ops (select_gauges' op list, a topological order) run as phase 1 (the centre
cone, in op order) then the rest (in op order); position[i] is an op's clock. Every gauged role B is untouched until
its first op, so its old-value read may run at any clock <= position[first op of B] (late read, #143 style). The
read order on each target must follow decreasing annihilators (reverse selection order), so deadlines are clamped
by late_tau (a later-selected gauge on a shared target never reads after an earlier-selected one).

A pair (A, B): B is born on the slot of a dead donor A (not a root role, last op before B's read, last frame
E_A inside sigma_B). A's chain becomes ... E_A -> sigma_B -> B's chain -> F, and B's own slot disappears:
  * no NDS:  W -1, B's tail child 3d removed, three (h - e_A) children become three (d_B - e_A);
  * NDS:     W -3/t_B, B's tail (m - t_B u_B at weight 3/t_B) removed, same donor change.
Rank and deficit are unchanged. A recipient may later donate (chains) if it is not a root role.

The cost is recounted from the actual role chains (F2 subspaces, nesting checked), not from a formula.
"""
import sys, os, pickle, math, time
from collections import Counter, defaultdict
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cube_word import build, compile_frames, select_gauges, children, basis, perp, contained
from matching import maxcard
from descent import terms_of
from gtypes import classify
from evalc import acert, D_of
T0 = time.time()
def log(*a): print('[%5.0fs]' % (time.time() - T0), *a, flush=True)

def ex(t, a, m): return t * math.expm1(a * math.log(m / t)) if t > 0 else 0.0

def load_word(p, pk=None, sel_override=None, gauge_kw=None):
    """rebuild cube-stack's best_p word: own max-card matching, descended annihilators from the pickle."""
    g = build(p); prof, wit = compile_frames(g, lambda D: maxcard(D)[0]); h = 2 * p
    B = pickle.load(open(pk or os.path.join(HERE, 'best_p%d.pkl' % p), 'rb'))
    ann = B['ann']; n = wit['n']
    rank = [h - len(ann[x]) if ann[x] is not None else 0 for x in range(n)]
    H = Counter()
    for w, c, k in terms_of(prof, wit): H[k + sum(a * rank[x] for x, a in c.items())] += w
    prof2 = dict(prof, H=+H); wit2 = dict(wit, ann=ann, rank=rank)
    pr, s = select_gauges(prof2, wit2, **(gauge_kw or {}))
    return dict(p=p, h=h, m=3 * h, g=g, prof=prof2, wit=wit2, pr=pr, s=s, best=B)

def gfirst_frames(W):
    """op frames, chronology, role op lists, chains' fixed ends."""
    wit, s, h = W['wit'], W['s'], W['h']; g = W['g']; ann = wit['ann']
    ops = s['ops']; ph = sorted(s['phase']); rest = [i for i in range(len(ops)) if i not in s['phase']]
    sched = ph + rest; pos = {i: k for k, i in enumerate(sched)}
    fr = {}; frames = []
    for (_, _, x) in ops:
        if x not in fr: fr[x] = perp(ann[x], h)
        frames.append(fr[x])
    R = W['pr']['R']; rops = [[] for _ in range(R)]
    for i in sched:
        a, b, _ = ops[i]; rops[a].append(i); rops[b].append(i)
    roots = g['roots']; inputs = g['inputs']
    # compile_frames shifts node ids by one when node 0 is a real input; s['ops'] / roots use the shifted ids
    G = wit['g']; inputs = G['inputs']; roots = G['roots']
    rootframe, rootkind = {}, {}
    for r, sr in zip(roots, s['rootroles']):
        if r['kind'] == 'center':
            args = G['args']; x = r['node']
            U = spans_of(G)[x]; rootframe[sr] = U; rootkind[sr] = 'center'
        else:
            rootframe[sr] = perp(basis(inputs[t] for t in r['targets']), h); rootkind[sr] = 'side'
    start = [()] * R
    for x, sr in s['sources'].items(): start[sr] = (inputs[x - 1],)
    sel = s['sel']; gsig = {}
    for z in sel: gsig[z['role']] = perp(z['A'], h); start[z['role']] = gsig[z['role']]
    return dict(ops=ops, sched=sched, pos=pos, frames=frames, rops=rops, rootframe=rootframe, rootkind=rootkind,
                start=start, gsig=gsig, cut=len(ph), phase=set(ph))

_SP = {}
def spans_of(G):
    # The cache is keyed by id(G), and Python reuses the id of a freed graph. A process that scores many graphs
    # (an annealer or a module search) could otherwise be handed the previous graph's spans, so each entry also
    # keeps the graph itself and is reused only when it is the same object.
    k = id(G)
    if k not in _SP or _SP[k][0] is not G:
        args = G['args']; inputs = G['inputs']; n = len(args); sp = [()] * n
        for x in range(1, n):
            sp[x] = (inputs[x - 1],) if args[x] is None else basis(sp[args[x][0]] + sp[args[x][1]])
        _SP[k] = (G, sp)
    return _SP[k][1]

def late_tau(W, F, strict_only=True):
    """latest legal read clock of every gauged role: <= its first op's clock, and target order kept. On a target the
    gauges' annihilators grow along selection order, so a later-selected gauge with a STRICTLY larger annihilator
    must not read after an earlier one; equal annihilators may read in any order (the step between them is 0)."""
    sel = W['s']['sel']; pos = F['pos']; rops = F['rops']; tau = {}
    lastA, gmin, pmin = {}, {}, {}
    INF = 10 ** 9
    for z in sel:
        sr = z['role']; A = tuple(z['A']); d = pos[rops[sr][0]] if rops[sr] else INF; t = d
        for x in z['targets']:
            if x not in lastA: continue
            if strict_only and lastA[x] == A: t = min(t, pmin[x])
            else: t = min(t, pmin[x], gmin[x])
        tau[sr] = t
        for x in z['targets']:
            if x in lastA and strict_only and lastA[x] == A: gmin[x] = min(gmin[x], t)
            else:
                pmin[x] = min(pmin.get(x, INF), gmin.get(x, INF)); gmin[x] = t; lastA[x] = A
    return tau

def nds_t(z, h, m):
    u, r, alt = classify(z['A'], h)
    return 3 if alt else max(3, m // (u + r)), u

def recount(W, F, pairs, tau, nds=True, check=True):
    """physical recount from role chains with hand-offs; returns (children Counter of Q, W per vertex, info)."""
    h, m, pr, s = W['h'], W['m'], W['pr'], W['s']; R = pr['R']; v = pr['v']
    donor = {a: b for a, b in pairs}; recip = {b: a for a, b in pairs}
    assert len(donor) == len(pairs) == len(recip)
    frames, rops, start, rootframe, rootkind, pos = F['frames'], F['rops'], F['start'], F['rootframe'], F['rootkind'], F['pos']
    full = tuple(1 << j for j in range(h - 1, -1, -1))
    sp = spans_of(W['wit']['g']); bad = Counter()
    for i, (_, _, x) in enumerate(F['ops']):
        if not contained(sp[x], frames[i]): bad['span'] += 1
    for a, b in pairs:
        if a in rootframe or not rops[a]: bad['donor kind'] += 1
        if b not in F['gsig']: bad['recipient kind'] += 1
        if not pos[rops[a][-1]] < tau[b]: bad['chronology'] += 1
        if rops[b] and not tau[b] <= pos[rops[b][0]]: bad['read after first touch'] += 1
        if not contained(frames[rops[a][-1]], F['gsig'][b]): bad['donor frame outside sigma'] += 1
    def chain(sr):
        return [start[sr]] + [frames[i] for i in rops[sr]] + ([rootframe[sr]] if sr in rootframe else [])
    local = Counter(); nsteps = 0
    for sr in range(R):
        if sr in recip: continue
        seq = []; cur = sr
        while cur is not None:
            seq += chain(cur); cur = donor.get(cur)
        seq.append(full)
        if check:
            for A_, B_ in zip(seq, seq[1:]):
                if not contained(A_, B_): bad['nesting'] += 1
        if start[sr] and len(start[sr]) == 1 and sr not in F['gsig']: local[1] += 1   # source load 0 -> <q>
        dims = [len(X) for X in seq]
        for d0, d1 in zip(dims, dims[1:]):
            if d1 > d0: local[d1 - d0] += 1
        cur = sr
        while cur is not None:
            if rootkind.get(cur) == 'center': local[len(rootframe[cur])] += 1   # copied-centre copy (loss)
            cur = donor.get(cur)
    # targets along the actual read clocks (sel reversed = read order within one clock)
    sel = s['sel']; G = W['wit']['g']; inputs = G['inputs']
    when = {z['role']: (tau[z['role']], k) for k, z in enumerate(reversed(sel))}
    cur = [full] * v; target = Counter()
    for z in sorted(sel, key=lambda z: when[z['role']]):
        A = tuple(z['A'])
        for t in z['targets']:
            if not contained(A, cur[t]): bad['target order'] += 1
            target[len(cur[t]) - len(A)] += 1; cur[t] = A
    for r, sr in zip(G['roots'], s['rootroles']):
        if r['kind'] == 'side':
            A = basis(inputs[t] for t in r['targets'])
            for t in r['targets']:
                if not contained(A, cur[t]): bad['side order'] += 1
                target[len(cur[t]) - len(A)] += 1; cur[t] = A
    for t, q in enumerate(inputs):
        if not contained((q,), cur[t]): bad['input'] += 1
        target[len(cur[t]) - 1] += 1
    C = Counter()
    for hist in (local, pr['source'], target):
        for r, n in hist.items():
            if r: C[r] += Q(3 * n)
    C[2] += 2 * v
    Wv = Q(2 * v + R - len(pairs)); deg = Counter()
    for z in sel:
        if z['role'] in recip: continue
        d = z['d']
        if nds:
            t, u = nds_t(z, h, m); deg[t] += 1
            Wv += Q(3, t) - 1
            if m - t * u: C[m - t * u] += Q(3, t)
        else: C[3 * d] += 1
    return +C, Wv, dict(bad=dict(bad), local=local, target=target, deg=deg)

def candidates(W, F, tau, a0, nds=True, chain_ok=True, kdim=int(os.environ.get('KDIM', 120)), extra=None):
    """compatible (recipient, donor, gain) edges; gain = decrease of sum ex at a0 (D is unchanged). The gain depends
    only on the donor's last-frame dim e, so each recipient lists up to kdim timely donors per dim e (earliest
    deaths first, rotated by recipient to spread the lists)."""
    h, m = W['h'], W['m']; R = W['pr']['R']; rops, frames, pos = F['rops'], F['frames'], F['pos']
    gsig = F['gsig']; sel = W['s']['sel'] if extra is None else extra
    elig = [sr for sr in range(R) if rops[sr] and sr not in F['rootframe'] and (chain_ok or sr not in gsig)]
    byE = defaultdict(list)
    for sr in elig: byE[frames[rops[sr][-1]]].append(sr)
    for E in byE: byE[E].sort(key=lambda sr: pos[rops[sr][-1]])
    byA = defaultdict(list)
    for z in sel: byA[tuple(z['A'])].append(z)
    edges = []
    for A, zs in byA.items():
        comp = defaultdict(list)
        for E in byE:
            if all(bin(a & f).count('1') % 2 == 0 for a in A for f in E): comp[len(E)].extend(byE[E])
        for e in comp: comp[e].sort(key=lambda sr: pos[rops[sr][-1]])
        for k_, z in enumerate(zs):
            B = z['role']; d = z['d']
            if nds:
                t, u = nds_t(z, h, m); tg = 3 / t * ex(m - t * u, a0, m)
            else: tg = ex(3 * d, a0, m)
            for e, L in comp.items():
                gain = 3 * (ex(h - e, a0, m) - ex(d - e, a0, m)) + tg
                if gain <= 0: continue
                ok = [sr for sr in L if pos[rops[sr][-1]] < tau[B] and sr != B]
                if len(ok) > kdim:
                    off = (k_ * kdim) % len(ok); ok = (ok[off:] + ok[:off])[:kdim]
                edges.extend((B, sr, gain) for sr in ok)
    return edges

def match(edges):
    import numpy as np
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import min_weight_full_bipartite_matching
    if not edges: return []
    recs = sorted({b for b, _, _ in edges}); dons = sorted({a for _, a, _ in edges})
    ri = {b: i for i, b in enumerate(recs)}; di = {a: j for j, a in enumerate(dons)}
    gmax = max(g for _, _, g in edges)
    r_ = [ri[b] for b, _, _ in edges] + list(range(len(recs)))
    c_ = [di[a] for _, a, _ in edges] + [len(dons) + i for i in range(len(recs))]
    v_ = [gmax - g + 1.0 for _, _, g in edges] + [gmax + 1.0] * len(recs)
    M = csr_matrix((v_, (r_, c_)), shape=(len(recs), len(dons) + len(recs)))
    rr, cc = min_weight_full_bipartite_matching(M)
    return [(dons[c], recs[r]) for r, c in zip(rr, cc) if c < len(dons)]

def evaluate(W, C, Wv):
    a, rf = acert(C, Wv, W['m']); return a

if __name__ == '__main__':
    p = int(sys.argv[1]); a0 = float(sys.argv[2]) if len(sys.argv) > 2 else 5.5e-4
    W = load_word(p); B = W['best']; pr = W['pr']
    log('p=%d rebuilt: R %d gauges %d  sel==pickle %s' % (p, pr['R'], len(W['s']['sel']),
        [z['role'] for z in W['s']['sel']] == [z['role'] for z in B['sel']]))
    F = gfirst_frames(W); tau = late_tau(W, F)
    first = {z['role']: F['pos'][F['rops'][z['role']][0]] for z in W['s']['sel'] if F['rops'][z['role']]}
    log('cut %d, gauged %d, clamped deadlines %d' % (F['cut'], len(tau), sum(tau[b] < first.get(b, 10**9) for b in tau)))
    for nds in (False, True):
        C, Wv, info = recount(W, F, [], tau, nds)
        log('no pairs nds=%s: bad %s  D %s  W %s  ==pickle %s' % (nds, info['bad'], D_of(C, Wv, W['m']), Wv,
            (dict(C) == {r: Q(n) for r, n in B['C'].items()} and Wv == B['W']) if nds else ''))
        a = evaluate(W, C, Wv); log('   a_c %s = %.10e' % (a, float(a)))
    res = {}
    for nds in (False, True):
        for chain_ok in (False, True):
            ed = candidates(W, F, tau, a0, nds, chain_ok); pairs = match(ed)
            C, Wv, info = recount(W, F, pairs, tau, nds)
            a = evaluate(W, C, Wv)
            nch = len({x for x, _ in pairs} & {y for _, y in pairs})
            log('pairs nds=%s chain=%s: edges %d pairs %d chained %d bad %s D %s  a_c %s = %.10e' % (nds, chain_ok, len(ed),
                len(pairs), nch, info['bad'], D_of(C, Wv, W['m']), a, float(a)))
            res[nds, chain_ok] = dict(pairs=pairs, a=a, C=dict(C), W=Wv)
    pickle.dump(dict(p=p, tau=tau, res=res), open(os.path.join(HERE, 'stage1_p%d.pkl' % p), 'wb'))
