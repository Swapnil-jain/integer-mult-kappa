"""Own rebuild of #144's selected paired-cube complex word (p=12, h=24), from its published DATA only:
the pinned PR117 DAG (sha 3c034d0a...), the frozen 6,074 matching arcs (sha 22c90ed0...), and the construction as
specified in notes/paired-cube-construction.tex. No upstream code is imported or run; this file re-derives the same
deterministic objects (signed channel graph, coordinate matching frames, carrier annihilators, chronological gauges)
and is checked by EQUALITY against certificates/paired-cube-complex-input.json.

build(p) -> graph dict; compile_frames(g, arcs) -> (profile, witness); select_gauges(g, prof, wit, cost) -> (profile, sel)
cost: callable(role_rank r, gauge_dim d, target list, limit dims) -> delta; default reproduces #144's trial filter.
"""
import gzip, json, math, os, pickle, sys
from itertools import combinations, product
from collections import Counter, defaultdict
from functools import lru_cache
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------- F2 linear algebra on int bitmasks ----------------
def basis(rows):
    """reduced row echelon basis, as a tuple sorted by pivot descending (canonical for the span)"""
    b = {}
    for x in rows:
        for p_, y in sorted(b.items(), reverse=True):
            if x >> p_ & 1: x ^= y
        if x:
            p_ = x.bit_length() - 1
            for k, y in list(b.items()):
                if y >> p_ & 1: b[k] = y ^ x
            b[p_] = x
    return tuple(b[p_] for p_ in sorted(b, reverse=True))

def perp(rows, h):
    rows = basis(rows); piv = {r.bit_length() - 1: r for r in rows}; out = []
    for j in range(h):
        if j in piv: continue
        x = 1 << j
        for p_, r in piv.items():
            if r >> j & 1: x |= 1 << p_
        out.append(x)
    return basis(out)

@lru_cache(maxsize=1 << 20)
def contained(A, B):
    for x in A:
        for y in B: x = min(x, x ^ y)
        if x: return False
    return True

# ---------------- the signed channel graph ----------------
class G:
    def __init__(s, p):
        s.p = p; s.h = 2 * p; s.cubes = list(combinations(range(p), 3))
        s.labels = [tuple(2 * i + b for i, b in zip(I, bits)) for I in s.cubes for bits in product(range(2), repeat=3)]
        s.a = [None] * len(s.labels); s.signs = [1] * len(s.labels); s.sup = [1 << i for i in range(len(s.labels))]
        s.intern = {}; s.src = {(I, bits): 8 * j + k for j, I in enumerate(s.cubes) for k, bits in enumerate(product(range(2), repeat=3))}
        s.F, s.A, s.Gc = {}, {}, {}; s.root = []; s.centers = []
    def add(s, a, b, sign=1):
        if a is None: assert sign == 1; return b
        if b is None: return a
        if sign == 1 and a > b: a, b = b, a
        key = (a, b, sign)
        if key not in s.intern:
            assert not s.sup[a] & s.sup[b]
            s.intern[key] = len(s.a); s.a.append([a, b]); s.signs.append(sign); s.sup.append(s.sup[a] | s.sup[b])
        return s.intern[key]
    def tsum(s, xs):
        xs = list(xs)
        if not xs: return None
        while len(xs) > 1:
            xs = [s.add(xs[i], xs[i + 1]) if i + 1 < len(xs) else xs[i] for i in range(0, len(xs), 2)]
        return xs[0]
    def module(s, data, inputs):
        args = data['args']; n = data.get('input_count', len(inputs))
        one = args[0] == [0, 0] and len(args) > n and args[n] == [0, 0]
        img = ([None] if one else []) + list(inputs)
        for a, b in args[(n + 1 if one else n):]: img.append(s.add(img[a], img[b]))
        return [img[x] for x in data['roots']]
    def local(s):
        for I in s.cubes:
            E = {}
            for ii, jj in combinations(range(3), 2):
                kk = 3 - ii - jj
                for a, b in product(range(2), repeat=2):
                    b0 = [0] * 3; b0[ii], b0[jj] = a, b; b1 = b0[:]; b1[kk] = 1
                    E[ii, jj, a, b] = s.add(s.src[I, tuple(b0)], s.src[I, tuple(b1)])
                s.Gc[I, I[ii], I[jj], 0] = s.add(E[ii, jj, 0, 0], E[ii, jj, 1, 1], -1)
                s.Gc[I, I[ii], I[jj], 1] = s.add(E[ii, jj, 0, 1], E[ii, jj, 1, 0], -1)
            for ii in range(3):
                jj = next(j for j in range(3) if j != ii)
                for a in range(2):
                    ids = [E[(ii, jj, a, b) if ii < jj else (jj, ii, b, a)] for b in range(2)]
                    s.A[I, I[ii], a] = s.add(*ids)
            s.F[I] = s.add(s.A[I, I[0], 0], s.A[I, I[0], 1])

def local_cfg(s, cfg, bug=None):
    """Configured cube-local channels (#168 v3's local_L1 semantics, own code from its docstring and data): the same
    13 outputs per cube, by a different circuit. A[i,a] = sum of the 4 ports with bit i = a, built as two edges along
    direction d ('e<d>') or two face diagonals ('fd'); G[j,k,mode] = [bits (j,k) = (0,mode)] - [bits (j,k) = (1,1-mode)]
    summed over the third bit r, built from two edges along r ('e'), two long-diagonal differences ('l') or two
    fixed-r differences ('s'); F = A[f,0] + A[f,1]. Signed pairs take the smaller port as minuend so equal
    differences are shared. bug (negative controls only): 'A' builds A[i,a] from bit 1-a, 'G' flips one G mode,
    'F' subtracts A[f,1] instead of adding it."""
    A_cfg = {tuple(int(c) for c in k.split(',')): v for k, v in cfg['A'].items()}
    G_cfg = {tuple(int(c) for c in k.split(',')): v for k, v in cfg['G'].items()}
    def signed(xa, xb):
        return (s.add(xa, xb, -1), 1) if xa < xb else (s.add(xb, xa, -1), -1)
    def combine(t1, t2):
        (n1, s1), (n2, s2) = t1, t2
        if s1 == 1: return s.add(n1, n2, s2)
        assert s2 == 1, 'negated local output'
        return s.add(n2, n1, -1)
    for I in s.cubes:
        def at(vals):
            bits = [0] * 3
            for q, b in vals.items(): bits[q] = b
            return s.src[I, tuple(bits)]
        def edge(d, fixed): return s.add(at({**fixed, d: 0}), at({**fixed, d: 1}))
        for j, k in combinations(range(3), 2):
            r = 3 - j - k
            for mode in range(2):
                u, v = 0, (1 - mode if bug == 'G' and (j, k) == (0, 1) else mode)
                kind = G_cfg[j, k, mode]
                if kind == 'e':
                    node = s.add(edge(r, {j: u, k: v}), edge(r, {j: 1 - u, k: 1 - v}), -1)
                else:
                    assert kind in ('l', 's')
                    far = (1, 0) if kind == 'l' else (0, 1)
                    node = combine(signed(at({j: u, k: v, r: 0}), at({j: 1 - u, k: 1 - v, r: far[0]})),
                                   signed(at({j: u, k: v, r: 1}), at({j: 1 - u, k: 1 - v, r: far[1]})))
                s.Gc[I, I[j], I[k], mode] = node
        for i in range(3):
            j, k = [q for q in range(3) if q != i]
            for a in range(2):
                kind = A_cfg[i, a]; aa = 1 - a if bug == 'A' and i == 0 else a
                if kind == 'fd':
                    n1 = s.add(at({i: aa, j: 0, k: 0}), at({i: aa, j: 1, k: 1}))
                    n2 = s.add(at({i: aa, j: 0, k: 1}), at({i: aa, j: 1, k: 0}))
                else:
                    d = int(kind[1]); assert kind[0] == 'e' and d != i
                    o = j if d == k else k
                    n1, n2 = edge(d, {i: aa, o: 0}), edge(d, {i: aa, o: 1})
                s.A[I, I[i], a] = s.add(n1, n2)
        f = cfg['F']
        s.F[I] = s.add(s.A[I, I[f], 0], s.A[I, I[f], 1], -1 if bug == 'F' else 1)

FUSE = None      # #168 v3's 'f8:<8 digits>' fused reads: each port folds its face2, edge02, edge12 singles
LOCAL = None     # a local_L1-style config dict, or None for #144's local channels
LOCAL_BUG = None
LOCALFN = None   # localmotif hook: callable(s, recipe, bug) -> {mode: +-1}

def fuse_roots(s, roots, code):
    """own code of #168 v3's merge_outputs 'f8:' variant: port t folds its three single-target channels (face2,
    edge02, edge12) in the order PERMS[code[t % 8]], one signed node read with the first channel's coefficient.
    Order of the result: the other side roots, then the fused ones by target, then the centres (as upstream)."""
    from fractions import Fraction
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    names = ('face2', 'edge02', 'edge12'); single = {}; keep = []
    side = [r for r in roots if r['kind'] == 'side']; centre = [r for r in roots if r['kind'] != 'side']
    for r in side:
        if r['channel'] in names:
            assert len(r['targets']) == 1; single[r['targets'][0], r['channel']] = (r['node'], Fraction(r['coefficients'][0]))
        else: keep.append(r)
    new = []
    for t in range(len(s.labels)):
        parts = [single[t, ch] for ch in names]; order = perms[int(code[t % 8])]
        node, c0 = parts[order[0]]
        for k in order[1:]:
            nd, c = parts[k]; sg = c / c0; assert sg in (1, -1)
            node = s.add(node, nd, int(sg))
        new.append(dict(node=node, targets=[t], coefficients=[str(c0)], kind='side', channel='fused'))
    return keep + new + centre

def _load_dag():
    return json.loads(gzip.decompress(open(os.path.join(HERE, 'up', 'pr117_dag.json.gz'), 'rb').read()))

def mod_triples(p, w):
    labels = list(combinations(range(p), 3)); index = {t: i + 1 for i, t in enumerate(labels)}
    args = [(0, 0)] + [(0, 0)] * len(labels); sup = [0] + [1 << i for i in range(len(labels))]
    intern = {x: i for i, x in enumerate(sup) if x}
    old = list(combinations(range(w['h']), 3)); img = [0] + [index.get(t, 0) for t in old]
    for aa, bb in zip(w['args'][::2], w['args'][1::2]):
        a, b = img[aa], img[bb]
        if not a or not b: img.append(a or b); continue
        x = sup[a] | sup[b]
        if x not in intern: intern[x] = len(args); args.append((a, b)); sup.append(x)
        img.append(intern[x])
    roots = [img[x] for t, x in zip(old, w['D']) if t[-1] < p]
    act = set(range(len(labels) + 1)); st = roots[:]
    while st:
        x = st.pop()
        if x in act: continue
        act.add(x); st.extend(args[x])
    order = sorted(act); rn = {x: i for i, x in enumerate(order)}
    return dict(input_count=len(labels), args=[[rn[a], rn[b]] for x in order for a, b in [args[x]]], roots=[rn[x] for x in roots])

def mod_pairs(n, w):
    oldt = list(combinations(range(24), 3)); pairs = list(combinations(range(n), 2)); pid = {q: i for i, q in enumerate(pairs)}
    args = [None] * len(pairs); sup = [1 << i for i in range(len(pairs))]; by = {x: i for i, x in enumerate(sup)}; img = [None]
    for t in oldt:
        q = tuple(x for x in t if x != 23); img.append(pid[q] if 23 in t and len(q) == 2 and q[-1] < n else None)
    for aa, bb in zip(w['args'][::2], w['args'][1::2]):
        a, b = img[aa], img[bb]
        if a is None or b is None: img.append(a if b is None else b); continue
        x = sup[a] | sup[b]
        if x not in by: by[x] = len(args); args.append([a, b]); sup.append(x)
        img.append(by[x])
    oldroots = {t: x for t, x in zip(oldt, w['D'])}
    roots = [img[oldroots[(i, j, 22)]] for i, j in pairs]
    act = set(range(len(pairs))); st = list(roots)
    while st:
        x = st.pop()
        if x in act: continue
        act.add(x)
        if args[x] is not None: st.extend(args[x])
    ids = sorted(act); rn = {x: i for i, x in enumerate(ids)}
    return dict(input_count=len(pairs), args=[None if args[x] is None else [rn[y] for y in args[x]] for x in ids], roots=[rn[x] for x in roots])

def mod_allbut(n):
    args = [None] * n; sup = [1 << i for i in range(n)]; by = {x: i for i, x in enumerate(sup)}
    def add(a, b):
        if a is None: return b
        if b is None: return a
        x = sup[a] | sup[b]
        if x not in by: by[x] = len(args); args.append([a, b]); sup.append(x)
        return by[x]
    def tree(lo, hi):
        if hi - lo == 1: return (lo, None, None)
        mid = (lo + hi) // 2; L = tree(lo, mid); R = tree(mid, hi); return (add(L[0], R[0]), L, R)
    T = tree(0, n); roots = [None] * n
    def walk(t, out):
        x, L, R = t
        if L is None: roots[x] = out; return
        walk(L, add(out, R[0])); walk(R, add(out, L[0]))
    walk(T, None)
    act = set(range(n)); st = list(roots)
    while st:
        x = st.pop()
        if x in act: continue
        act.add(x)
        if args[x] is not None: st.extend(args[x])
    ids = sorted(act); rn = {x: i for i, x in enumerate(ids)}
    return dict(input_count=n, args=[None if args[x] is None else [rn[y] for y in args[x]] for x in ids], roots=[rn[x] for x in roots])

PERM = None
RELAXED = False
ORDER_RANK = os.environ.get('ORDER_RANK') == '1'
MERGE = None
def build(p=12):
    w = _load_dag(); s = G(p)
    negm = {0: 1, 1: 1}
    if LOCAL is None: s.local()
    elif LOCALFN is not None: negm = LOCALFN(s, LOCAL, LOCAL_BUG)   # localmotif: generalised recipes
    else: local_cfg(s, LOCAL, LOCAL_BUG)
    D = dict(zip(s.cubes, s.module(mod_triples(p, w), [s.F[I] for I in s.cubes])))
    Pm, Qm = {}, {}; pm = mod_pairs(p - 1, w); am = mod_allbut(p - 2)
    for i in range(p):
        others = [a for a in range(p) if a != i]; prs = list(combinations(others, 2))
        for bit in range(2):
            for K, node in zip(prs, s.module(pm, [s.A[tuple(sorted((i, *K))), i, bit] for K in prs])):
                Pm[tuple(sorted((i, *K))), i, bit] = node
    for i, j in combinations(range(p), 2):
        others = [a for a in range(p) if a not in (i, j)]
        for mode in range(2):
            for k, node in zip(others, s.module(am, [s.Gc[tuple(sorted((i, j, k))), i, j, mode] for k in others])):
                Qm[tuple(sorted((i, j, k))), i, j, mode] = node
    def emit(I, bl, node, num=1, ch=''):
        s.root.append(dict(node=node, targets=[s.src[I, b] for b in bl], coefficients=['%d/2' % num] * len(bl), kind='side', channel=ch))
    allb = list(product(range(2), repeat=3))
    for I in s.cubes:
        # PERM hook (smallh-search): read order of the cube's three positions; default (0, 1, 2) is #144's schedule
        x0, x1, x2 = PERM(I) if PERM else (0, 1, 2)
        def grp(fix): return [b for b in allb if all(b[k] == v for k, v in fix)]
        def edge(xa, xb, bits):
            lo, hi = (xa, xb) if xa < xb else (xb, xa)
            # a G mode whose local output is negated (uniformly, see lmlocal) is read with the opposite sign
            return Qm[I, I[lo], I[hi], bits[lo] ^ bits[hi]], (2 * bits[lo] - 1) * negm[bits[lo] ^ bits[hi]]
        emit(I, allb, D[I], ch='disjoint')
        for a in range(2): emit(I, grp([(x0, a)]), Pm[I, I[x0], 1 - a], ch='face0')
        for a, b in product(range(2), repeat=2):
            bl = grp([(x0, a), (x1, b)]); emit(I, bl, Pm[I, I[x1], 1 - b], ch='face1')
        for a, b in product(range(2), repeat=2):
            bl = grp([(x0, a), (x1, b)]); node, sg = edge(x0, x1, bl[0]); emit(I, bl, node, sg, 'edge01')
        if MERGE is None:
            for bits in allb: emit(I, [bits], Pm[I, I[x2], 1 - bits[x2]], ch='face2')
            for xa in (x0, x1):
                for bits in allb:
                    node, sg = edge(xa, x2, bits); emit(I, [bits], node, sg, 'edge%d2' % xa)
        else:
            # MERGE hook (smallh-search): face2 + edge(xm, x2) summed once per (bit xm, bit x2) pair and read by
            # its two targets; the other singleton edge channel is read as before
            xm = (x0, x1)[MERGE]; xo = (x0, x1)[1 - MERGE]
            for a, b in product(range(2), repeat=2):
                bl = grp([(xm, a), (x2, b)]); en, sg = edge(xm, x2, bl[0])
                node = s.add(Pm[I, I[x2], 1 - b], en, sg)
                for bits in bl: emit(I, [bits], node, 1, 'merge')
            for bits in allb:
                node, sg = edge(xo, x2, bits); emit(I, [bits], node, sg, 'edge%d2' % xo)
    for i in range(p):
        for e in range(2): s.centers.append(s.tsum(s.A[I, i, e] for I in s.cubes if i in I))
    if FUSE is not None:
        # fused nodes are created AFTER the centre nodes, as upstream's merge_outputs runs after finish(); this keeps
        # node ids identical to #168's graph (checked by its graph sha), which frozen arcs and frames refer to
        s.root = fuse_roots(s, s.root, FUSE)
    for c, node in enumerate(s.centers):
        s.root.append(dict(node=node, targets=list(range(len(s.labels))),
                           coefficients=['1/3' if c in t else '-1/6' for t in s.labels], kind='center', coordinate=c))
    return dict(p=p, h=s.h, v=len(s.labels), labels=s.labels, inputs=[sum(1 << x for x in t) for t in s.labels],
                args=s.a, signs=s.signs, roots=s.root, centers=s.centers, matching_frames='coordinate')

# ---------------- frames, matching, rank ledger ----------------
def compile_frames(g, frozen):
    if g['args'][len(g['inputs'])] is not None:
        g = dict(g); g['args'] = [None] + [None if a is None else [x + 1 for x in a] for a in g['args']]
        g['roots'] = [dict(r, node=r['node'] + 1) for r in g['roots']]
    h = g['h']; inputs = g['inputs']; v = len(inputs); args = g['args']; n = len(args); roots = g['roots']; q = len(roots)
    spans = [()] * n
    for i, u in enumerate(inputs): spans[i + 1] = (u,)
    for x in range(v + 1, n):
        a, b = args[x]; spans[x] = basis(spans[a] + spans[b])
    rframe, rann = [], []; Y = [()] * v; TH = Counter(); ell = 0
    for r in roots:
        x = r['node']
        if r['kind'] == 'center':
            U = spans[x]; A = perp(U, h); ell += len(U)
        else:
            A = basis(inputs[t] for t in r['targets']); U = perp(A, h)
            for t in r['targets']:
                assert contained(Y[t], U); TH[len(U) - len(Y[t])] += 1; Y[t] = U
        assert contained(spans[x], U); rframe.append(U); rann.append(A)
    for t in range(v):
        assert contained(Y[t], perp((inputs[t],), h)); TH[h - 1 - len(Y[t])] += 1
    active = set(range(1, v + 1)); todo = [r['node'] for r in roots]
    while todo:
        x = todo.pop()
        if x in active: continue
        active.add(x)
        if args[x]: todo.extend(args[x])
    initial = [()] * n; pre = [None] * n; succ0 = [[] for _ in args]; cons = [[] for _ in args]
    for x in sorted(active):
        if args[x]:
            for y in args[x]: succ0[y].append(x)
    for j, r in enumerate(roots): cons[r['node']].extend(rann[j])
    for x in sorted(active, reverse=True):
        pre[x] = basis(cons[x] + [z for y in succ0[x] for z in pre[y]])
        cover = 0
        for u in spans[x]: cover |= u
        initial[x] = perp(pre[x] + tuple(1 << j for j in range(h) if not cover >> j & 1), h)
        assert contained(spans[x], initial[x])
    order = sorted(active, key=lambda x: (len(initial[x]), x)); pos = {x: i for i, x in enumerate(order)}
    uses = [[] for _ in args]; uval = []; utgt = []; ufr = []; ucode = []
    for x in order:
        if args[x]:
            for j, y in enumerate(args[x]):
                uses[y].append(len(uval)); uval.append(y); utgt.append(x); ufr.append(initial[x]); ucode.append(2 * x + j)
    for j, r in enumerate(roots):
        x = r['node']; uses[x].append(len(uval)); uval.append(x); utgt.append(n + j); ufr.append(rframe[j]); ucode.append((1 << 31) | j)
    donors = [x for x in order if args[x]]; bycode = {c: u for u, c in enumerate(ucode)}
    if callable(frozen):   # own matching: pass the eligibility data, get back [(donor, use code)]
        frozen = frozen(dict(order=order, pos=pos, args=args, uval=uval, utgt=utgt, ufr=ufr, ucode=ucode,
                             initial=initial, n=n, uses=uses, h=h, rframe=rframe, rann=rann, roots=roots, v=v, spans=spans))
    arcs = {x: bycode[c] for x, c in frozen}; assert len(arcs) == len(frozen) == len(set(arcs.values()))
    for x, u in arcs.items():
        assert args[x] and uval[u] in args[x]; t = utgt[u]
        if not RELAXED: assert t >= n or pos[t] > pos[x]; assert contained(initial[x], ufr[u])
    if RELAXED:
        # #162-style acceptance (cube-birth): an arc needs only (1) DAG + arc edges acyclic and (2) every value span
        # inside its full backward-intersection frame (checked after the frames below). The coordinate-padded
        # pre-test above is sufficient, not necessary. Re-order the ops as a topological order of DAG + arc edges,
        # keyed (dim initial frame, node) as before, so every arc runs forward.
        import heapq
        outs = defaultdict(list); indeg = Counter()
        for x in active:
            if args[x]:
                for y in set(args[x]): outs[y].append(x); indeg[x] += 1
        for x, u in arcs.items():
            t = utgt[u]
            if t < n: outs[x].append(t); indeg[t] += 1
        # ORDER_RANK: compile_closure's plain heap on node ids (its generalized topological index)
        hq = [((0 if ORDER_RANK else len(initial[x])), x) for x in active if not indeg[x]]; heapq.heapify(hq); order = []
        while hq:
            _, x = heapq.heappop(hq); order.append(x)
            for t in outs[x]:
                indeg[t] -= 1
                if not indeg[t]: heapq.heappush(hq, ((0 if ORDER_RANK else len(initial[t])), t))
        assert len(order) == len(active), 'arc graph has a cycle'
        pos = {x: i for i, x in enumerate(order)}
        uses = [[] for _ in args]
        for x in order:
            if args[x]:
                for j, y in enumerate(args[x]): uses[y].append(bycode[2 * x + j])
        for j, r in enumerate(roots): uses[r['node']].append(bycode[(1 << 31) | j])
    succ = [[] for _ in args]; direct = [[] for _ in args]
    for x in order:
        if args[x]:
            for y in args[x]: succ[y].append(x)
    for j, r in enumerate(roots): direct[r['node']].extend(rann[j])
    for x, u in arcs.items():
        t = utgt[u]
        if t < n: succ[x].append(t)
        else: direct[x].extend(rann[t - n])
    ann = [None] * n; rank = [0] * n
    for x in reversed(order):
        ann[x] = basis(direct[x] + [z for y in succ[x] for z in ann[y]]); rank[x] = h - len(ann[x])
    if ORDER_RANK:
        # #162/#168 compile_closure's op order: (frame rank, generalized topological index). It is a linear extension
        # (a successor's annihilator is inside its predecessor's, so its rank is >= ), and it changes which use of
        # each value is served first, hence the role chains. The annihilators do not depend on the order.
        idx = {x: i for i, x in enumerate(order)}; order = sorted(order, key=lambda x: (rank[x], idx[x]))
        pos = {x: i for i, x in enumerate(order)}
        for x in order:
            for y in succ[x]: assert pos[y] > pos[x], 'rank order is not a linear extension'
    if RELAXED:
        for x in order:
            assert all(bin(a & b).count('1') % 2 == 0 for a in spans[x] for b in ann[x]), 'value span outside its frame'
    H = Counter()
    for x in order:
        r = rank[x]; H[r] += len(uses[x]) - 1
        if args[x]:
            H[h - r] += 1
            for y in args[x]: H[r - rank[y]] += 1
        else: H[1] += 1; H[r - 1] += 1
    for j, root in enumerate(roots):
        x = root['node']; r = rank[x]; rt = len(rframe[j])
        if root['kind'] == 'center': H[r] += 1; H[h - r] += 1
        else: H[rt - r] += 1; H[h - rt] += 1
    for x, u in arcs.items():
        val = uval[u]; t = utgt[u]; rv = rank[val]; rd = rank[x]; rt = rank[t] if t < n else len(rframe[t - n])
        H[h - rd] -= 1; H[rv] -= 1; H[rt - rv] -= 1; H[rt - rd] += 1
    assert min(H.values()) >= 0
    R = len(donors) + q - len(arcs); assert sum(r * c for r, c in H.items()) == h * R + ell
    prof = dict(h=h, v=v, R=R, q=q, c=len(donors), matched=len(arcs), loss=ell, H=H, TH=TH,
                source={2: v, h - 4: v, 1: v})
    wit = dict(arcs=arcs, ann=ann, order=order, rank=rank, uses=uses, uval=uval, utgt=utgt, ucode=ucode, rframe=rframe,
               rann=rann, g=g, n=n)
    return prof, wit

def excess_fn(a, m):
    return lambda t: t * math.expm1(a * math.log(m / t)) if t else 0.

def select_gauges(prof, wit, tail=None, trial=0.00065, accept_all=False, variants=None):
    """#144's chronological partial-gauge selection (ownership of every use slot rebuilt from the frozen arcs).
    tail(d) -> list of (rank, weight) children replacing #144's single 3d tail child per gauged role (NDS hook);
    returns the per-vertex child multiset pieces and the selected gauges (role, annihilator A, dim d, targets)."""
    g = wit['g']; h = prof['h']; v = prof['v']; args = g['args']; roots = g['roots']; order = wit['order']
    ann = wit['ann']; arcs = wit['arcs']; ucode = wit['ucode']; m = 3 * h
    arcs_c = {x: ucode[u] for x, u in arcs.items()}; incoming = set(arcs_c.values())
    uses = [[] for _ in args]
    for x in order:
        if args[x]:
            for j, y in enumerate(args[x]): uses[y].append(2 * x + j)
    for j, r in enumerate(roots): uses[r['node']].append((1 << 31) | j)
    assign = {}; ops = []; R = 0; sources = {}
    for x in order:
        if args[x]:
            aa, bb = args[x]; dest, ctl = assign[2 * x], assign[2 * x + 1]
            if x in arcs_c:
                u = arcs_c[x]; val = roots[u & 0x7fffffff]['node'] if u >> 31 else args[u // 2][u & 1]
                if val == aa: dest, ctl = ctl, dest
                else: assert val == bb
                assign[u] = ctl
            ops.append((dest, ctl, x))
        else: dest = R; R += 1; sources[x] = dest
        free = [u for u in uses[x] if u not in incoming]; assert free
        for j, u in enumerate(free):
            if j == 0: assign[u] = dest
            else: t = R; R += 1; assign[u] = t; ops.append((t, dest, x))
    assert R == prof['R']
    rootroles = [assign[(1 << 31) | j] for j in range(len(roots))]
    prev = [-1] * R; pred = []; first = [None] * R
    for i, (a, b, x) in enumerate(ops):
        pred.append((prev[a], prev[b])); prev[a] = prev[b] = i
        if first[a] is None: first[a] = x
        if first[b] is None: first[b] = x
    st = [prev[s] for r, s in zip(roots, rootroles) if r['kind'] == 'center' and prev[s] >= 0]; phase = set()
    while st:
        i = st.pop()
        if i in phase: continue
        phase.add(i); st.extend(j for j in pred[i] if j >= 0)
    touched = set(sources.values())
    for i in phase: touched.update(ops[i][:2])
    co = [0] * R
    for r, s in zip(roots, rootroles): co[s] |= sum(1 << t for t in r['targets'])
    for a, b, x in reversed(ops): co[b] |= co[a]
    inputs = g['inputs']; limit = [None] * v
    for r in roots:
        if r['kind'] == 'center': continue
        A = basis(inputs[t] for t in r['targets'])
        for t in r['targets']:
            if limit[t] is None: limit[t] = A
    for t in range(v):
        if limit[t] is None: limit[t] = (inputs[t],)
    H = Counter(prof['H']); gauges = Counter(); sel = []
    cands = sorted((s for s in range(R) if s not in touched and first[s] is not None), key=lambda s: (len(ann[first[s]]), bin(co[s]).count("1"), s))
    ex = excess_fn(trial, m)
    def tailcost(d, A):
        if tail is None: return ex(3 * d)
        return sum(float(w) * ex(r) for r, w in tail(d, A))
    rej = 0
    for s in cands:
        A = tuple(ann[first[s]]); tg = []; bits = co[s]
        while bits:
            low = bits & -bits; t = low.bit_length() - 1; bits ^= low; tg.append(t)
            A = basis(A + tuple(limit[t]))
            if len(A) == h: break
        if len(A) == h: continue
        r = h - len(ann[first[s]])
        best = None
        for A1 in ([A] if variants is None else variants(A)):
            d = h - len(A1)
            delta = 3 * (ex(r - d) - ex(r)) + tailcost(d, A1)
            delta += 3 * math.fsum(ex(d) + ex(h - len(limit[t]) - d) - ex(h - len(limit[t])) for t in tg)
            if best is None or delta < best[0] - 1e-12: best = (delta, A1, d)
        delta, A, d = best
        if delta >= -1e-12 and not accept_all: rej += 1; continue
        H[r] -= 1; H[r - d] += 1; gauges[d] += 1
        for t in tg: limit[t] = A
        sel.append(dict(role=s, A=A, d=d, targets=tg, r=r))
    Y = Counter(); cur = [()] * v
    for z in reversed(sel):
        U = perp(z['A'], h)
        for t in z['targets']:
            assert contained(cur[t], U); Y[len(U) - len(cur[t])] += 1; cur[t] = U
    for r in roots:
        if r['kind'] == 'center': continue
        U = perp(basis(inputs[t] for t in r['targets']), h)
        for t in r['targets']:
            assert contained(cur[t], U); Y[len(U) - len(cur[t])] += 1; cur[t] = U
    for t in range(v): Y[h - 1 - len(cur[t])] += 1
    out = dict(prof, H=H, TH=Y, gauges=gauges, nsel=len(sel), rejected=rej, phase1=len(phase), ops=len(ops),
               untouched=R - len(touched))
    return out, dict(sel=sel, ops=ops, phase=phase, sources=sources, rootroles=rootroles, first=first, touched=touched)

def _wsave(d, A, tail, ex):
    # placeholder for a W-aware gauge filter (set by callers through tail); default none
    return 0.0

def children(pr, extra=None):
    """per-vertex child multiset in #144's normalisation: 3(H+source+target) + 2v[2] + gauge tails"""
    v = pr['v']; C = Counter()
    for hist in (pr['H'], pr['source'], pr['TH']):
        for r, c in hist.items():
            if r: C[r] += 3 * c
    C[2] += 2 * v
    for d, c in pr.get('gauges', {}).items(): C[3 * d] += c
    return +C

def published():
    return json.load(open(os.path.join(HERE, '..', 'cover-bit', 'pr144', 'certificates__paired-cube-complex-input.json')))

if __name__ == '__main__':
    import time; t0 = time.time()
    g = build(12)
    import hashlib
    bind = {k: g[k] for k in ('inputs', 'labels', 'args', 'signs', 'roots', 'centers')}
    gh = hashlib.sha256(json.dumps(bind, separators=(',', ':')).encode()).hexdigest()
    print('graph sha', gh, 'pinned 7a7a58df...', gh.startswith('7a7a58df'), '%.1fs' % (time.time() - t0), flush=True)
    frozen = json.load(open(os.path.join(HERE, 'up', 'references__paired-cube__selected-module__matching-arcs.json')))
    prof, wit = compile_frames(g, frozen); print('compiled', prof['R'], prof['c'], prof['q'], '%.1fs' % (time.time() - t0), flush=True)
    pr, sel = select_gauges(prof, wit); P = published()
    C = children(pr)
    print('R', pr['R'], 'c', pr['c'], 'q', pr['q'], 'gauges', dict(pr['gauges']), 'untouched', pr['untouched'], 'phase1', pr['phase1'], 'ops', pr['ops'])
    print('internal == pub', [pr['H'][r] for r in range(25)] == P['remaining_internal_histogram'])
    print('target == pub', {str(k): x for k, x in pr['TH'].items()} == P['target_data_histogram'])
    print('children == pub', {str(k): x for k, x in C.items()} == {k: x for k, x in P['child_histogram'].items() if x})
    pickle.dump(dict(prof=prof, pr=pr, sel=sel['sel'], first=sel['first'], ops=sel['ops'], phase=sel['phase'],
                     touched=sel['touched'], ann=wit['ann'], order=wit['order'], rank=wit['rank'], arcs=wit['arcs'],
                     rootroles=sel['rootroles'], sources=sel['sources']), open(os.path.join(HERE, 'word144.pkl'), 'wb'))
    print('done %.1fs' % (time.time() - t0))
