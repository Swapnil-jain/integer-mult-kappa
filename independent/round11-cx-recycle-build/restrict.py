"""Restrictions of the pinned PR117 h=24 disjoint-sum DAG to label subsets (own code, generalising cube_word's
mod_triples / mod_pairs, which use labels 0..p-1 and star label 23 / reader label 22).
tri(labels): triple-exclusion module over the given labels (in the given order).
pair(labels, star, reader): pair-exclusion module over labels, inputs x_{c,d,star}, roots D_{(i,j,reader)}."""
import gzip, json, os
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__))
W = json.loads(gzip.decompress(open(os.path.join(HERE, 'up', 'pr117_dag.json.gz'), 'rb').read()))
OLD = list(combinations(range(24), 3)); ROOT = dict(zip(OLD, W['D']))

def _restrict(inmap, n_in, roots_old):
    """inmap: old triple -> new input index or None; returns module dict (simple format)"""
    args = [None] * n_in; sup = [1 << i for i in range(n_in)]; by = {x: i for i, x in enumerate(sup)}
    img = [None] + [inmap(t) for t in OLD]
    a = W['args']
    for aa, bb in zip(a[::2], a[1::2]):
        x, y = img[aa], img[bb]
        if x is None or y is None: img.append(x if y is None else y); continue
        s = sup[x] | sup[y]
        if s not in by: by[s] = len(args); args.append([x, y]); sup.append(s)
        img.append(by[s])
    roots = [img[ROOT[t]] for t in roots_old]
    act = set(range(n_in)); st = [r for r in roots if r is not None]
    while st:
        x = st.pop()
        if x in act: continue
        act.add(x)
        if args[x] is not None: st.extend(args[x])
    ids = sorted(act); rn = {x: i for i, x in enumerate(ids)}
    return dict(input_count=n_in, args=[None if args[x] is None else [rn[y] for y in args[x]] for x in ids],
                roots=[rn[x] for x in roots])

def tri(labels):
    labels = list(labels); pos = {l: i for i, l in enumerate(labels)}
    new = list(combinations(range(len(labels)), 3)); idx = {t: i for i, t in enumerate(new)}
    def inmap(t):
        if all(l in pos for l in t): return idx[tuple(sorted(pos[l] for l in t))]
        return None
    # root for new triple (i,j,k) = old D of the old labels (sorted)
    roots_old = [tuple(sorted(labels[i] for i in t)) for t in new]
    return _restrict(inmap, len(new), roots_old)

def pair(labels, star=23, reader=22):
    labels = list(labels); assert star not in labels and reader not in labels and star != reader
    pos = {l: i for i, l in enumerate(labels)}
    new = list(combinations(range(len(labels)), 2)); idx = {t: i for i, t in enumerate(new)}
    def inmap(t):
        if star in t:
            q = [l for l in t if l != star]
            if all(l in pos for l in q): return idx[tuple(sorted(pos[l] for l in q))]
        return None
    roots_old = [tuple(sorted((labels[i], labels[j], reader))) for i, j in new]
    return _restrict(inmap, len(new), roots_old)

def nadd(m): return sum(a is not None for a in m['args'])

if __name__ == '__main__':
    import sys, random
    sys.path.insert(0, os.path.join(HERE, 'lib')); import cube_word as cw
    for p in (11, 12, 13):
        print('p', p, 'tri default', nadd(tri(range(p))), 'cw', nadd(cw.mod_triples(p, W)), ' pair default', nadd(pair(range(p - 1))), 'cw', nadd(cw.mod_pairs(p - 1, W)))
    rnd = random.Random(1)
    for p in (12,):
        vals = []
        for _ in range(200):
            L = rnd.sample(range(24), p); vals.append((nadd(tri(L)), L))
        vals.sort(); print('tri p=12 random subsets: min', vals[0][0], 'median', vals[100][0], 'max', vals[-1][0])
        vals = []
        for _ in range(200):
            st, rd = rnd.sample(range(24), 2); L = rnd.sample([x for x in range(24) if x not in (st, rd)], p - 1)
            vals.append((nadd(pair(L, st, rd)), L, st, rd))
        vals.sort(); print('pair n=11 random: min', vals[0][0], 'median', vals[100][0], 'max', vals[-1][0])

def tri_from(path, keep):
    """Triple-exclusion module restricted from any #117-format witness (gzip JSON with h, args, D) to the kept points
    (own code; same zero-restriction as tri(), witness given as data)."""
    w = json.loads(gzip.decompress(open(os.path.join(HERE, path), 'rb').read()))
    keep = sorted(keep); pos = {l: i for i, l in enumerate(keep)}
    old = list(combinations(range(w['h']), 3)); new = list(combinations(range(len(keep)), 3)); idx = {t: i for i, t in enumerate(new)}
    args = [None] * len(new); sup = [1 << i for i in range(len(new))]; by = {x: i for i, x in enumerate(sup)}
    img = [None] + [idx[tuple(pos[l] for l in t)] if all(l in pos for l in t) else None for t in old]
    a = w['args']
    for aa, bb in zip(a[::2], a[1::2]):
        x, y = img[aa], img[bb]
        if x is None or y is None: img.append(x if y is None else y); continue
        s = sup[x] | sup[y]; assert not sup[x] & sup[y]
        if s not in by: by[s] = len(args); args.append([x, y]); sup.append(s)
        img.append(by[s])
    R = dict(zip(old, w['D'])); roots = [img[R[tuple(keep[i] for i in t)]] for t in new]
    for J, r in zip(new, roots):   # contract: root J sums exactly the inputs disjoint from J
        assert r is not None and sup[r] == sum(1 << i for i, I in enumerate(new) if not set(I) & set(J))
    act = set(range(len(new))); st = list(roots)
    while st:
        x = st.pop()
        if x in act: continue
        act.add(x)
        if args[x] is not None: st.extend(args[x])
    ids = sorted(act); rn = {x: i for i, x in enumerate(ids)}
    return dict(input_count=len(new), args=[None if args[x] is None else [rn[y] for y in args[x]] for x in ids],
                roots=[rn[x] for x in roots])

def allbut_prefix(n):
    """Nested-prefix all-but-one module (#168's construction, own code): prefix sums p_{k+1} = p_k + x_k, suffix sums
    s_k = x_k + s_{k+1}; y_0 = s_1, y_{n-1} = p_{n-1}, y_{i+1} = p_i + (x_i + s_{i+2}). 3n - 6 additions; each prefix
    p_i (i >= 2) feeds both p_{i+1} and y_{i+1}, with nested supports, which is what lets more carrier arcs in."""
    args = [None] * n; sup = [1 << i for i in range(n)]; by = {x: i for i, x in enumerate(sup)}
    def add(a, b):
        s = sup[a] | sup[b]; assert not sup[a] & sup[b]
        if s not in by: by[s] = len(args); args.append([a, b]); sup.append(s)
        return by[s]
    P = {1: 0}
    for k in range(1, n - 1): P[k + 1] = add(P[k], k)
    S = {n - 1: n - 1}
    for k in range(n - 2, 0, -1): S[k] = add(k, S[k + 1])
    roots = [None] * n; roots[0] = S[1]; roots[n - 1] = P[n - 1]
    for i in range(n - 2):
        inner = add(i, S[i + 2])
        roots[i + 1] = add(P[i], inner) if i >= 1 else inner
    full = (1 << n) - 1
    assert all(sup[r] == full ^ (1 << i) for i, r in enumerate(roots))
    return dict(input_count=n, args=args, roots=roots)

def allbut_mixed(n, pattern):
    """#173's mixed-gap variant of the nested prefix (own code, from its description): each middle output y_j
    (j = 1..n-2) is either L: p_{j-1} + (x_{j-1} + s_{j+1}) (#168's form; p_0 empty) or R: (p_j + x_{j+1}) + s_{j+2}
    (s_n empty). pattern: string of n-2 letters L/R. All-L equals allbut_prefix(n)."""
    assert len(pattern) == n - 2 and set(pattern) <= set('LR')
    args = [None] * n; sup = [1 << i for i in range(n)]; by = {x: i for i, x in enumerate(sup)}
    def add(a, b):
        if a is None: return b
        if b is None: return a
        s = sup[a] | sup[b]; assert not sup[a] & sup[b]
        if s not in by: by[s] = len(args); args.append([a, b]); sup.append(s)
        return by[s]
    P = {0: None, 1: 0}
    for k in range(1, n - 1): P[k + 1] = add(P[k], k)
    S = {n: None, n - 1: n - 1}
    for k in range(n - 2, 0, -1): S[k] = add(k, S[k + 1])
    roots = [None] * n; roots[0] = S[1]; roots[n - 1] = P[n - 1]
    for j in range(1, n - 1):
        if pattern[j - 1] == 'L': roots[j] = add(P[j - 1], add(j - 1, S[j + 1]))
        else: roots[j] = add(add(P[j], j + 1), S[j + 2])
    full = (1 << n) - 1
    assert all(sup[r] == full ^ (1 << i) for i, r in enumerate(roots))
    act = set(range(n)); st = list(roots)
    while st:
        x = st.pop()
        if x in act: continue
        act.add(x)
        if args[x] is not None: st.extend(args[x])
    ids = sorted(act); rn = {x: i for i, x in enumerate(ids)}
    return dict(input_count=n, args=[None if args[x] is None else [rn[y] for y in args[x]] for x in ids], roots=[rn[x] for x in roots])
