"""Generalised cube-local channel circuits (own code). A recipe fixes how each cube computes its 13 local outputs
  F       = sum of the 8 ports
  A[i,a]  = sum of the 4 ports with bit i = a                               (6 outputs)
  G[j,k,m]= [bits (j,k) = (0,m)] - [bits (j,k) = (1,1-m)], summed over r    (6 outputs)
The VALUES are fixed by #144's decoder identity; only the addition circuit changes. Legality of a recipe is therefore
(a) every addition combines disjoint supports with sign +-1 (asserted by cube_word.G.add), (b) each output has its
required value, up to one sign per (G mode) that is uniform over the three position pairs (a negated mode is read
with the opposite coefficient; the all-but-one module only sums), (c) checked independently by the identity test.

Recipe dict (JSON-friendly):
  A: {"i,a": kind}, kind in e<d> (d != i: two edges along d), fd (two face diagonals),
                    c<P> (a chain ((x+y)+z)+w over the face's ports in the order P, a permutation code 0..23)
  G: {"j,k,m": kind}, kind in e (two edges along r), l (long-diagonal differences), s (fixed-r differences),
                    c<P> (chain over the 4 signed ports in order P)
  F: f in 0,1,2 (A[f,0] + A[f,1]) or "p<f>" (even + odd parity halves from axis-f face diagonals)
Orientation: signed pairs take the smaller port as minuend (shares equal differences), as in local_L1."""
from itertools import combinations, permutations

PERM4 = list(permutations(range(4)))


def gen_local(s, cfg, bug=None):
    A_cfg = {tuple(int(c) for c in k.split(',')): v for k, v in cfg['A'].items()}
    G_cfg = {tuple(int(c) for c in k.split(',')): v for k, v in cfg['G'].items()}
    negmode = {}

    def signed(xa, xb):
        return (s.add(xa, xb, -1), 1) if xa < xb else (s.add(xb, xa, -1), -1)

    def combine(t1, t2):
        """returns (node, sign) with node value = sign * (v1 + v2) where t = (node, sign) means value sign*node"""
        (n1, s1), (n2, s2) = t1, t2
        if s1 == 1: return s.add(n1, n2, s2), 1
        if s2 == 1: return s.add(n2, n1, -1), 1
        return s.add(n1, n2, 1), -1

    def chain(terms):
        """terms: list of (port node, coefficient +-1) folded left; returns (node, sign)"""
        (n, sg) = terms[0]
        cur = (n, sg)
        for (m_, c) in terms[1:]:
            cur = combine(cur, (m_, c))
        return cur

    for I in s.cubes:
        def at(vals):
            bits = [0] * 3
            for q, b in vals.items(): bits[q] = b
            return s.src[I, tuple(bits)]

        def edge(d, fixed): return s.add(at({**fixed, d: 0}), at({**fixed, d: 1}))
        Gsign = {}
        for j, k in combinations(range(3), 2):
            r = 3 - j - k
            for mode in range(2):
                u, v = 0, (1 - mode if bug == 'G' and (j, k) == (0, 1) else mode)
                kind = G_cfg[j, k, mode]
                if kind == 'e':
                    node, sg = s.add(edge(r, {j: u, k: v}), edge(r, {j: 1 - u, k: 1 - v}), -1), 1
                elif kind in ('l', 's'):
                    far = (1, 0) if kind == 'l' else (0, 1)
                    node, sg = combine(signed(at({j: u, k: v, r: 0}), at({j: 1 - u, k: 1 - v, r: far[0]})),
                                       signed(at({j: u, k: v, r: 1}), at({j: 1 - u, k: 1 - v, r: far[1]})))
                else:
                    assert kind[0] == 'c'
                    P = PERM4[int(kind[1:])]
                    terms = [(at({j: u, k: v, r: 0}), 1), (at({j: u, k: v, r: 1}), 1),
                             (at({j: 1 - u, k: 1 - v, r: 0}), -1), (at({j: 1 - u, k: 1 - v, r: 1}), -1)]
                    terms = [terms[x] for x in P]
                    if terms[0][1] == -1: terms = [(n_, -c) for n_, c in terms]; flip = -1
                    else: flip = 1
                    node, sg = chain(terms); sg *= flip
                s.Gc[I, I[j], I[k], mode] = node; Gsign[j, k, mode] = sg
        for mode in range(2):
            sgs = {Gsign[j, k, mode] for j, k in combinations(range(3), 2)}
            assert len(sgs) == 1, 'G mode %d sign not uniform over position pairs' % mode
            sg = sgs.pop()
            assert negmode.setdefault(mode, sg) == sg
        for i in range(3):
            j, k = [q for q in range(3) if q != i]
            for a in range(2):
                kind = A_cfg[i, a]; aa = 1 - a if bug == 'A' and i == 0 else a
                if kind == 'fd':
                    n1 = s.add(at({i: aa, j: 0, k: 0}), at({i: aa, j: 1, k: 1}))
                    n2 = s.add(at({i: aa, j: 0, k: 1}), at({i: aa, j: 1, k: 0}))
                    s.A[I, I[i], a] = s.add(n1, n2)
                elif kind[0] == 'e':
                    d = int(kind[1]); assert d != i
                    o = j if d == k else k
                    s.A[I, I[i], a] = s.add(edge(d, {i: aa, o: 0}), edge(d, {i: aa, o: 1}))
                else:
                    assert kind[0] == 'c'
                    P = PERM4[int(kind[1:])]
                    ports = [at({i: aa, j: b1, k: b2}) for b1 in range(2) for b2 in range(2)]
                    node = ports[P[0]]
                    for x in P[1:]: node = s.add(node, ports[x])
                    s.A[I, I[i], a] = node
        f = cfg['F']
        if isinstance(f, str) and f[0] == 'p':
            f = int(f[1]); j, k = [q for q in range(3) if q != f]
            d = {}
            for a in range(2):
                d[a, 0] = s.add(at({f: a, j: 0, k: 0}), at({f: a, j: 1, k: 1}))
                d[a, 1] = s.add(at({f: a, j: 0, k: 1}), at({f: a, j: 1, k: 0}))
            ev, od = s.add(d[0, 0], d[1, 1]), s.add(d[0, 1], d[1, 0])
            s.F[I] = s.add(ev, od, -1 if bug == 'F' else 1)
        else:
            s.F[I] = s.add(s.A[I, I[f], 0], s.A[I, I[f], 1], -1 if bug == 'F' else 1)
    return {m: negmode[m] for m in negmode}


# ---- configured-space encoding (the #168 v3 space plus F parity), index <-> recipe
AKEYS = [(i, a) for i in range(3) for a in range(2)]
GKEYS = [(j, k, m) for j, k in combinations(range(3), 2) for m in range(2)]
FOPTS = [0, 1, 2, 'p0', 'p1', 'p2']


def akinds(i): return ['e%d' % d for d in range(3) if d != i] + ['fd']


def recipe(ach, gch, f):
    """ach: 6 ints in 0..2 (index into akinds(i)), gch: 6 ints in 0..2 (e,l,s), f: FOPTS entry"""
    return dict(A={'%d,%d' % key: akinds(key[0])[c] for key, c in zip(AKEYS, ach)},
                G={'%d,%d,%d' % key: 'els'[c] for key, c in zip(GKEYS, gch)}, F=f)


def encode(cfg):
    ach = [akinds(i).index(cfg['A']['%d,%d' % (i, a)]) for i, a in AKEYS]
    gch = ['els'.index(cfg['G']['%d,%d,%d' % key]) for key in GKEYS]
    return ach, gch, cfg['F']


def key(cfg):
    return '|'.join('%s=%s' % (k, cfg['A'][k]) for k in sorted(cfg['A'])) + '#' + \
        '|'.join('%s=%s' % (k, cfg['G'][k]) for k in sorted(cfg['G'])) + '#F=%s' % cfg['F']
