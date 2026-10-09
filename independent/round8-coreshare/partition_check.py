"""Independent check of PR #128's 87-group triple partition (read as data, their code is not run).
Checks: exact cover of all C(24,3) triples (zero-based lexicographic indexing), group sizes, binary Gram = I in every
group, each partial-group complement basis spans exactly the F2-orthogonal complement of the group's span, and the
quadratic-phase identity wt(P_U x) summed over an orthonormal basis == wt(x) mod 4 (exhaustive over F2^6 blocks is not
enough, so we check the identity C_U C_V = C_{U+V} through the general lemma: for orthonormal u_i, sum_i wt(P_i x) -
wt(sum_i P_i x) = 2 sum_{i<j} (x.u_i)(x.u_j)|u_i & u_j|, which is 0 mod 4 iff every pairwise support overlap is even).
Usage: python3 -I partition_check.py <pr-network-mrp24-signed.json>"""
import json, sys, itertools, random

h = 24
d = json.load(open(sys.argv[1]))
T = list(itertools.combinations(range(h), 3))
mask = [sum(1 << p for p in t) for t in T]
G = d['groups']
pc = lambda x: bin(x).count('1')
dot = lambda a, b: pc(a & b) & 1
flat = [i for g in G for i in g]
assert sorted(flat) == list(range(len(T))) and len(T) == 2024, 'cover'
sizes = sorted(len(g) for g in G); from collections import Counter; print('sizes', Counter(sizes))
for g in G:
    for a in g:
        assert dot(mask[a], mask[a]) == 1
        for b in g:
            if a != b: assert dot(mask[a], mask[b]) == 0, 'gram'
            if a != b: assert pc(mask[a] & mask[b]) % 2 == 0      # even overlaps -> phases add exactly
print('cover + Gram I + even overlaps: PASS')


def rank(vs):
    piv = {}
    for v in vs:
        while v:
            p = v & -v
            if p in piv: v ^= piv[p]
            else: piv[p] = v; break
    return len(piv), piv


def inspan(v, piv):
    while v:
        p = v & -v
        if p not in piv: return False
        v ^= piv[p]
    return True


parts = [g for g in G if len(g) < h]
comps = d['partial_complement_basis_supports']
assert len(parts) == len(comps) == 4
for g, cb in zip(parts, comps):
    cm = [sum(1 << p for p in s) for s in cb]
    assert len(cm) == 16 and rank(cm)[0] == 16
    for c_ in cm:
        for a in g: assert dot(c_, mask[a]) == 0, 'complement not orthogonal'
    # orthonormal complement with even pairwise overlaps (so C on it is the product of its rank-one pieces)
    for i, a in enumerate(cm):
        assert dot(a, a) == 1
        for b in cm[i + 1:]: assert dot(a, b) == 0 and pc(a & b) % 2 == 0
    r, piv = rank([mask[a] for a in g] + cm); assert r == h
    nwt = Counter(pc(c_) for c_ in cm)
    # phase sign of the complement directions: weight 3 -> i^3 = -i (C^-1 direction), weight 1 -> +i (C direction)
    print('partial group: size', len(g), 'complement weights', dict(nwt))
# direct exhaustive phase check on random x for every group + complement: sum_i wt((x.u_i) u_i) = wt(sum) mod 4
rnd = random.Random(5)
for g in G:
    us = [mask[a] for a in g]
    for _ in range(200):
        x = rnd.getrandbits(h)
        proj = 0; tot = 0
        for u in us:
            if dot(x, u): proj ^= u; tot += pc(u)
        assert (tot - pc(proj)) % 4 == 0
print('phase additivity on random points: PASS')
