"""H-pair exact-ledger bound, disjoint version: each dead copy pair (two retiring rows holding the same node) hosts at
most one fresh birth; both rows are used once. Max-weight bipartite matching (pairs x births), exact recount + root.
Excludes rows whose dirty content depends on a late-read (gauged) value. Search/bound evidence."""
import sys
from wl import *
from collections import Counter, defaultdict
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
w = Word(sys.argv[1]); h, m = w.h, w.m
src = set(w.sources.values()); fresh = [sr for sr in range(w.R) if sr not in src and sr not in w.gsig and sr not in w.sinks and w.rops[sr]]
dep = [set() for _ in range(w.R)]
for b in w.gsig: dep[b].add(b)
for i in w.sched:
    d, c, _ = w.ops[i]; dep[d] |= dep[c]
bynode = defaultdict(list)
for slot in w.slots():
    last = slot[-1]
    if w.rootkind.get(last) == 'side' or last in w.sinks or dep[last]: continue
    if w.rootkind.get(last) == 'center': t = w.cut; U = w.rootframe[last]
    else: t = w.pos[w.rops[last][-1]] + 1; U = w.frames[w.rops[last][-1]]
    bynode[w.final[last]].append((t, U, last))
prs = [v for v in bynode.values() if len(v) == 2]
print('dead same-node pairs', len(prs), 'larger groups', sum(len(v) > 2 for v in bynode.values()))
C, W = w.ledger(); a = froot(C, W, m)
G = lambda r: g(r, a, m)
def gain(u0, u1, V): return 3 * (G(h - u0) + G(V) - G(V - u0) + G(h - u1) - G(V - u1) - G(h - V))
births = [(w.pos[w.rops[s][0]], w.frames[w.rops[s][0]], s) for s in fresh]
E = []
for pi, ((t0, U0, r0), (t1, U1, r1)) in enumerate(prs):
    for bi, (tb, V, s) in enumerate(births):
        if max(t0, t1) <= tb and inside(U0, V) and inside(U1, V):
            gg0 = gain(len(U0), len(U1), len(V)); gg1 = gain(len(U1), len(U0), len(V))
            gg, piv = (gg0, 0) if gg0 >= gg1 else (gg1, 1)
            if gg > 0: E.append((pi, bi, gg, piv))
print('compatible (pair, birth) edges', len(E), 'pairs with any edge', len({e[0] for e in E}))
if E:
    P_ = sorted({e[0] for e in E}); B_ = sorted({e[1] for e in E}); pi_ = {p: i for i, p in enumerate(P_)}; bi_ = {b: i for i, b in enumerate(B_)}
    gm = max(e[2] for e in E)
    r = [pi_[e[0]] for e in E] + list(range(len(P_))); c = [bi_[e[1]] for e in E] + [len(B_) + i for i in range(len(P_))]
    v = [gm - e[2] + 1 for e in E] + [gm + 1] * len(P_)
    M = csr_matrix((v, (r, c)), shape=(len(P_), len(B_) + len(P_))); rr, cc = min_weight_full_bipartite_matching(M)
    best = {(e[0], e[1]): e for e in E}; picks = [best[(P_[x], B_[y])] for x, y in zip(rr, cc) if y < len(B_)]
    C2 = Counter(C)
    for pi, bi, gg, piv in picks:
        (t0, U0, _), (t1, U1, _) = prs[pi] if piv == 0 else prs[pi][::-1]; V = len(births[bi][1]); u0, u1 = len(U0), len(U1)
        for r_, s_ in ((h - u0, -3), (V, -3), (V - u0, 3), (h - u1, -3), (V - u1, 3), (h - V, 3)):
            if r_: C2[r_] += s_
    assert all(x >= 0 for x in C2.values())
    a2 = froot(+C2, W - len(picks), m)
    print('hosts %d: a_c %.6e -> %.6e (%+.3f%%)' % (len(picks), a, a2, 100 * (a2 / a - 1)))
    print('u0,u1,V of picks', Counter((len((prs[p][piv])[1]), len((prs[p][1-piv])[1]), len(births[b][1])) for p, b, _, piv in picks).most_common(8))
# export: hosts list + formula-predicted claim (the checker recounts from chains independently)
import gzip, json
from fractions import Fraction
if E and len(sys.argv) > 2:
    J = w.J; hosts = []
    for pi, bi, gg, piv in picks:
        A = prs[pi] if piv == 0 else prs[pi][::-1]
        hosts.append(dict(pivot=A[0][2], control=A[1][2], birth=births[bi][2]))
    J['hosts'] = hosts
    J['claim_hosts'] = dict(C={str(k): int(x) for k, x in sorted((+C2).items())}, W=W - len(picks), a_float=a2)
    json.dump(J, gzip.open(sys.argv[2], 'wt')); print('wrote', sys.argv[2], len(hosts))
