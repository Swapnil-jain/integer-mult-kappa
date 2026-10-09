"""Own carrier matching for the cube word (replaces #144's frozen arcs): eligible arcs donor x -> use u of an operand
value y of x, with u's target later than x and initial[x] contained in u's frame (as #144's compiler requires);
each donor <= 1 arc, each use <= 1 arc, each value keeps >= 1 non-arc use. Max-cardinality by Dinic max-flow, or
max-weight by min-cost flow (weights: a callable on the arc)."""
import sys
from collections import deque
sys.path.insert(0, '.')
from cube_word import contained
def eligible(D):
    order, pos, args, uval, utgt, ufr, initial, n = (D[k] for k in ('order', 'pos', 'args', 'uval', 'utgt', 'ufr', 'initial', 'n'))
    byval = {}
    for u, y in enumerate(uval): byval.setdefault(y, []).append(u)
    E = []
    for x in order:
        if not args[x]: continue
        for y in set(args[x]):
            for u in byval.get(y, []):
                t = utgt[u]
                if t < n and pos[t] <= pos[x]: continue
                if contained(initial[x], ufr[u]): E.append((x, u))
    return E, byval
class Dinic:
    def __init__(s, N): s.N = N; s.g = [[] for _ in range(N)]
    def add(s, a, b, c):
        s.g[a].append([b, c, len(s.g[b])]); s.g[b].append([a, 0, len(s.g[a]) - 1])
    def flow(s, S, T):
        F = 0
        while True:
            lv = [-1] * s.N; lv[S] = 0; q = deque([S])
            while q:
                a = q.popleft()
                for b, c, _ in s.g[a]:
                    if c and lv[b] < 0: lv[b] = lv[a] + 1; q.append(b)
            if lv[T] < 0: return F
            it = [0] * s.N
            def dfs(a, f):
                if a == T: return f
                while it[a] < len(s.g[a]):
                    e = s.g[a][it[a]]; b, c, r = e
                    if c and lv[b] == lv[a] + 1:
                        d = dfs(b, min(f, c))
                        if d: e[1] -= d; s.g[b][r][1] += d; return d
                    it[a] += 1
                return 0
            while True:
                f = dfs(S, 1 << 30)
                if not f: break
                F += f
def maxcard(D):
    sys.setrecursionlimit(100000)
    E, byval = eligible(D); ucode = D['ucode']
    donors = sorted({x for x, _ in E}); usesE = sorted({u for _, u in E}); vals = sorted({D['uval'][u] for u in usesE})
    di = {x: i + 2 for i, x in enumerate(donors)}; ui = {u: i + 2 + len(donors) for i, u in enumerate(usesE)}
    vi = {y: i + 2 + len(donors) + len(usesE) for i, y in enumerate(vals)}
    G = Dinic(2 + len(donors) + len(usesE) + len(vals))
    for x in donors: G.add(0, di[x], 1)
    for x, u in E: G.add(di[x], ui[u], 1)
    for u in usesE: G.add(ui[u], vi[D['uval'][u]], 1)
    for y in vals: G.add(vi[y], 1, len(byval[y]) - 1)
    F = G.flow(0, 1)
    arcs = []
    for x in donors:
        for b, c, _ in G.g[di[x]]:
            if b in ui.values() and c == 0 and b != 0:
                pass
    inv = {i: u for u, i in ui.items()}
    for x in donors:
        for b, c, _ in G.g[di[x]]:
            if b in inv and c == 0: arcs.append((x, ucode[inv[b]]))
    assert len(arcs) == F
    return arcs, len(E)

import heapq, math
def mincost(D, weight, scale=10**6):
    """max-weight matching under the same constraints: min-cost flow, successive shortest paths (Dijkstra+potentials),
    stopping when the next augmenting path has nonnegative cost. weight(x, u) -> float gain (only > 0 arcs kept)."""
    E, byval = eligible(D); ucode = D['ucode']
    E = [(x, u, int(round(weight(x, u) * scale))) for x, u in E]; E = [e for e in E if e[2] > 0]
    donors = sorted({x for x, _, _ in E}); usesE = sorted({u for _, u, _ in E}); vals = sorted({D['uval'][u] for u in usesE})
    di = {x: i + 2 for i, x in enumerate(donors)}; ui = {u: i + 2 + len(donors) for i, u in enumerate(usesE)}
    vi = {y: i + 2 + len(donors) + len(usesE) for i, y in enumerate(vals)}
    N = 2 + len(donors) + len(usesE) + len(vals); g = [[] for _ in range(N)]
    def add(a, b, c, w):
        g[a].append([b, c, w, len(g[b])]); g[b].append([a, 0, -w, len(g[a]) - 1])
    for x in donors: add(0, di[x], 1, 0)
    for x, u, w in E: add(di[x], ui[u], 1, -w)
    for u in usesE: add(ui[u], vi[D['uval'][u]], 1, 0)
    for y in vals: add(vi[y], 1, len(byval[y]) - 1, 0)
    # initial potentials: Bellman-Ford on the DAG-layered graph (source->donor->use->value->sink)
    pot = [0] * N
    for x in donors: pot[di[x]] = 0
    for x, u, w in E: pot[ui[u]] = min(pot[ui[u]], -w)
    for u in usesE: pot[vi[D['uval'][u]]] = min(pot[vi[D['uval'][u]]], pot[ui[u]])
    pot[1] = min([pot[vi[y]] for y in vals] + [0])
    total = 0; flow = 0
    while True:
        dist = [None] * N; dist[0] = 0; prev = [None] * N; pq = [(0, 0)]
        while pq:
            d, a = heapq.heappop(pq)
            if d > dist[a]: continue
            for k, (b, c, w, r) in enumerate(g[a]):
                if c <= 0: continue
                nd = d + w + pot[a] - pot[b]
                if dist[b] is None or nd < dist[b]: dist[b] = nd; prev[b] = (a, k); heapq.heappush(pq, (nd, b))
        if dist[1] is None: break
        real = dist[1] - pot[0] + pot[1]
        if real >= 0: break
        for i in range(N):
            if dist[i] is not None: pot[i] += dist[i]
        b = 1
        while b != 0:
            a, k = prev[b]; e = g[a][k]; e[1] -= 1; g[b][e[3]][1] += 1; b = a
        total += real; flow += 1
    inv = {i: u for u, i in ui.items()}; arcs = []
    for x in donors:
        for b, c, w, _ in g[di[x]]:
            if b in inv and c == 0: arcs.append((x, ucode[inv[b]]))
    return arcs, flow, -total / scale
