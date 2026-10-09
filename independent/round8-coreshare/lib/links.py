"""Carrier links (maximum matching), the complex analogue of PR #41 links, as PR #104 uses them ("carrier
matching", PR #24 lineage; read as data, reimplemented here).

Time order of the gates: (phase, label dim, id), a topological order since an argument's label lies in its
node's label. A donor gate x reads two argument slots; after x, the NON-pivot argument slot still holds its value a.
A link hands that slot to a later use u of a whose frame contains x's frame (the slot moves M_x -> M_u, a nested
step), so u needs no fresh copy. If the linked value is args[x][0], x takes args[x][1]'s slot as its pivot
(orientation change). Each donor gives one slot, each use receives at most one: maximum bipartite matching
(Hopcroft-Karp). Slots: R = #additions + #roots - #links."""
from collections import defaultdict, deque
PIVOT_SMALL_DEAD = [False]


def hopcroft_karp(adj, nl):
    """maximum bipartite matching; adj[u] = right vertices of left u. Returns ml (right vertex or -1)."""
    ml = [-1] * nl; mr = {}
    while True:
        dist = {}; q = deque()
        for u in range(nl):
            if ml[u] < 0: dist[u] = 0; q.append(u)
        found = False
        while q:
            u = q.popleft()
            for v in adj[u]:
                w = mr.get(v)
                if w is None: found = True
                elif w not in dist: dist[w] = dist[u] + 1; q.append(w)
        if not found: return ml
        its = {}
        for u0 in range(nl):
            if ml[u0] >= 0 or dist.get(u0) != 0: continue
            stack = [u0]; vs = []
            while stack:
                x = stack[-1]
                if x not in its: its[x] = iter(adj[x])
                nxt = None; aug = False
                for v in its[x]:
                    w = mr.get(v)
                    if w is None: vs.append(v); aug = True; break
                    if dist.get(w) == dist[x] + 1: nxt = (v, w); break
                if aug:
                    for k, y in enumerate(stack): ml[y] = vs[k]; mr[vs[k]] = y
                    break
                if nxt: vs.append(nxt[0]); stack.append(nxt[1])
                else:
                    dist[x] = -1; stack.pop()
                    if vs: vs.pop()


def compile_links(c, lab, X, tperp, tkey, verbose=True):
    args = c.args
    uses = defaultdict(list)                     # value -> [use]; use = ('g', g, pos) | ('p', i) | ('r', name)
    for g in c.active:
        if args[g]:
            for pos, a in enumerate(args[g]): uses[a].append(('g', g, pos))
    for i, (S, n, cf) in enumerate(c.pieces): uses[n].append(('p', i))
    for nm, n in c.retained: uses[n].append(('r', nm))
    def utime(u):
        if u[0] == 'g': return tkey(u[1])
        if u[0] == 'r': return (1.5, 0, 0)
        return (9, 0, u[1])
    def uframe(u):
        return lab[u[1]] if u[0] == 'g' else tperp[c.pieces[u[1]][0]]
    donors = [x for x in sorted(c.active, key=tkey) if args[x]]
    uid = {}; ulist = []
    adj = []
    for x in donors:
        e = []
        for pos, a in enumerate(args[x]):
            for u in uses[a]:
                if u[0] == 'r' or u == ('g', x, pos): continue
                if not (tkey(x) < utime(u)): continue
                F = uframe(u)
                if X.le(lab[x], F) and X.step_ok(lab[x], F):
                    key = (u, pos)
                    if u not in uid: uid[u] = len(ulist); ulist.append(u)
                    e.append(uid[u] * 2 + pos)        # encode the arg position with the use
        adj.append(e)
    # a use may be targeted through either position only once: match on use ids (strip pos)
    adj_u = [[k // 2 for k in e] for e in adj]
    ml = hopcroft_karp(adj_u, len(donors))
    link = {}                                    # target use -> (donor x, pos of the linked value in x)
    for i, x in enumerate(donors):
        if ml[i] >= 0:
            u = ulist[ml[i]]
            pos = next(k % 2 for k in adj[i] if k // 2 == ml[i])
            link[u] = (x, pos)
    # unlinked donors: the larger argument stays in place, so the slot left behind holds the smaller signal
    # (easier to clear exactly by dependence reclamation); PIVOT_SMALL_DEAD = False keeps the original order
    piv = {x: (1 if PIVOT_SMALL_DEAD[0] and bin(c.sup[args[x][1]]).count('1') > bin(c.sup[args[x][0]]).count('1') else 0) for x in donors}
    for u, (x, pos) in link.items():
        piv[x] = 1 - pos                         # the linked value's slot must be the non-pivot one
    # slot assignment
    slot = {}; size = 0; gates = []; src = {}; pout = {}; rout = {}
    order = sorted(c.active, key=tkey)
    fresh = {}
    for a in order:
        us = sorted(uses[a], key=utime)
        fresh[a] = [u for u in us if u not in link]
    def slot_of(u):
        if u in slot: return slot[u]
        x, pos = link[u]
        s_ = slot_of(('g', x, pos)); slot[u] = s_; return s_
    for n in order:
        if args[n]:
            p = piv[n]; ins = (slot_of(('g', n, p)), slot_of(('g', n, 1 - p))); pv = ins[0]
        else:
            pv = size; size += 1; ins = (pv,); src[c.triples[n - 1]] = pv
        fr = fresh[n]
        outs = [pv]
        slot[fr[0]] = pv
        for u in fr[1:]:
            slot[u] = size; outs.append(size); size += 1
        gates.append((n, ins, tuple(outs)))
    for u in list(link): slot_of(u)
    for a in order:
        for u in uses[a]:
            if u[0] == 'p': pout[u[1]] = slot[u]
            elif u[0] == 'r': rout[u[1]] = slot[u]
    nadd = sum(1 for n in c.active if args[n])
    assert size == nadd + len(c.pieces) + len(c.retained) - len(link), (size, nadd, len(link))
    if verbose: print('links', len(link), 'R', size, flush=True)
    return dict(size=size, gates=gates, src=src, pout=pout, rout=rout, links=link, piv=piv)
