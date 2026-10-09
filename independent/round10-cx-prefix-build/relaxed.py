"""#162-style relaxed carrier matching (own code; #162 read as data): our max-card matching under the pre-test, then
a greedy extension by arcs that only need (1) DAG + arc edges acyclic and (2) every value span orthogonal to its
backward-intersection annihilator. Adding arc x -> t enlarges ann for x and every node that reaches x, so the test is
span(z) orthogonal to ann(t) for all those z. Usage: compile_frames(g, relaxed_arcs) with cube_word.RELAXED = True."""
import sys, os
from collections import defaultdict, deque
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cube_word import basis
from matching import maxcard

def orth(S, A): return all(bin(a & b).count('1') % 2 == 0 for a in S for b in A)

def relaxed_arcs(D, verbose=True):
    arcs0 = maxcard(D)[0]
    order, args, uval, utgt, ucode, n, h, rann, roots, spans = (D[k] for k in
        ('order', 'args', 'uval', 'utgt', 'ucode', 'n', 'h', 'rann', 'roots', 'spans'))
    bycode = {c: u for u, c in enumerate(ucode)}
    arcs = {x: bycode[c] for x, c in arcs0}; used_u = set(arcs.values())
    active = set(order)
    succ = defaultdict(list); pred = defaultdict(list); direct = defaultdict(list)
    for x in order:
        if args[x]:
            for y in set(args[x]): succ[y].append(x); pred[x].append(y)
    for j, r in enumerate(roots): direct[r['node']].extend(rann[j])
    for x, u in arcs.items():
        t = utgt[u]
        if t < n: succ[x].append(t); pred[t].append(x)
        else: direct[x].extend(rann[t - n])
    # ann by a reverse topological pass of the current DAG + arc graph
    indeg = defaultdict(int)
    for x in active:
        for t in succ[x]: indeg[t] += 1
    q = deque(x for x in active if not indeg[x]); topo = []
    while q:
        x = q.popleft(); topo.append(x)
        for t in succ[x]:
            indeg[t] -= 1
            if not indeg[t]: q.append(t)
    assert len(topo) == len(active)
    ann = {}
    for x in reversed(topo): ann[x] = basis(list(direct[x]) + [z for t in succ[x] for z in ann[t]])
    free = defaultdict(int)                      # non-arc uses left per value
    for u, y in enumerate(uval):
        if u not in used_u: free[y] += 1
    def reaches(a, b):                           # does a reach b in the current graph
        st = [a]; seen = {a}
        while st:
            x = st.pop()
            if x == b: return True
            for t in succ[x]:
                if t not in seen: seen.add(t); st.append(t)
        return False
    def ancestors(x):
        st = [x]; seen = {x}
        while st:
            z = st.pop()
            for y in pred[z]:
                if y not in seen: seen.add(y); st.append(y)
        return seen
    added = 0
    for x in order:
        if not args[x] or x in arcs: continue
        for y in set(args[x]):
            if x in arcs: break
            if free[y] < 2: continue
            for u in [u for u, yy in enumerate(uval) if yy == y] if False else D['uses'][y]:
                if u in used_u or utgt[u] == x: continue
                t = utgt[u]; At = ann[t] if t < n else rann[t - n]
                if t < n and reaches(t, x): continue
                anc = ancestors(x)
                if not all(orth(spans[z], At) for z in anc): continue
                arcs[x] = u; used_u.add(u); free[y] -= 1; added += 1
                if t < n: succ[x].append(t); pred[t].append(x)
                else: direct[x].extend(At)
                for z in anc: ann[z] = basis(list(ann[z]) + list(At))
                break
    if verbose: print('relaxed matching: %d pre-test arcs + %d extended' % (len(arcs0), added), flush=True)
    return [(x, ucode[u]) for x, u in arcs.items()]
