"""Load a frozen addition DAG (PR #117's complex producer, read as data: args pairs, D disjoint outputs, P pair-star
outputs, A centres) into the NStar3 interface, recomputing and checking every support with our own code.
Node j+1 is input triple j (combinations order); node v+1+k = args[2k] + args[2k+1]."""
import gzip, json
from itertools import combinations


class DagProducer:
    def __init__(s, path):
        d = json.load(gzip.open(path))
        h = s.h = d['h']; s.allE = True
        s.triples = list(combinations(range(h), 3)); s.tid = {t: i for i, t in enumerate(s.triples)}
        v = len(s.triples); assert v == d['v']
        s.sup = [0] + [1 << i for i in range(v)]; s.args = [None] * (v + 1)
        a = d['args']; assert len(a) % 2 == 0
        for k in range(len(a) // 2):
            x, y = a[2 * k], a[2 * k + 1]; n = len(s.sup)
            assert 0 < x < n and 0 < y < n, 'args must be earlier nodes'
            assert not s.sup[x] & s.sup[y], 'cancellation-free'
            s.sup.append(s.sup[x] | s.sup[y]); s.args.append((x, y))
        s.cov = {}
        pm = [sum(1 << j for j, t in enumerate(s.triples) if p in t) for p in range(h)]
        full = (1 << v) - 1
        s.pieces = []
        for j, (S, n) in enumerate(zip(s.triples, d['D'])):
            assert s.sup[n] == full & ~(pm[S[0]] | pm[S[1]] | pm[S[2]]), 'disjoint output'
            s.pieces.append((S, n, 1))
        seen = set()
        for n in d['P']:
            sp = s.sup[n]; pts = [p for p in range(h) if sp & pm[p] == sp]     # points in every triple
            assert len(pts) == 2, 'pair-star output'
            a_, b_ = pts; rest = [i for i in range(h) if i not in pts and not (sp & pm[i] & pm[a_] & pm[b_]) ]
            # the missing third point i: {a,b,i} not in support
            miss = [i for i in range(h) if i not in pts and not sp >> s.tid[tuple(sorted((a_, b_, i)))] & 1]
            assert len(miss) == 1; S = tuple(sorted((a_, b_, miss[0])))
            assert sp == (pm[a_] & pm[b_]) & ~(1 << s.tid[S]); assert (a_, b_, miss[0]) not in seen
            seen.add((a_, b_, miss[0])); s.pieces.append((S, n, -1))
        assert len(seen) == 3 * v
        s.retained = []
        for i, n in enumerate(d['A']):
            assert s.sup[n] == full & ~pm[i], 'centre E_i'
            s.retained.append((('E', i), n))
        st = [n for _, n, _ in s.pieces] + [n for _, n in s.retained]; s.active = set()
        while st:
            n = st.pop()
            if n in s.active: continue
            s.active.add(n)
            if s.args[n]: st.extend(s.args[n])
        s.additions = sum(1 for n in s.active if s.args[n])
        s.roles = s.additions + len(s.pieces) + len(s.retained)
    def x(s, *p): return s.tid[tuple(sorted(p))] + 1
    def cover(s, n):
        if n not in s.cov:
            sp = s.sup[n]; m = 0
            while sp:
                lo = sp & -sp; m |= sum(1 << p for p in s.triples[lo.bit_length() - 1]); sp ^= lo
            s.cov[n] = m
        return s.cov[n]


if __name__ == '__main__':
    import sys
    c = DagProducer(sys.argv[1]); print(dict(h=c.h, additions=c.additions, pieces=len(c.pieces), roles=c.roles))
