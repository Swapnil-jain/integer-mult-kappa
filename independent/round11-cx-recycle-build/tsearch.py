"""Local search of the triple-exclusion restriction (p labels of the pinned DAG) at p = 11, 13, default pair
restriction, merged reads; score = caps-level float a* (#144 gauge rules, no NDS), as smallh-search's ctxtri.
Usage: python3 tsearch.py p iters seed"""
import sys, json, random, time, math, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cube_word as cw
from matching import maxcard
from make_word import build_cfg
p, iters, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); rnd = random.Random(seed); m = 6 * p
def froot(C, W):
    items = [(float(n), r / m) for r, n in C.items() if n]; W = float(W); lo, hi = 0.0, 0.05
    for _ in range(60):
        mid = (lo + hi) / 2; s = math.fsum(n * x ** (1 - mid) for n, x in items) / W
        lo, hi = (mid, hi) if s < 1 else (lo, mid)
    return lo
def score(L):
    try:
        g = build_cfg(p, dict(merge=0, tri=L)); prof, wit = cw.compile_frames(g, lambda D: maxcard(D)[0])
        pr, s = cw.select_gauges(prof, wit); C = cw.children(pr); return froot(C, 2 * pr['v'] + pr['R'])
    except (AssertionError, TypeError, KeyError, IndexError): return 0.0
cur = list(range(p)) if seed == 1 else sorted(rnd.sample(range(24), p))
cs = score(cur); best = (cs, cur); t0 = time.time()
print(json.dumps(dict(it=-1, cur=cs, best=cs, L=cur)), flush=True)
for it in range(iters):
    nx = list(cur); free = [x for x in range(24) if x not in nx]
    nx[rnd.randrange(p)] = rnd.choice(free); nx.sort(); ns = score(nx)
    if ns >= cs: cur, cs = nx, ns
    if cs > best[0]: best = (cs, cur)
    print(json.dumps(dict(it=it, cur=cs, best=best[0], L=best[1], t=round(time.time() - t0))), flush=True)
