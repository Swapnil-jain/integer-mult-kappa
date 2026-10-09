"""Triple-restriction local search under the cube-prefix objective (nested prefix, merged reads, G37 pair module,
relaxed arcs, pair-aware gauges + birth reuse, caps frames). SRC=pr117: labels of the PR117 h=24 DAG (restrict.tri);
SRC=h20: kept points of eumemic's h20_g1 witness (restrict.tri_from, data). Usage: python3 tsearch2.py SRC iters seed"""
import sys, json, random, time, math, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cube_word as cw
from make_word import build_cfg, G37
from relaxed import relaxed_arcs
from resel_br import reselect, with_sel, gfirst_frames, late_tau, candidates, match, recount
from restrict import tri, tri_from
src, iters, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); rnd = random.Random(seed); p = 11; m = 6 * p
NL = 24 if src == 'pr117' else 20
def froot(C, W):
    items = [(float(n), r / m) for r, n in C.items() if n]; W = float(W); lo, hi = 0.0, 0.05
    for _ in range(60):
        mid = (lo + hi) / 2; s = math.fsum(n * x ** (1 - mid) for n, x in items) / W
        lo, hi = (mid, hi) if s < 1 else (lo, mid)
    return lo
def score(L):
    try:
        cfg = dict(merge=0, pairmod=G37, allbut='prefix')
        if src == 'pr117': cfg['tri'] = L
        else: cfg['trisrc'] = ('data/h20_g1.json.gz', L)
        g = build_cfg(p, cfg); cw.RELAXED = True
        try: prof, wit = cw.compile_frames(g, lambda D: relaxed_arcs(D, verbose=False))
        finally: cw.RELAXED = False
        pr, s = cw.select_gauges(prof, wit)
        W = dict(p=p, h=2 * p, m=m, g=g, prof=prof, wit=wit, pr=pr, s=s)
        F0 = gfirst_frames(W); sel, _, _ = reselect(W, F0, 5.6e-4, set(), True)
        W2 = with_sel(W, sel); F = gfirst_frames(W2); tau = late_tau(W2, F)
        pairs = match(candidates(W2, F, tau, 5.6e-4, True, True))
        C, Wv, inf = recount(W2, F, pairs, tau, True)
        return 0.0 if inf['bad'] else froot(C, Wv)
    except (AssertionError, TypeError, KeyError, IndexError): return 0.0
if seed == 1: cur = [1, 3, 4, 7, 8, 11, 13, 17, 19, 20, 23] if src == 'pr117' else [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 16]
else: cur = sorted(rnd.sample(range(NL), p))
cs = score(cur); best = (cs, cur); t0 = time.time()
print(json.dumps(dict(it=-1, cur=cs, best=cs, L=cur)), flush=True)
for it in range(iters):
    nx = list(cur); free = [x for x in range(NL) if x not in nx]
    nx[rnd.randrange(p)] = rnd.choice(free); nx.sort(); ns = score(nx)
    if ns >= cs: cur, cs = nx, ns
    if cs > best[0]: best = (cs, cur)
    print(json.dumps(dict(it=it, cur=cs, best=best[0], L=best[1], t=round(time.time() - t0))), flush=True)
