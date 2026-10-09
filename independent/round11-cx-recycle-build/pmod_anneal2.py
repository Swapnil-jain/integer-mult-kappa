"""Own annealer for the pair-disjoint module at n = p - 1 labels (output for pair q = sum of x_e over edges e of K_n
avoiding q), scored on the cover itself (caps-level float a*, merged reads, searched triple restriction), as the
coordinator suggested; smallh-search showed that fewest additions is not the objective.

State: an ordered list of seed label sets U; seed U denotes the sub-sum over edges avoiding U. The module is
synthesised greedily: seeds first (in order), then the 45 roots; each support S is built from the largest pool node
inside what is left of S, repeatedly, and the pieces are summed left to right (every new partial sum joins the pool).
All supports are exact sub-sums of exclusion sums, so every node is frame-legal. Moves: add a seed (|U| = 3..6, or the
union of two roots' pairs), drop a seed, replace a seed's label, swap two seeds. Metropolis on log a*.
Usage: python3 pmod_anneal.py p iters seed [T0]

cube-prefix version (pmod_anneal2.py): scored on the stronger objective -- nested-prefix all-but-one, merged reads,
searched triple, RELAXED carrier arcs, pair-aware gauges and birth reuse (exact matching + chain recount), caps frames
(no descent, for speed) -- and started from eumemic's G37 module (its internal supports as the seed list).
Env SCORE=desc adds node descent to every score (slow)."""
import sys, os, json, math, random, time
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cube_word as cw
from matching import maxcard
from make_word import build_cfg, CFG

p, iters, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); T0 = float(sys.argv[4]) if len(sys.argv) > 4 else 0.004
n = p - 1; E = list(combinations(range(n), 2)); NE = len(E); rnd = random.Random(seed); m = 6 * p
TRI = CFG['tri%d' % p]['tri'] if ('tri%d' % p) in CFG else None
def avoid(U): return sum(1 << k for k, e in enumerate(E) if not set(e) & set(U))

def synth(seeds):
    args = [None] * NE; sup = [1 << k for k in range(NE)]; by = {s: i for i, s in enumerate(sup)}
    def build(S):
        if S in by: return by[S]
        pieces = []; rest = S
        while rest:
            # largest pool node inside rest (pool is small enough for a linear scan)
            best = max((i for i, s in enumerate(sup) if s & rest == s), key=lambda i: bin(sup[i]).count('1'))
            pieces.append(best); rest &= ~sup[best]
        cur = pieces[0]
        for b in pieces[1:]:
            s = sup[cur] | sup[b]
            if s in by: cur = by[s]; continue
            by[s] = len(args); args.append([cur, b]); sup.append(s); cur = by[s]
        return cur
    for S in seeds:
        if S: build(S)
    roots = [build(avoid(q)) for q in E]
    act = set(range(NE)); st = list(roots)
    while st:
        x = st.pop()
        if x in act: continue
        act.add(x)
        if args[x] is not None: st.extend(args[x])
    ids = sorted(act); rn = {x: i for i, x in enumerate(ids)}
    return dict(input_count=NE, args=[None if args[x] is None else [rn[y] for y in args[x]] for x in ids], roots=[rn[x] for x in roots])

def froot(C, W):
    items = [(float(k), r / m) for r, k in C.items() if k]; W = float(W); lo, hi = 0.0, 0.05
    for _ in range(60):
        mid = (lo + hi) / 2; s = math.fsum(k * x ** (1 - mid) for k, x in items) / W
        lo, hi = (mid, hi) if s < 1 else (lo, mid)
    return lo
def score(M):
    try:
        from resel_br import reselect, with_sel, gfirst_frames, late_tau, candidates, match, recount
        from relaxed import relaxed_arcs
        cfg = dict(merge=0, pairmod=M, allbut='prefix')
        if TRI: cfg['tri'] = TRI
        g = build_cfg(p, cfg); cw.RELAXED = True
        try: prof, wit = cw.compile_frames(g, lambda D: relaxed_arcs(D, verbose=False))
        finally: cw.RELAXED = False
        if os.environ.get('SCORE') == 'desc':
            from descent import run as descend, apply as dapply
            ann, rank, H = descend(prof, wit, sweeps=4, verbose=False); prof, wit = dapply(prof, wit, ann, rank, H)
        pr, s = cw.select_gauges(prof, wit)
        W = dict(p=p, h=2 * p, m=m, g=g, prof=prof, wit=wit, pr=pr, s=s)
        F0 = gfirst_frames(W); sel, _, _ = reselect(W, F0, 5.6e-4, set(), True)
        W2 = with_sel(W, sel); F = gfirst_frames(W2); tau = late_tau(W2, F)
        pairs = match(candidates(W2, F, tau, 5.6e-4, True, True))
        C, Wv, inf = recount(W2, F, pairs, tau, True)
        if inf['bad']: return 0.0
        return froot(C, Wv)
    except (AssertionError, KeyError, IndexError, TypeError): return 0.0

def rand_seed(S=None):
    r = rnd.random()
    if r < 0.4:
        q1, q2 = rnd.sample(E, 2); U = sorted(set(q1) | set(q2))
        if len(U) < 3: U = sorted(set(U) | {rnd.choice([x for x in range(n) if x not in U])})
        return avoid(U)
    if r < 0.7 or not S or len(S) < 2: return avoid(rnd.sample(range(n), rnd.randint(3, 6)))
    for _ in range(50):                      # union of two disjoint current seeds inside some root support
        a, b = rnd.sample(S, 2)
        if not a & b and any((a | b) & ~avoid(q) == 0 for q in E): return a | b
    return avoid(rnd.sample(range(n), rnd.randint(3, 6)))
def mutate(S):
    S = list(S); r = rnd.random()
    if r < 0.35 or not S: S.insert(rnd.randint(0, len(S)), rand_seed(S))
    elif r < 0.6: S.pop(rnd.randrange(len(S)))
    elif r < 0.85:                            # move one seed to a random new position
        i = rnd.randrange(len(S)); x = S.pop(i); S.insert(rnd.randint(0, len(S)), x)
    elif len(S) > 1:
        i, j = rnd.sample(range(len(S)), 2); S[i], S[j] = S[j], S[i]
    return S
# start: the internal node supports of a PR117 restriction (default labels; other chains a random restriction)
from restrict import pair as _pair
def supports(M):
    sp = []
    for i, a in enumerate(M['args']):
        sp.append(1 << i if a is None else sp[a[0]] | sp[a[1]])
    return [sp[i] for i, a in enumerate(M['args']) if a is not None]
base = json.load(open(os.path.join(HERE, 'data', 'pmod_G37_w02_5.6098194e-4.json')))
cur = supports(base)
Mc = synth(cur); cs = score(Mc); best = (cs, cur, Mc); t0 = time.time()
print(json.dumps(dict(it=-1, cur=cs, best=cs, adds=sum(a is not None for a in Mc['args']))), flush=True)
PROG = os.path.join(HERE, 'pmod2_p%d_s%d_best.json' % (p, seed))
for it in range(iters):
    T = T0 * (1 - it / iters) + 1e-6
    nx = mutate(cur); Mn = synth(nx); ns = score(Mn)
    if ns > 0 and (ns >= cs or rnd.random() < math.exp((math.log(ns) - math.log(cs)) / T)): cur, cs = nx, ns
    if cs > best[0]:
        best = (cs, cur, synth(cur)); json.dump(dict(score=best[0], seeds=best[1], module=best[2]), open(PROG, 'w'))
    print(json.dumps(dict(it=it, cur=cs, best=best[0], nseeds=len(cur), t=round(time.time() - t0))), flush=True)
json.dump(dict(score=best[0], seeds=best[1], module=best[2]), open(os.path.join(HERE, 'pmod2_p%d_s%d.json' % (p, seed)), 'w'))
