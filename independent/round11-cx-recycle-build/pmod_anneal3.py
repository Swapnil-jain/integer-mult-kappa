"""Pair-module annealer acting on the module DAG itself (own code), started from eumemic's G37 module (data).
pmod_anneal2's seed-list synthesis loses the decomposition (G37 resynthesised scores 5.602e-4 at caps against 5.918e-4
for G37 itself), and the decomposition is what the carrier arcs see. Moves, all keeping every support exact:
  rewire  node x = a + b becomes y + z for two other earlier nodes with sup y | sup z = sup x (disjoint);
  split   node x = a + b becomes y + w with w a NEW node built for sup x minus sup y from earlier nodes;
then unused nodes are pruned and the module is renumbered topologically. Score: score_cp (caps, the whole cube-prefix
pipeline). Metropolis on log score. Usage: python3 pmod_anneal3.py iters seed T0"""
import sys, os, json, math, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from score_cp import score
from make_word import TRI11
iters, seed, T0 = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]); rnd = random.Random(seed); p = 11
KIND = os.environ.get('KIND', 'pair')   # pair: anneal G37; allbut: anneal the nested-prefix all-but-one module
if KIND == 'pair': M0 = json.load(open(os.path.join(HERE, 'data', 'pmod_G37_w02_5.6098194e-4.json')))
else:
    from restrict import allbut_prefix; M0 = allbut_prefix(p - 2)
if os.environ.get('START'): M0 = json.load(open(os.path.join(HERE, os.environ['START']))); M0 = M0.get('module', M0)
NI = M0['input_count']
def sups(args):
    sp = []
    for i, a in enumerate(args): sp.append(1 << i if a is None else sp[a[0]] | sp[a[1]])
    return sp
def normal(args, roots):
    """prune unused nodes, renumber (inputs first, then a topological order by the current index)"""
    act = set(range(NI)); st = list(roots)
    while st:
        x = st.pop()
        if x in act: continue
        act.add(x)
        if args[x] is not None: st.extend(args[x])
    # topological order: index order is topological except for split-created nodes appended at the end
    order = []; seen = set()
    def visit(x):
        stack = [(x, False)]
        while stack:
            y, done = stack.pop()
            if done: order.append(y); continue
            if y in seen: continue
            seen.add(y); stack.append((y, True))
            if args[y] is not None:
                for c in reversed(args[y]): stack.append((c, False))
    for x in range(NI): visit(x)
    for x in sorted(act): visit(x)
    rn = {x: i for i, x in enumerate(order)}
    return dict(input_count=NI, args=[None if args[x] is None else [rn[c] for c in args[x]] for x in order], roots=[rn[r] for r in roots])
def check(M):
    sp = sups(M['args'])
    for i, a in enumerate(M['args']):
        if a is not None: assert a[0] < i and a[1] < i and not sp[a[0]] & sp[a[1]]
    return sp
def mutate(M):
    args = [None if a is None else list(a) for a in M['args']]; roots = list(M['roots']); sp = sups(args)
    by = {}
    for i, s in enumerate(sp): by.setdefault(s, i)
    xs = [i for i, a in enumerate(args) if a is not None]
    for _ in range(200):
        x = rnd.choice(xs); S = sp[x]
        ys = [y for y in range(x) if sp[y] & S == sp[y] and sp[y] != S and y not in args[x]]
        if not ys: continue
        y = rnd.choice(ys); rest = S & ~sp[y]
        z = by.get(rest)
        if z is not None and z < x:
            args[x] = [y, z]; return normal(args, roots)
        if rnd.random() < 0.5: continue
        # split: build the remainder from the largest earlier pieces, as new nodes appended (indices > x is fine:
        # normal() re-sorts topologically; supports stay disjoint)
        pieces = []; r_ = rest
        while r_:
            c = max((c for c in range(x) if sp[c] & r_ == sp[c]), key=lambda c: bin(sp[c]).count('1'))
            pieces.append(c); r_ &= ~sp[c]
        cur = pieces[0]
        for c in pieces[1:]:
            args.append([cur, c]); sp.append(sp[cur] | sp[c]); cur = len(args) - 1
        args[x] = [y, cur]; return normal(args, roots)
    return M
cur = normal([None if a is None else list(a) for a in M0['args']], list(M0['roots'])); check(cur)
TRI = TRI11
from make_word import G37
def sc(M):
    base = dict(merge=0, tri=TRI) if not os.environ.get('H20') else dict(merge=0, trisrc=('data/h20_g1.json.gz', [0, 1, 2, 5, 7, 8, 11, 12, 14, 16, 19]))
    if KIND == 'pair': return score(p, dict(base, pairmod=M, allbut=os.environ.get('QMOD', 'prefix')))
    return score(p, dict(base, pairmod=G37, allbut=M))
cs = sc(cur); best = (cs, cur); t0 = time.time(); PROG = os.path.join(HERE, 'pa3_%s%s_s%d_best.json' % (KIND, os.environ.get('TAG', ''), seed))
print(json.dumps(dict(it=-1, cur=cs, adds=len(cur['args']) - NI)), flush=True)
for it in range(iters):
    T = T0 * (1 - it / iters) + 1e-7
    nx = mutate(cur)
    try: check(nx)
    except AssertionError: continue
    ns = sc(nx)
    if ns > 0 and (ns >= cs or rnd.random() < math.exp((math.log(ns) - math.log(cs)) / T)): cur, cs = nx, ns
    if cs > best[0]: best = (cs, cur); json.dump(dict(score=cs, module=cur), open(PROG, 'w'))
    print(json.dumps(dict(it=it, cur=cs, best=best[0], adds=len(cur['args']) - NI, t=round(time.time() - t0))), flush=True)
