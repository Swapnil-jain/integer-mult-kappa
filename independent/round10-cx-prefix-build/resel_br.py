"""Pair-aware gauge re-selection. #144's chronological gauge filter (same candidate order, same target-limit chain)
but a candidate that can be born on a dead donor's slot is charged the hand-off instead of its tail:
  paired:   3[ex(r-d) - ex(r)] + 3[ex(d-e) - ex(h-e)] + target splits      (W -1, D unchanged)
  unpaired: 3[ex(r-d) - ex(r)] + (3/t) ex(m - t u) + target splits         (NDS tail; #144 tail 3d if nds off)
Donors are reserved greedily during the scan (largest last frame inside sigma, latest death before the candidate's
first op); the final pairs come from the exact max-weight matching in cbirth, then the recount."""
import sys, os, bisect, pickle
import numpy as np
from collections import defaultdict, Counter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cbirth import *
from cube_word import basis, perp

DBG = {int(t) for t in os.environ.get("DBG", "").split(",") if t}
def compat_table(F, R, gsig_excl=()):
    rops, frames, pos = F['rops'], F['frames'], F['pos']
    elig = [sr for sr in range(R) if rops[sr] and sr not in F['rootframe']]
    byE = defaultdict(list)
    for sr in elig: byE[frames[rops[sr][-1]]].append((pos[rops[sr][-1]], sr))
    Es = list(byE); [byE[E].sort() for E in Es]
    return Es, byE

def reselect(W, F, a0, nopair=frozenset(), nds=True, use_pairs=True, reserve=os.environ.get('RESERVE', '0') == '1'):
    G = W['wit']['g']; h, m = W['h'], W['m']; pr, s = W['pr'], W['s']; v = pr['v']; R = pr['R']
    roots = G['roots']; inputs = G['inputs']; ann = W['wit']['ann']; ops = s['ops']; first = s['first']
    touched = s['touched']; pos, rops = F['pos'], F['rops']
    co = [0] * R
    for r, sr in zip(roots, s['rootroles']): co[sr] |= sum(1 << t for t in r['targets'])
    for a, b, x in reversed(ops): co[b] |= co[a]
    limit = [None] * v
    for r in roots:
        if r['kind'] == 'center': continue
        A = basis(inputs[t] for t in r['targets'])
        for t in r['targets']:
            if limit[t] is None: limit[t] = A
    for t in range(v):
        if limit[t] is None: limit[t] = (inputs[t],)
    cands = sorted((x for x in range(R) if x not in touched and first[x] is not None),
                   key=lambda x: (len(ann[first[x]]), bin(co[x]).count("1"), x))
    # donor classes by exact last frame; compatibility E <= perp(A) via a parity table over all vectors of A's seen
    Es, byE = compat_table(F, R)
    md = max(len(E) for E in Es); emat = np.zeros((len(Es), md), dtype=np.int64)
    for i, E in enumerate(Es): emat[i, :len(E)] = E
    edim = np.array([len(E) for E in Es])
    avail = {i: list(byE[E]) for i, E in enumerate(Es)}
    vcache = {}
    def okvec(a):
        if a not in vcache:
            vcache[a] = ~np.any(np.bitwise_count(emat & a) & 1, axis=1)
        return vcache[a]
    E_order_cache = {}
    sel = []; reserved = {}; npair = 0; nun = 0; rej = 0
    for x in cands:
        A = tuple(ann[first[x]]); tg = []; bits = co[x]
        while bits:
            low = bits & -bits; t = low.bit_length() - 1; bits ^= low; tg.append(t)
            A = basis(A + tuple(limit[t]))
            if len(A) == h: break
        if len(A) == h: continue
        r = h - len(ann[first[x]]); d = h - len(A)
        base = 3 * (ex(r - d, a0, m) - ex(r, a0, m))
        base += 3 * sum(ex(d, a0, m) + ex(h - len(limit[t]) - d, a0, m) - ex(h - len(limit[t]), a0, m) for t in tg)
        z = dict(role=x, A=A, d=d, targets=tg, r=r)
        if nds:
            t_, u = nds_t(z, h, m); unp = 3 / t_ * ex(m - t_ * u, a0, m)
        else: unp = ex(3 * d, a0, m)
        best = (base + unp, None)
        if use_pairs and x not in nopair and rops[x]:
            dl = pos[rops[x][0]]
            ok = np.ones(len(Es), dtype=bool)
            for a in A: ok &= okvec(a)
            idx = np.nonzero(ok)[0]; idx = idx[np.argsort(-edim[idx], kind='stable')]
            for i in idx:
                L = avail[i]; k = bisect.bisect_left(L, (dl, -1)) - 1
                while k >= 0 and L[k][1] == x: k -= 1
                if k < 0: continue
                e = int(edim[i]); val = base + 3 * (ex(d - e, a0, m) - ex(h - e, a0, m))
                if val < best[0]: best = (val, (i, k))
                break
        if DBG and x in DBG: print("DBG", x, r, d, base/a0, unp/a0, best[0]/a0, best[1], int(ok.sum()) if use_pairs else None, flush=True)
        if best[0] >= -1e-12: rej += 1; continue
        if best[1] is not None:
            i, k = best[1]; dd, sr = avail[i].pop(k) if reserve else avail[i][k]; reserved[x] = sr; npair += 1
        else: nun += 1
        for t in tg: limit[t] = A
        sel.append(z)
    return sel, reserved, dict(paired=npair, unpaired=nun, rejected=rej, cands=len(cands))

def with_sel(W, sel):
    W2 = dict(W); W2['s'] = dict(W['s'], sel=sel); return W2

def run(p, a0, nds=True, iters=2, tag='', name=None):
    if name:
        from make_word import make
        W = make(p, name); tag = '_' + name + tag
    else: W = load_word(p)
    F0 = gfirst_frames(W)
    for nd in (False, True):
        C, Wv, inf0 = recount(W, F0, [], late_tau(W, F0), nd)
        a = evaluate(W, C, Wv)
        log('p=%d %s base (descent, #144 gauges, no pairs) nds=%s: R %d gauges %d bad %s D %s a_c %s = %.10e' % (p, name, nd,
            W['pr']['R'], len(W['s']['sel']), inf0['bad'], D_of(C, Wv, W['m']), a, float(a)))
    out = {}
    nopair = set()
    for it in range(iters):
        sel, res, info = reselect(W, F0, a0, nopair, nds)
        W2 = with_sel(W, sel); F = gfirst_frames(W2); tau = late_tau(W2, F)
        ed = candidates(W2, F, tau, a0, nds, True); pairs = match(ed)
        C, Wv, inf2 = recount(W2, F, pairs, tau, nds)
        a = evaluate(W2, C, Wv)
        rec = {b for _, b in pairs}
        log('p=%d it %d: gauges %d (greedy paired %d, unpaired %d, rejected %d of %d cands), matched pairs %d, chained %d, '
            'bad %s D %s W %s a_c %s = %.10e' % (p, it, len(sel), info['paired'], info['unpaired'], info['rejected'],
            info['cands'], len(pairs), len({x for x, _ in pairs} & rec), inf2['bad'], D_of(C, Wv, W['m']), Wv, a, float(a)))
        out[it] = dict(sel=sel, pairs=pairs, tau=tau, a=a, C=dict(C), W=Wv, deg=dict(inf2['deg']))
        nopair |= set(res) - rec      # greedy assumed a donor the exact matching could not give
    pickle.dump(dict(p=p, a0=a0, nds=nds, out=out), open(os.path.join(HERE, 'resel_p%d%s.pkl' % (p, tag)), 'wb'))
    return out

if __name__ == '__main__':
    p = int(sys.argv[1]); a0 = float(sys.argv[2]) if len(sys.argv) > 2 else 5.6e-4
    nds = (sys.argv[3] != '0') if len(sys.argv) > 3 else True
    run(p, a0, nds, tag=os.environ.get('TAG', ''), name=sys.argv[4] if len(sys.argv) > 4 else None)
