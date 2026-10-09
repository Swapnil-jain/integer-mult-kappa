"""Independent recount of #144's p=12 complex aux histogram by SIMULATING the physical word role by role.

compile_frames() reproduces #144's histogram by a per-node ledger formula. Here instead every one of the R physical
roles is followed chronologically through the op list built by select_gauges (dest += ctl gates, fan-out copies,
carried arcs, root reads), recording the actual F2 subspace it sits at; each transition must be NESTED (checked
on subspaces, not dimensions) and costs dim(new) - dim(old). Legality checks: every gate frame contains the span of
the value computed there; both gate registers share the identical frame; side root reads sit at the target cap
q_T^perp (which must contain the role's current frame); gauged roles enter at sigma = A^perp contained in their first
frame. The resulting histogram, plus sources/targets/gauges, is compared with the published certificate."""
import json, sys, pickle
from collections import Counter
import os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cube_word import build, compile_frames, select_gauges, published, basis, perp, contained, children
from evalc import acert, D_of

def run():
    g = build(12); frozen = json.load(open(os.path.join(HERE, '..', '..', 'certificates', 'round9', 'pr144', 'matching-arcs.json')))
    prof, wit = compile_frames(g, frozen); pr, s = select_gauges(prof, wit)
    G = wit['g']; h = 24; v = 1760; n = wit['n']; args = G['args']; roots = G['roots']; inputs = G['inputs']
    ann = wit['ann']; U = {}
    def frame(x):
        if x not in U: U[x] = perp(ann[x], h)
        return U[x]
    # spans computed afresh (support legality)
    spans = [()] * n
    for i, u in enumerate(inputs): spans[i + 1] = (u,)
    for x in range(v + 1, n):
        if args[x]: spans[x] = basis(spans[args[x][0]] + spans[args[x][1]])
    R = pr['R']; chain = [[()] for _ in range(R)]; bad = Counter()
    gsig = {z['role']: perp(z['A'], h) for z in s['sel']}
    srcrole = {r: x for x, r in s['sources'].items()}
    for r, x in srcrole.items(): chain[r] += [(inputs[x - 1],), frame(x)]   # load source at <q_S> (rank 1), climb to its node frame
    for (a, b, x) in s['ops']:
        Ux = frame(x)
        if not contained(spans[x], Ux): bad['support'] += 1
        for role in (a, b):
            if len(chain[role]) == 1 and role in gsig: chain[role].append(gsig[role])   # entrance gauge
            chain[role].append(Ux)
    loss = 0; centre_copies = []
    for j, (rt, rr) in enumerate(zip(roots, s['rootroles'])):
        x = rt['node']
        if rt['kind'] == 'center':
            Uc = spans[x]
            if basis(Uc) != basis(frame(x)): bad['center frame'] += 1
            chain[rr].append(frame(x)); loss += len(Uc); centre_copies.append(len(Uc))   # copied-centre transform pays dim U_c
        else:
            cap = perp(basis(inputs[t] for t in rt['targets']), h)
            chain[rr].append(cap)
    full = tuple(1 << j for j in range(h - 1, -1, -1))
    H = Counter();
    for r in range(R):
        ch = chain[r] + [full]
        for A_, B_ in zip(ch, ch[1:]):
            if not contained(A_, B_): bad['retreat'] += 1
            d = len(basis(B_)) - len(basis(A_))
            if r in gsig and A_ == () and B_ == gsig[r]: continue          # gauge entrance = the 3d tail, not internal
            H[d] += 1
        if len(basis(ch[-2])) == 0 and r not in gsig: bad['never moved'] += 1
    for d in centre_copies: H[d] += 1      # the copied-centre copies are internal children (the loss l)
    mass = sum(d * c for d, c in H.items()) - loss
    return pr, H, loss, bad, mass

if __name__ == '__main__':
    pr, H, loss, bad, mass = run(); P = published()
    pub = P['remaining_internal_histogram']
    print('violations', dict(bad))
    print('sim internal (r>0)', {r: H[r] for r in range(1, 25) if H[r]})
    print('sim == published remaining_internal_histogram (r>0):', all(H[r] == pub[r] for r in range(1, 25)))
    diff = {r: (H[r], pub[r]) for r in range(1, 25) if H[r] != pub[r]}
    print('diffs', diff)
    print('copied-centre loss', loss, '(528)', ' sim internal mass + loss', mass + loss, ' Rh - sum d', pr['R'] * 24 - 20 * 4840 + loss)
    C = Counter()
    for r, c in H.items():
        if r: C[r] += 3 * c
    for hist in (pr['source'], pr['TH']):
        for r, c in hist.items():
            if r: C[r] += 3 * c
    C[2] += 2 * pr['v']; C[60] += 4840
    W = 2 * pr['v'] + pr['R']; rank = sum(r * c for r, c in C.items())
    print('W0', W, 'rank', rank, 'D', W * 72 - rank, ' children == published:', {str(k): x for k, x in C.items() if x} == {k: x for k, x in P['child_histogram'].items() if x})
    a, rf = acert(C, W); print('a_c* certified floor', a, float(a), 'float root', rf)
