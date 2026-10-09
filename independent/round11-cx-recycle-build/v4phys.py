"""Apply #168 v4's frozen physical layer (data: references/paired-cube/physical/{frames,pairs}.json at 4a3c769) to our
rebuild of its word (config v4f: same graph sha, frozen arcs, same ops), recount it with our code, then run our
levers on top: per-op pdescent with the pairs fixed, re-matching, and pair-aware reselection as a variant.
Writes plan2_v4f_p11<TAG>.pkl (sel, pairs, tau, frames) for export_word.py. Usage: python3 v4phys.py [resel]"""
import sys, os, json, pickle
from resel_br import *
import make_word, pdescent
p = 11; W = make_word.make(p, 'v4f'); F = gfirst_frames(W); sel = W['s']['sel']
fr = json.load(open('data/phys_frames_v4.json'))['frames']; pr_ = json.load(open('data/phys_pairs_v4.json'))['pairs']
frames = list(F['frames'])
for i, B in fr: frames[i] = tuple(B)
F = dict(F, frames=frames)
pairs = [(a, b) for a, b, _ in pr_]; tau = late_tau(W, F)
for a, b, t in pr_: tau[b] = F['pos'][t] if t is not None else F['cut']
C, Wv, inf = recount(W, F, pairs, tau, True); a, _ = acert(C, Wv, W['m'])
ref = json.load(open('pr168v4/paired-cube-physical-input.json'))
same = {k: n for k, n in C.items()} == {int(k): n for k, n in ref['child_histogram'].items() if n} and Wv == ref['W_per_vertex']
log('v4 frozen physical layer on our word: pairs %d bad %s W %s a_c %s = %.10e; equals #168 physical certificate: %s' % (len(pairs), inf['bad'], Wv, a, float(a), same))
best = (a, F, pairs)
Fc = F
for it in range(int(os.environ.get('PD_IT', 3))):
    Fc = pdescent.with_frames(Fc, pdescent.run(W, Fc, pairs, a=6.55e-4, log=log))
    C, Wv, inf = recount(W, Fc, pairs, tau, True); a, _ = acert(C, Wv, W['m'])
    log('round %d pdescent: bad %s W %s a_c %s = %.10e' % (it, inf['bad'], Wv, a, float(a)))
    if not inf['bad'] and a > best[0]: best = (a, Fc, pairs)
    p2 = match(candidates(W, Fc, tau, 5.6e-4, True, True))
    C, Wv, inf = recount(W, Fc, p2, tau, True); a, _ = acert(C, Wv, W['m'])
    log('round %d rematch: pairs %d bad %s W %s a_c %s = %.10e' % (it, len(p2), inf['bad'], Wv, a, float(a)))
    if not inf['bad'] and a > best[0]: best = (a, Fc, p2)
    pairs = p2
a, F, pairs = best
log('BEST v4f a_c %s = %.10e' % (a, float(a)))
pickle.dump(dict(sel=sel, pairs=pairs, tau=tau, frames=F['frames']), open('plan2_v4f_p11.pkl', 'wb'))
