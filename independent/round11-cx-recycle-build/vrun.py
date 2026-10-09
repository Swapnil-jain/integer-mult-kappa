"""Validation at small p: our word (descent) + pair-aware gauges + birth-read pairs; recount legality, exact replay,
controls (scalar and frame)."""
import sys, os
from resel_br import *
from replay import validate
import make_word
p = int(sys.argv[1]); name = sys.argv[2] if len(sys.argv) > 2 else 'default'; a0 = 5.6e-4
W = make_word.make(p, name, cache=False); F0 = gfirst_frames(W)
log('p=%d %s: R %d ops %d cut %d gauges(#144) %d' % (p, name, W['pr']['R'], len(W['s']['ops']), F0['cut'], len(W['s']['sel'])))
ok0 = validate(W, F0, [], late_tau(W, F0), log)
sel, res, info = reselect(W, F0, a0, set(), True)
W2 = with_sel(W, sel); F = gfirst_frames(W2); tau = late_tau(W2, F)
pairs = match(candidates(W2, F, tau, a0, True, True))
C, Wv, inf = recount(W2, F, pairs, tau, True)
log('gauges %d pairs %d chained %d; recount bad %s; D %s' % (len(sel), len(pairs), len({a for a, _ in pairs} & {b for _, b in pairs}), inf['bad'], D_of(C, Wv, W['m'])))
ok, cok, r, c = validate(W2, F, pairs, tau, log)
# frame control: re-pair one recipient with a donor whose last frame is NOT inside its sigma (timely)
rec = {b for _, b in pairs}; don = {a for a, _ in pairs}
bad = None
for z in sel:
    B = z['role']
    if B not in rec: continue
    for A in range(W['pr']['R']):
        if A in don or A in rec or A in F['rootframe'] or not F['rops'][A]: continue
        if F['pos'][F['rops'][A][-1]] < tau[B] and not contained(F['frames'][F['rops'][A][-1]], F['gsig'][B]):
            bad = (A, B); break
    if bad: break
if bad:
    pr2 = [(a, b) for a, b in pairs if b != bad[1]] + [bad]
    _, _, inf2 = recount(W2, F, pr2, tau, True)
    log('frame control (donor last frame outside sigma): recount bad %s' % inf2['bad'])
# dropped birth: recount with the recipient removed from pairs but its slot also removed is not expressible; the
# scalar 'omit recipient read' control above is the dropped birth read.
log('RESULT p=%d: replay %s, scalar controls rejected %s' % (p, 'PASS' if ok and ok0[0] else 'FAIL', cok))
