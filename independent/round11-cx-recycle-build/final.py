"""Final certified number + real-size replay for one configuration.
Usage: python3 final.py p name [replay]   (name: a make_word.CFG key, or 'best' for cube-stack's best_p pickle)"""
import sys, os, json, math, pickle
from fractions import Fraction as Q
from resel_br import *
import make_word
p = int(sys.argv[1]); name = sys.argv[2]; do_replay = len(sys.argv) > 3; a0 = float(os.environ.get('A0', 5.6e-4))
W = load_word(p) if name == 'best' else make_word.make(p, name)
F0 = gfirst_frames(W)
sel, res, info = reselect(W, F0, a0, set(), True)
W2 = with_sel(W, sel); F = gfirst_frames(W2); tau = late_tau(W2, F)
pairs = match(candidates(W2, F, tau, a0, True, True))
out = {}
for nds in (True, False):
    C, Wv, inf = recount(W2, F, pairs, tau, nds)
    a, rf = acert(C, Wv, W['m'])
    out['nds' if nds else 'plain'] = dict(a=str(a), af=float(a), W=str(Wv), D=str(D_of(C, Wv, W['m'])), bad=inf['bad'],
                                          deg=dict(inf['deg']), C={str(r): str(n) for r, n in sorted(C.items())})
    log('p=%d %s nds=%s: gauges %d pairs %d chained %d unpaired gauges %d bad %s D %s W %s a_c %s = %.10e' % (p, name, nds,
        len(sel), len(pairs), len({x for x, _ in pairs} & {y for _, y in pairs}), len(sel) - len(pairs), inf['bad'],
        D_of(C, Wv, W['m']), Wv, a, float(a)))
# finite semantic guard (cube-stack/run_p.py formulas, #144's assembly note) with this word's counts
C, Wv, inf = recount(W2, F, pairs, tau, True); h = W['h']; m = W['m']; pr = W['pr']
k = m // 2; V = 2 ** (2 * k - 1)
for i in range(1, k): V *= 2 ** (2 * i) - 1
V *= 2 ** ((k - 1) ** 2)
lcm = math.lcm(3, Wv.denominator, *[Q(n).denominator for n in C.values()])
wc = int(Wv * lcm); sc = int(sum(Q(r) * n for r, n in C.items()) * lcm); vv = pr['v'] * lcm
Wf, sf, Nf = V * wc, V * sc, V * vv; M = len(W['s']['ops']); Rp = pr['R'] - len(pairs)
L = 4 * (pr['c'] + pr['v']) + 10 * pr['v'] + 4 * h * pr['v'] + 4 * h * h + 8 * h + 8 + 2 * h + 8 * pr['R'] * pr['v'] * (M + 16) + 32 * pr['v']
K = 3 * V * L + 8 * Wf + 4 * Nf + 8 * m * pr['R'] * V; Gg = 64 * (m + 1) ** 3 * (K + 1) * (Wf + 1) ** 2
E = 64 * (Wf + m + Gg + 1) ** 3; B = sf + E; C0 = 32 * m * B * B; mx = max(r for r, n in C.items() if n)
guard = (2 * Gg * Wf * Wf + 8 * sf + 4 * Wf + 4 + 32 * m < E, 2 * B * (m - mx) >= sf + E, 2 * B + 18 < C0)
log('guard (m=%d, max child %d, lcm %d): %s' % (m, mx, lcm, guard)); out['guard'] = list(guard)
out.update(p=p, name=name, gauges=len(sel), pairs=len(pairs), R=pr['R'], physical_slots=Rp, ops=M)
if do_replay:
    from replay import validate
    ok, cok, r, c = validate(W2, F, pairs, tau, log)
    out['replay'] = dict(pass_=ok, controls_rejected=cok, runs={k_: list(v_) for k_, v_ in r.items()}, controls={k_: list(v_) for k_, v_ in c.items()})
json.dump(out, open(os.path.join(HERE, 'final_%s_p%d%s.json' % (name, p, '' if do_replay else '_norep')), 'w'), indent=1)
pickle.dump(dict(sel=sel, pairs=pairs, tau=tau), open(os.path.join(HERE, 'plan_%s_p%d.pkl' % (name, p)), 'wb'))
