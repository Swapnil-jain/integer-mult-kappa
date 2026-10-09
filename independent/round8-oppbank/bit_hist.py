"""Round-eight bit histogram: round-seven witness 2 with every projector residual compiled as ONE child in opposite
bank orders (PR #104's mechanism: conjugating by diag(I, C_n) turns the generic anti-diagonal Bruhat cell of each
partial swap into one contiguous run, so a residual of rank r is one child of width r).

Rebuilt from the frozen round-seven schedule (certificates/round7/, scripts/certificate_round7.build); per class:
  slot aux edge, rank m - r_u        [m - 2 r_u] + inner(r_u)       ->  [m - r_u]
  chain / y_T / X_S / centre step r  inner(r)                         ->  [r]
  data entrance, rank m - 2h + 1     [m - 4h + 2] + corner(2h - 1)   ->  [m - 2h + 1]
  copy correction                    [1]                              ->  [1]
The rank multiset per edge is unchanged, so W and s equal round seven's (asserted). The certified round-seven saving
of the same witness (round-six entrance corner) is reproduced first: it is a_old in the stopped mix.
independent/round8-oppbank/oppbank_check.py checks, on the real edges, that one common rational basis makes every
residual a single reversed run.

Usage: python3 independent/round8-oppbank/bit_hist.py [OUT.json]       # default: compare with the frozen file"""
import gzip, json, os, sys
from collections import Counter
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path[:0] = [os.path.join(ROOT, 'scripts'), os.path.join(ROOT, 'independent', 'two-stage-bit'),
                os.path.join(ROOT, 'independent', 'deferred-readout'), os.path.join(ROOT, 'independent', 'joint-frame-stack')]
import certificate_round7 as c7

THETA = Q(1, 1000)


def onerun_hist(S, rk, ylev, xd):
    h, v, R = S.h, S.v, S.R; m = h * h; N = v * v; H = Counter()
    corners = Counter(h - S.f[s] for s in range(R))
    for _stage in range(2):
        for r, c in corners.items(): H[m - r] += v * c          # aux edge with its corner: one child m - r_u
        for r, c in rk.items(): H[r] += v * c                   # side chain steps
        H[h - 1] += v * h                                       # copied centres U_c -> 0
        for t in range(v):                                      # y_T steps between deferred readout levels
            ds = sorted(set([0] + ylev[t] + [h - 1]))
            for a, b in zip(ds, ds[1:]): H[b - a] += v
        for xs in xd:                                           # X_S chains
            for x in xs: H[x] += v
    H[m - 2 * h + 1] += 2 * N                                   # data entrance, one run
    H[1] += N                                                   # copy corrections
    W = 2 * N + 2 * v * R; L = 2 * v * h * (h - 1); s = W * m - N + L
    H = {w: c for w, c in sorted(H.items()) if c}
    assert sum(w * c for w, c in H.items()) == s, 'rank sum'
    return dict(m=m, W=W, s=s, hist=H)


def main(argv):
    wit = c7.HEADLINE
    S, rk, ylev, xd = c7.build(wit)
    a_old, base = c7.saving(S, rk, ylev, xd, c7.dr.ROUND6_CORNER(S.h))
    assert a_old == wit['plain'][0], a_old
    c = onerun_hist(S, rk, ylev, xd)
    assert (c['W'], c['s']) == (base['W'], base['s']), 'W and s must equal round seven'
    out = dict(side='bit', source='opposite-bank one-run compile of round-7.1 witness 2 (h=%d)' % S.h,
               m=c['m'], W=c['W'], s=c['s'], hist=sorted(c['hist'].items()),
               bound='inv', den=10**9,                    # the bound and grid of the gated a_b
               stop=dict(theta=str(THETA), a_old=str(a_old)))
    print('h=%d R=%d W=%d s=%d widths=%d max=%d a_old=%s' % (S.h, S.R, c['W'], c['s'], len(c['hist']), max(c['hist']), a_old))
    if argv:
        json.dump(out, open(argv[0], 'w')); print('wrote', argv[0]); return
    with gzip.open(os.path.join(ROOT, 'certificates', 'round8', 'bit_hist.json.gz'), 'rt') as f:
        fz = json.load(f)
    same = all(fz[k] == out[k] for k in ('m', 'W', 's', 'stop', 'bound', 'den')) and [list(p) for p in out['hist']] == fz['hist']
    print('matches certificates/round8/bit_hist.json.gz:', same)
    assert same


if __name__ == '__main__':
    main(sys.argv[1:])
