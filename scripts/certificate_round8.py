"""Round-eight witness: kappa from the frozen bit and complex child-width histograms.

Bit side: round-seven witness 2 (joint frame compiler, lifted frames, late copies, deferred readouts, V leaves) with
every projector residual compiled as ONE child in opposite bank orders (PR #104's mechanism, notes and
independent/round8-oppbank/). The reversed family costs a linear adapter toll, so the recursion is stopped at atom
width e^theta and the usable bit saving is the mix
    a_b = (1 - theta) a* + theta a_old,   valid iff theta > a_b,
with a* the moment saving of the one-run histogram and a_old an already certified saving of the same witness
(round seven, round-six entrance corner).
Complex side: our two-stage complex word on eumemic's frozen addition DAG (PR #117) with signed dependence
reclamation (James Chang, PR #112) and completed-core sharing (Andrey Mas, PR #128), one fix-up child per group and
stream (independent/round8-coreshare/, gated by independent/round8-complex-gate/).

Savings are certified with rigorous rational bounds: L_w >= ln(m/w) (moment.ln_upper, rounded up), and either
(m/w)^a <= 1/(1 - a L_w) ('inv', the bound of rounds five to seven; the bit side, grid 10^-9, as gated) or
(m/w)^a <= exp_up(a L_w) with exp(x) <= 1 + x + x^2/2 + x^3/(6(1 - x/4)) ('exp'; the complex side, grid 10^-12).
lean/Round8.lean checks the same bounds in the kernel (momentOK, momentExpOK). Assembly: certificate_round3.evaluate(
a_b, a_c, 1/1000, 'crude', m_c, s_c) with s_c the complex rank sum; extra_rows() is the hook for further rows.

A histogram file is JSON (optionally gzipped): {side, source, m, W, s, hist: [[width, count], ...]}, optionally
bound ('exp' default, or 'inv' = (m/w)^a <= 1/(1 - a L_w)) and den (grid, default 10^12); the bit file also carries
stop = {theta, a_old} as 'p/q' strings. Optional 'claim' (bit: a_star, complex: a_c) is asserted.

Usage:
  python3 scripts/certificate_round8.py                       # check the frozen files in certificates/round8/
  python3 scripts/certificate_round8.py BIT.json CX.json      # evaluate other files (nothing written)
  python3 scripts/certificate_round8.py --freeze BIT.json CX.json
        # pin the computed savings as claims, write certificates/round8/{bit,cx}_hist.json.gz and
        # lean/round8-histograms.json"""
import gzip, json, math, os, sys
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
sys.path[:0] = [HERE, os.path.join(ROOT, 'independent', 'two-stage-bit')]
import moment
from certificate_round3 import evaluate

FROZEN = os.path.join(ROOT, 'certificates', 'round8')
BETA = Q(1, 1000)
DEN = 10**12              # grid of the certified savings
TARGET = Q(1, 2**13)


def extra_rows(r):
    """Assembly rows beyond certificate_round3.evaluate, as {name: margin}; every margin must be > 0, and kappa must
    stay below every entry of the optional 'kappa_caps' key. Empty today.
    Open item O8: certificate_round3 has no row-reserve row for the bit family. PR #104's assembly needs a stock t
    with m^t > 2 r_max^t (here r_max = m - 1 = 528, from the slots with r_u = 1). Our guard only needs every width
    below m, which holds; if a row-reserve lemma (or a row-stock cost such as the one PR #129 pays) turns out to be
    required, add its margin here, or cap r_u >= 2 in the histogram, and re-run --freeze."""
    return {}


def load(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        d = json.load(f)
    d['hist'] = {int(w): int(n) for w, n in d['hist']}
    return d


def ln_up_grid(x, den=10**40):
    """moment.ln_upper rounded UP to a 1/den grid: still a rigorous upper bound, with small denominators."""
    l = moment.ln_upper(x)
    return Q(-((-l.numerator * den) // l.denominator), den)


def exp_up(x):
    """exp(x) <= 1 + x + x^2/2 + x^3/(6(1 - x/4)) for 0 <= x < 4 (Taylor tail dominated by a geometric series)."""
    assert 0 <= x < 4
    return 1 + x + x * x / 2 + x ** 3 / (6 * (1 - x / 4))


def F_exp_upper(c, a, lu):
    """upper bound of (1/W) sum_w n_w (w/m)^(1-a), with (m/w)^a = exp(a ln(m/w)) <= exp_up(a L_w)."""
    m = c['m']
    return sum(Q(n * w, m) * exp_up(a * lu[w]) for w, n in c['hist'].items()) / c['W']


def F_inv_upper(c, a, lu):
    """moment.F_upper: (m/w)^a <= 1/(1 - a L_w), the bound of rounds five to seven (momentOK in Lean)."""
    return moment.F_upper(c, a, lu)


BOUNDS = {'exp': F_exp_upper, 'inv': F_inv_upper}


def certify_saving(c, den=DEN, bound='exp'):
    """Largest a on the 1/den grid with BOUNDS[bound](c, a) < 1, by bisection on the exact upper bound.
    Lean checks the same bound (momentExpOK for 'exp', momentOK for 'inv'). Returns (a, float root)."""
    F = BOUNDS[bound]
    lu = {w: ln_up_grid(Q(c['m'], w)) for w in c['hist']}
    lo, hi = 0.0, 0.05
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if moment.F_float(c, mid) < 1 else (lo, mid)
    ok = lambda n: n == 0 or F(c, Q(n, den), lu) < 1
    a, b = 0, int(lo * den) + 2                  # ok(a) holds; ok(b) fails (the float root bounds the exact one)
    assert not ok(b)
    while b - a > 1:
        mid = (a + b) // 2
        a, b = (mid, b) if ok(mid) else (a, mid)
    assert a > 0 and ok(a) and not ok(a + 1)
    return Q(a, den), lo


def check_hist(d):
    """rank sum exact, every width below m, W m > s (positive deficit)."""
    c = dict(m=d['m'], W=d['W'], s=d['s'], hist=d['hist'])
    assert sum(w * n for w, n in c['hist'].items()) == c['s'], 'rank sum'
    assert min(c['hist']) >= 1 and max(c['hist']) < c['m'], 'width range'
    assert c['W'] * c['m'] > c['s'], 'deficit'
    return c


def certify(bit, cx):
    cb, cc = check_hist(bit), check_hist(cx)
    astar, root_b = certify_saving(cb, int(bit.get('den', DEN)), bit.get('bound', 'exp'))
    a_c, root_c = certify_saving(cc, int(cx.get('den', DEN)), cx.get('bound', 'exp'))
    theta, a_old = Q(bit['stop']['theta']), Q(bit['stop']['a_old'])
    a_b = (1 - theta) * astar + theta * a_old
    assert 0 < a_old < astar and theta > a_b, 'stopping: theta must exceed a_b'
    for key, val in (('a_star', astar),):
        if 'claim' in bit: assert Q(bit['claim'][key]) == val, (key, val)
    if 'claim' in cx: assert Q(cx['claim']['a_c']) == a_c, ('a_c', a_c)
    k = evaluate(a_b, a_c, BETA, 'crude', m_c=cc['m'], s_c=cc['s'])
    rows = dict(extra_rows(dict(a_b=a_b, a_c=a_c, k=k, cb=cb, cc=cc)))
    caps = rows.pop('kappa_caps', [])
    bad = [name for name, g in rows.items() if not g > 0] + ['kappa cap %s' % c for c in caps if not k['kappa'] < c]
    if bad: k = dict(k, ok=False, bad=k['bad'] + bad)
    return dict(astar=astar, root_b=root_b, a_b=a_b, a_c=a_c, root_c=root_c, theta=theta, a_old=a_old, k=k, cb=cb, cc=cc)


def report(r, bit, cx):
    k = r['k']
    print('bit:     %s' % bit.get('source', ''))
    print('  m=%d W=%d s=%d widths=%d max=%d; a* = %s (%.7e, root %.7e)' % (
        r['cb']['m'], r['cb']['W'], r['cb']['s'], len(r['cb']['hist']), max(r['cb']['hist']), r['astar'], float(r['astar']), r['root_b']))
    print('  stopped a_b = (1-%s) a* + %s * %s = %s (%.7e)' % (r['theta'], r['theta'], r['a_old'], r['a_b'], float(r['a_b'])))
    print('complex: %s' % cx.get('source', ''))
    print('  m=%d W=%d s=%d widths=%d max=%d; a_c = %s (%.7e, root %.7e)' % (
        r['cc']['m'], r['cc']['W'], r['cc']['s'], len(r['cc']['hist']), max(r['cc']['hist']), r['a_c'], float(r['a_c']), r['root_c']))
    print('kappa = %s (%.7e = 2^%.4f) ok=%s binding=%s; above 2^-13: %s' % (
        k['kappa'], float(k['kappa']), math.log2(float(k['kappa'])), k['ok'], k['binding'], k['kappa'] > TARGET))


def fr(q): return [Q(q).numerator, Q(q).denominator]


def freeze(r, bit, cx):
    os.makedirs(FROZEN, exist_ok=True)
    bit = dict(bit, claim=dict(a_star=str(r['astar'])), hist=sorted(bit['hist'].items()))
    cx = dict(cx, claim=dict(a_c=str(r['a_c'])), hist=sorted(cx['hist'].items()))
    for name, d in (('bit_hist.json.gz', bit), ('cx_hist.json.gz', cx)):
        # mtime=0 keeps the gzip bytes reproducible
        with open(os.path.join(FROZEN, name), 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as f:
            f.write(json.dumps(d, sort_keys=True).encode())
    k = r['k']
    lean = dict(sides=[['bit', 'bit interchange, opposite bank orders (one-run histogram, moment at a*)'],
                       ['cx', 'complex interchange, completed-core sharing']],
                kappas=[['kappa', 'bit']],
                bit=dict(a=fr(r['astar']), bound=bit.get('bound', 'exp'), m=r['cb']['m'], W=r['cb']['W'], s=r['cb']['s'], hist=bit['hist'],
                         stop=dict(theta=fr(r['theta']), a_old=fr(r['a_old']), ab=fr(r['a_b']))),
                cx=dict(a=fr(r['a_c']), bound=cx.get('bound', 'exp'), m=r['cc']['m'], W=r['cc']['W'], s=r['cc']['s'], hist=cx['hist']),
                kappa=dict(beta=fr(BETA), eps=fr(k['eps']), x=fr(k['x']), c1=fr(k['c1']), kappa=fr(k['kappa']), mc=r['cc']['m']))
    json.dump(lean, open(os.path.join(ROOT, 'lean', 'round8-histograms.json'), 'w'))
    print('froze certificates/round8/{bit,cx}_hist.json.gz and lean/round8-histograms.json')


def main(argv):
    do_freeze = '--freeze' in argv
    args = [a for a in argv if a != '--freeze']
    if args:
        bpath, cpath = args
    else:
        bpath, cpath = os.path.join(FROZEN, 'bit_hist.json.gz'), os.path.join(FROZEN, 'cx_hist.json.gz')
    bit, cx = load(bpath), load(cpath)
    assert bit['side'] == 'bit' and cx['side'] == 'complex'
    r = certify(bit, cx)
    report(r, bit, cx)
    assert r['k']['ok'], r['k']['bad']
    if do_freeze:
        freeze(r, bit, cx)
    return r


if __name__ == '__main__':
    main(sys.argv[1:])
