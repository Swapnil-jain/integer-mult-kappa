"""Round-ten witness: kappa through PR #144's paired-cube / three-stage cover assembly (icekylinx), generalised to
fractional child weights (nullity-dependent sharing, NDS) and to per-side measured inventories (birth reuse, frame
descent, relaxed carrier arcs).

Own implementation; no outside code is imported or run. It reuses round nine's rigorous bounds, savings search,
stopped mix and 47-row assembly (scripts/certificate_round9.py) unchanged, and replaces only the cover profile and
the finite bridge:

  1. Cover profile, per side, built in layers; each layer is checked, and the last one is certified.
       base   either PR #144's parts (gauges, chains, data roles; rebuilt exactly as in round nine and compared with
              the saved histogram), or a plain measured histogram (child_histogram, W_per_vertex) from a gate replay;
       nds    optional, rebuilt here: a gauge of exit nullity u, radical r and sharing degree t
              (key "u,r,t") replaces its exterior child 3(h - u) of weight 1 by one child of rank m - t u at weight
              3/t, and W changes by 3/t - 1. The sharing degree must satisfy t = 3 or t (u + r) <= m;
       reuse  optional, measured (birth reuse: a recipient's slot is spliced onto its donor's chain, so its tail and
              at most one unit of W disappear). The gate recounts it from the actual role chains; here the deficit
              must be unchanged, W must drop by at most one per pair, and every child must stay below m.
     Every layer must keep the telescoping deficit D = W m - s = 2v - 3 loss. Counts and W may be fractions, or
     integers with a common 'scale' k (the saved counts and W are k times the per-vertex values).
  2. Savings, stopped mix, assembly: round nine's (rigorous ln / exp enclosures on a 2^-120 grid, the full bit
     fallback, a_bit = (1 - theta) a* + theta a_old, beta = 10^-6, eta = 10^-8, kappa on the 10^-10 grid).
  3. Finite bridge (complex): round nine's router charge, semantic guard and external row reserve, with the
     per-vertex stock scaled by L = the lcm of all weight denominators, so every count is an integer; the logical
     size uses the unscaled roles and the scaled stock, as the round-ten complex gate's certcheck.py does.

Regression: on PR #144's published inventories and on the frozen round-nine inventories, this script reproduces
their kappa (4609169/10^10 and round nine's frozen value), both through the parts path and through the plain path
(the parts stripped, only the saved histogram kept).

Usage:
  python3 scripts/certificate_round10.py                       # regressions + the frozen round-ten files
  python3 scripts/certificate_round10.py BIT.json CX.json      # evaluate other inventories (nothing written)
  python3 scripts/certificate_round10.py BIT.json CX.json --ladder=2  # the stopped-product ladder (OFF by default)
  python3 scripts/certificate_round10.py --freeze BIT.json CX.json
        # pin the computed savings as claims; write certificates/round10/{bit,cx}_inventory.json.gz and
        # lean/round10-histograms.json"""
import math, os, sys, json
from collections import Counter
from fractions import Fraction as Q

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import certificate_round9 as c9
from certificate_round9 import load, fr, _dump_gz, BAD, BETA, ETA, CX_GAP, KGRID, ROW_RESERVE, ROW_DEGREE

ROOT = c9.ROOT
FROZEN = os.path.join(ROOT, 'certificates', 'round10')
FROZEN9 = c9.FROZEN
GATE = Q(6, 10**4)              # the release gate: publish only at kappa >= 6.0e-4
PARTS = 'selected_rank_histogram'
PR144_KAPPA = Q(4609169, 10**10)


# ---------------------------------------------------------------- cover profile in layers
def qhist(hist):
    """{rank: count} with exact rational counts; zero ranks and zero counts dropped"""
    out = {}
    for r, n in hist.items():
        r, n = int(r), Q(n)
        assert n >= 0, ('negative multiplicity', r)
        if r and n:
            out[r] = out.get(r, 0) + n
    return dict(sorted(out.items()))


def _check_layer(name, H, W, m, v, loss):
    assert H and all(0 < r < m and n > 0 for r, n in H.items()), (name, 'every child narrower than m')
    s = sum(r * n for r, n in H.items())
    assert W * m - s == 2 * v - 3 * loss, (name, 'telescoping deficit 2v - 3 loss', W * m - s, 2 * v - 3 * loss)
    return s


def nds_layer(H, W, h, m, degrees):
    """rebuild NDS from the gauges' (u, r, t) types: the exterior child 3(h - u) at weight 1 becomes one child of rank
    m - t u at weight 3/t; W += 3/t - 1 per gauge"""
    H = Counter(H)
    for key, n in degrees.items():
        u, r, t = (int(x) for x in key.split(','))
        n = Q(n)
        assert 0 < u < h and r >= 0 and n > 0, ('gauge type', key)
        assert t == 3 or (t > 3 and t * (u + r) <= m), ('sharing degree', key)
        assert t * u <= m, ('tail rank', key)
        assert H[3 * (h - u)] >= n, ('gauge exteriors present', key)
        H[3 * (h - u)] -= n
        H[m - t * u] += n * Q(3, t)          # a zero tail rank is a full share and leaves no child
        W += n * (Q(3, t) - 1)
    return qhist(H), W


def profile10(inv, side):
    h, v, loss = int(inv['h']), int(inv['v']), int(inv['loss'])
    m = 3 * h
    assert int(inv.get('m', m)) == m, 'm = 3h'
    layers = []
    if PARTS in inv:
        saved = qhist(inv['child_histogram'])
        s_saved = sum(r * n for r, n in saved.items())
        # a gate inventory may omit the derived totals; they follow from the saved histogram, which the rebuild
        # below is compared with
        inv = dict(dict(rank_per_vertex=s_saved, deficit_per_vertex=inv['W_per_vertex'] * m - s_saved), **inv)
        p = c9.profile(inv, side)             # round nine's exact rebuild from the parts, compared with the save
        H, W = {r: Q(n) for r, n in p['hist'].items()}, Q(p['W'])
        layers.append('parts')
    else:
        # a gate may save NDS histograms as integers times a common scale k (counts and W both multiplied by k)
        k = int(inv.get('scale', 1))
        assert k > 0
        H, W = {r: n / k for r, n in qhist(inv['child_histogram']).items()}, Q(inv['W_per_vertex']) / k
        layers.append('measured')
    _check_layer(layers[-1], H, W, m, v, loss)
    if 'nds' in inv:
        nd = inv['nds']
        H, W = nds_layer(H, W, h, m, nd['degrees'])
        assert H == qhist(nd['child_histogram']) and W == Q(nd['W_per_vertex']), 'NDS rebuilt == saved'
        layers.append('nds')
        _check_layer('nds', H, W, m, v, loss)
    if 'reuse' in inv:
        ru = inv['reuse']
        k = int(ru.get('scale', 1))
        H2, W2 = {r: n / k for r, n in qhist(ru['child_histogram']).items()}, Q(ru['W_per_vertex']) / k
        pairs = int(ru['pairs'])
        assert 0 <= W - W2 <= pairs, ('reuse removes at most one unit of W per pair', W - W2, pairs)
        H, W = H2, W2
        layers.append('reuse')
        _check_layer('reuse', H, W, m, v, loss)
    s = sum(r * n for r, n in H.items())
    L = 1
    for x in list(H.values()) + [W]:
        L = math.lcm(L, x.denominator)
    return dict(m=m, h=h, v=v, R=int(inv['R']), loss=loss, W=W, s=s, D=W * m - s, hist=H, maxchild=max(H),
                edges=sum(H.values()), lcm=L, layers=layers)


# ---------------------------------------------------------------- sides (round nine's, on the layered profile)
def bit_side(inv):
    p = profile10(inv, 'bit')
    cert = inv['cert']
    astar, g, g_next, logs = c9.saving(p, cert, True)
    theta, a_old = Q(cert['theta']), Q(cert['a_old'])
    # stopped-product ladder (switch, default OFF: cert 'ladder' = k >= 2 rungs). Rung j uses rung j - 1 as its leaf,
    # a_j = (1 - theta) a* + theta a_{j-1}, so a* - a_k = theta^k (a* - a_old). PLAUSIBLE only (recursion Prop D):
    # each extra rung nests the stopped interchange once more and adds a row reserve that is NOT derived; see
    # LADDER_OPEN. k = 1 is the ordinary stop.
    k = int(cert.get('ladder', 1))
    assert k >= 1
    rungs = [a_old]
    for _ in range(k):
        rungs.append((1 - theta) * astar + theta * rungs[-1])
        assert rungs[-1] < theta < 1 - rungs[-1], 'every rung: a_j < theta < 1 - a_j'
    ab = rungs[-1]
    m = p['m']
    assert 0 < a_old < astar, 'stopped mix needs a_old < a*'
    assert ab < theta < 1 - ab, 'subordinate adapter and row tolls: a_bit < theta < 1 - a_bit'
    assert Q(2 * m ** 3, 2 ** 80) < BAD, 'fixed prime bad-class allowance'
    assert p['s'] + BAD * 32 * m * m * p['edges'] < p['W'] * m, 'contaminated rank mass'
    return dict(p=p, astar=astar, gap=g, gap_next=g_next, theta=theta, a_old=a_old, a_bit=ab, ladder=k, rungs=rungs,
                added=c9.fallback_upper(p, astar, logs), moment=c9.moment_upper(p, astar, logs))


def group_order(m):
    n = m // 2
    return 2 ** (m - 1 + (n - 1) ** 2) * math.prod(4 ** i - 1 for i in range(1, n))


def complex_side(inv):
    p = profile10(inv, 'complex')
    a, g, g_next, logs = c9.saving(p, inv['cert'], False)
    return dict(p=p, a_c=a, gap=g, gap_next=g_next, moment=c9.moment_upper(p, a, logs), V=group_order(p['m']))


def finite_bridge(c, inv):
    """round nine's finite bridge with the stock scaled by L (an integer number of copies of the cover per vertex)"""
    p, V = c['p'], c['V']
    m, h, v, R, L, r = p['m'], p['h'], p['v'], p['R'], p['lcm'], p['maxchild']
    W0, s0 = p['W'] * L, p['s'] * L
    assert W0.denominator == 1 and s0.denominator == 1
    W0, s0 = int(W0), int(s0)
    W, s, N = V * W0, V * s0, V * v * L
    mops = int(inv['total_M_operations'])
    local = 4 * (inv['c'] + v) + 10 * v + 4 * h * v + 4 * h * h + 8 * h + 8 + 2 * h
    local += 8 * R * v * (mops + 16) + 32 * v
    logical = 3 * V * local + 8 * W + 4 * N + 8 * m * R * V
    G = 64 * (m + 1) ** 3 * (logical + 1) * (W + 1) ** 2
    E = 64 * (W + m + G + 1) ** 3
    charge = 2 * G * W * W + 8 * s + 4 * W + 4 + 32 * m
    B = s + E
    C0 = 32 * m * B * B
    assert charge < E and 2 * B * (m - r) >= s + E and 2 * B + 18 < C0, 'finite semantic guard'
    dc, wc = c9.halving(m, r), W.bit_length()
    coefficient = dc * wc + ROW_RESERVE
    degree_gap = Q(ROW_DEGREE) - Q(51 * coefficient, 25)
    assert degree_gap > 0, 'external row reserve'
    return dict(local=local, logical=logical, G=G, E=E, charge=charge, literal_gap=E - charge, B=B, C0=C0,
                halving_degree=dc, wire_bits=wc, coefficient=coefficient, degree=ROW_DEGREE, degree_gap=degree_gap,
                W0=W0, s0=s0, lcm=L)


def inherited_row_reserve(b):
    """Proof-interface flag. ROW_RESERVE keeps PR #144's bit row constants (m = 69, largest child 66, W = 32408).
    They carry over unchanged when our bit side has the same m, a largest child no wider and a W no larger;
    otherwise the bit row term has to be re-derived before the reserve can be quoted."""
    p = b['p']
    return p['m'] == 69 and p['maxchild'] <= 66 and p['W'] <= 32408


# What the ladder still needs before it can be quoted (recursion agent, Prop D):
#   1. the row stock: every extra rung adds a row reserve to ROW_RESERVE in ROW_DEGREE - (51/25) coefficient > 0.
#      The provisional check below charges one full ROW_RESERVE per rung; that charge is a guess, not a derivation;
#   2. the adapter condition (e/w)^{a*} <= w^{1 - a_j} at every rung (needs theta >~ a*, checked here as a_j < theta);
#   3. the constant outer wrapper of each stopped interchange is an ordinary interchange of exponent 1 - a_j;
#   4. Lean: lean/gen.py checks each rung's stopped mix, but not items 1 to 3.
LADDER_OPEN = 'stopped-product ladder: per-rung row reserve, adapter condition and outer wrapper not derived'


def certify(bit_inv, cx_inv, kappa=None):
    b, c = bit_side(bit_inv), complex_side(cx_inv)
    bridge = finite_bridge(c, cx_inv)
    if b['ladder'] > 1:     # provisional: one full row reserve per rung (see LADDER_OPEN)
        coef = bridge['halving_degree'] * bridge['wire_bits'] + b['ladder'] * ROW_RESERVE
        bridge = dict(bridge, ladder_degree_gap=Q(ROW_DEGREE) - Q(51 * coef, 25))
    a = min(b['a_bit'], (1 - BETA) * c['a_c'] - CX_GAP)
    k = c9.assembly(a, c['a_c'], bridge, kappa)
    for side, d, key in (('bit', b, 'astar'), ('complex', c, 'a_c')):
        inv = bit_inv if side == 'bit' else cx_inv
        for src in ('claim', 'cert'):           # a gate may record its certified point in cert as well
            if key in inv.get(src, {}):
                assert Q(inv[src][key]) == d[key], (side, src, key, d[key])
    if 'claim' in bit_inv and 'a_bit' in bit_inv['claim']:
        assert Q(bit_inv['claim']['a_bit']) == b['a_bit']
    return dict(bit=b, cx=c, bridge=bridge, k=k, bit_binds=b['a_bit'] < (1 - BETA) * c['a_c'] - CX_GAP,
                inherited_reserve=inherited_row_reserve(b))


# ---------------------------------------------------------------- regressions
def plain(inv):
    """the same inventory with its parts stripped: only the saved histogram and W are kept"""
    keep = ('side', 'source', 'h', 'v', 'R', 'loss', 'm', 'W_per_vertex', 'child_histogram', 'cert', 'claim', 'c',
            'total_M_operations')
    return {k: inv[k] for k in keep if k in inv}


def regressions(verbose=True):
    bit, cx = c9.pr144_inputs()
    out = {}
    for path, f in (('parts', lambda d: d), ('plain', plain)):
        out['PR #144 kappa, ' + path] = certify(f(bit), f(cx))['k']['kappa'] == PR144_KAPPA
    b9 = load(os.path.join(FROZEN9, 'bit_inventory.json.gz'))
    x9 = load(os.path.join(FROZEN9, 'cx_inventory.json.gz'))
    k9 = c9.certify(b9, x9)['k']['kappa']
    for path, f in (('parts', lambda d: d), ('plain', plain)):
        out['round 9 kappa, ' + path] = certify(f(b9), f(x9))['k']['kappa'] == k9
    if verbose:
        print('regressions (PR #144 kappa %s, round 9 kappa %s):' % (PR144_KAPPA, k9))
        for name, ok in out.items():
            print('  %-50s %s' % (name, 'PASS' if ok else 'FAIL'))
    assert all(out.values()), [n for n, ok in out.items() if not ok]
    return out


# ---------------------------------------------------------------- report and freeze
def report(r, bit_inv, cx_inv):
    b, c, k = r['bit'], r['cx'], r['k']
    for name, d, inv in (('bit', b, bit_inv), ('complex', c, cx_inv)):
        p = d['p']
        print('%s: %s' % (name, inv.get('source', '')))
        print('  layers %s; h=%d v=%d R=%d loss=%d; m=%d W=%s s=%s D=%s; widths=%d max=%d lcm=%d' % (
            '+'.join(p['layers']), p['h'], p['v'], p['R'], p['loss'], p['m'], p['W'], p['s'], p['D'], len(p['hist']),
            p['maxchild'], p['lcm']))
    print('  bit a* = %s (%.12e), gap %.3e, next point %s' % (b['astar'], float(b['astar']), float(b['gap']),
          'pinned' if b['gap_next'] is None else 'fails (gap %.3e)' % float(b['gap_next'])))
    print('  stopped a_bit = (1-%s) a* + %s * %s = %s (%.12e)' % (b['theta'], b['theta'], b['a_old'], b['a_bit'], float(b['a_bit'])))
    print('  complex a_c = %s (%.12e), gap %.3e, next point %s' % (c['a_c'], float(c['a_c']), float(c['gap']),
          'pinned' if c['gap_next'] is None else 'fails (gap %.3e)' % float(c['gap_next'])))
    print('assembly: a = %s (%s binds); minimum margin eps q = %.11e' % (k['a'], 'bit' if r['bit_binds'] else 'complex', float(k['minimum'])))
    print('kappa = %s (%.7e = 2^%.4f), ok=%s, 47 constraints and 7 margins; >= 6.0e-4: %s' % (
        k['kappa'], float(k['kappa']), math.log2(float(k['kappa'])), k['ok'], k['kappa'] >= GATE))
    print('row reserve inherited from PR #144 applies to this bit side: %s' % r['inherited_reserve'])
    if r['bit']['ladder'] > 1:
        print('LADDER ON (k = %d rungs, PLAUSIBLE): %s; provisional row gap with one reserve per rung %.1f' % (
            r['bit']['ladder'], LADDER_OPEN, float(r['bridge']['ladder_degree_gap'])))
    if k['bad']:
        print('  failing rows:', k['bad'])


def scaled(p):
    """the histogram, W, s and edges times the lcm L, all integers (the moment inequality is homogeneous)"""
    L = p['lcm']
    hist = [(r, int(n * L)) for r, n in sorted(p['hist'].items())]
    assert all(n * L == x for (r, n), (_, x) in zip(sorted(p['hist'].items()), hist))
    return dict(m=p['m'], W=int(p['W'] * L), s=int(p['s'] * L), edges=int(p['edges'] * L), hist=hist, lcm=L)


def lean_input(r, bit_inv, cx_inv):
    """lean/round10-histograms.json: both moments on the L-scaled histograms (bit with the fallback), the stopped mix,
    the bridge inputs and the assembly"""
    b, c, k, br = r['bit'], r['cx'], r['k'], r['bridge']
    return dict(kind='cover10',
                bit=dict(scaled(b['p']), a=fr(b['astar']), label=bit_inv.get('source', 'bit'), fallback=True,
                         stop=dict(theta=fr(b['theta']), a_old=fr(b['rungs'][-2]), ab=fr(b['a_bit']),   # the last rung's leaf
                                   rungs=[fr(x) for x in b['rungs']])),
                cx=dict(scaled(c['p']), a=fr(c['a_c']), label=cx_inv.get('source', 'complex'), fallback=False,
                        h=c['p']['h'], v=c['p']['v'], R=c['p']['R'], cAdd=cx_inv['c'],
                        mops=int(cx_inv['total_M_operations']), maxchild=c['p']['maxchild'],
                        halving=br['halving_degree'], wire_bits=br['wire_bits']),
                bad=fr(BAD), beta=fr(BETA), eta=fr(ETA), cx_gap=fr(CX_GAP), row_reserve=ROW_RESERVE,
                row_degree=ROW_DEGREE, kappa=fr(k['kappa']))


def freeze(r, bit_inv, cx_inv):
    os.makedirs(FROZEN, exist_ok=True)
    bit_inv = dict(bit_inv, claim=dict(astar=str(r['bit']['astar']), a_bit=str(r['bit']['a_bit'])))
    cx_inv = dict(cx_inv, claim=dict(a_c=str(r['cx']['a_c'])))
    _dump_gz(os.path.join(FROZEN, 'bit_inventory.json.gz'), bit_inv)
    _dump_gz(os.path.join(FROZEN, 'cx_inventory.json.gz'), cx_inv)
    with open(os.path.join(ROOT, 'lean', 'round10-histograms.json'), 'w') as f:
        json.dump(lean_input(r, bit_inv, cx_inv), f)
    with open(os.path.join(FROZEN, 'summary.json'), 'w') as f:
        k = r['k']
        json.dump(dict(kappa=str(k['kappa']), kappa_float=float(k['kappa']), log2=math.log2(float(k['kappa'])),
                       a_bit=str(r['bit']['a_bit']), astar=str(r['bit']['astar']), a_c=str(r['cx']['a_c']),
                       bit_binds=r['bit_binds'], inherited_reserve=r['inherited_reserve'],
                       bit_layers=r['bit']['p']['layers'], cx_layers=r['cx']['p']['layers'],
                       theta=str(r['bit']['theta']), ladder=r['bit']['ladder'],
                       ladder_open=LADDER_OPEN if r['bit']['ladder'] > 1 else None), f, indent=1, sort_keys=True)
    print('froze certificates/round10/{bit,cx}_inventory.json.gz, summary.json and lean/round10-histograms.json')


def main(argv):
    do_freeze = '--freeze' in argv
    ladder = [int(a.split('=', 1)[1]) for a in argv if a.startswith('--ladder=')]
    args = [a for a in argv if not a.startswith('--')]
    c9.reproduce_pr144()
    regressions()
    if args:
        bpath, cpath = args
    else:
        bpath, cpath = os.path.join(FROZEN, 'bit_inventory.json.gz'), os.path.join(FROZEN, 'cx_inventory.json.gz')
    bit, cx = load(bpath), load(cpath)
    assert bit['side'] == 'bit' and cx['side'] == 'complex'
    if ladder:
        bit = dict(bit, cert=dict(bit['cert'], ladder=ladder[0]))
    if do_freeze and int(bit['cert'].get('ladder', 1)) > 1 and '--allow-ladder' not in argv:
        sys.exit('the ladder is PLAUSIBLE only (%s); pass --allow-ladder to freeze it anyway' % LADDER_OPEN)
    r = certify(bit, cx)
    report(r, bit, cx)
    assert r['k']['ok'], r['k']['bad']
    if do_freeze:
        freeze(r, bit, cx)
    return r


if __name__ == '__main__':
    main(sys.argv[1:])
