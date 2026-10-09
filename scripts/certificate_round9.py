"""Round-nine witness: kappa through IceKylin's paired-cube / three-stage cover assembly (icekylinx, PR #144, merged
upstream through PR #149), with our witness 2 on the bit side.

This is our own implementation of PR #144's exact assembly, written from its published description and certificate
files. None of its code is imported or run. It reads frozen inventories (JSON, optionally gzipped) and checks:

  1. Cover moments, per side. The child histogram of the shared-core cover is rebuilt from the inventory's parts.
     There are m = 3h coordinates and W = 2v + R roles per vertex, and the parts contribute:
       - a retained gauge of rank r: one exterior child of width 3r;
       - every chain histogram (bit: auxiliary, source, target, copied centre; complex: remaining internal, source,
         target): three children per entry, one per stage;
       - each of the 2v data roles: one finishing child of width 2.
     The rebuilt histogram must equal the saved one, and the deficit must telescope: W m - s = 2v - 3 loss.
  2. Savings, with PR #144's rigorous rational bounds:
       - L_w >= ln(m/w), from ln x = k ln 2 + 2 atanh((x-1)/(x+1)) (32 terms plus a geometric tail), rounded up to
         a 2^-120 grid;
       - exp u <= 1 + u + u^2/(2(1 - u/3)) for 0 <= u < 3.
     The moment is M(a) = sum_w n_w (w/(W m)) exp_up(a L_w), and a saving a is certified iff M(a) < 1.
     The bit side also pays the full rare-class fallback: a fraction 10^-16 of edges each spawn 32 m^2 width-1
     children, priced at width 1 (the "added bad moment"), so the bit condition is 1 - M(a) - added > 0.
     A saving is either pinned in the inventory (PR #144 pins both of its own) or searched: the largest point on a
     1/den grid, with the next grid point checked to fail.
  3. Stopped mix (bit). a_bit = (1 - theta) a* + theta a_old, valid for a_old < a* and a_bit < theta < 1 - a_bit.
     PR #144 uses theta = 1/1000; PR #148 proposes 1/2000. theta is a parameter of the bit inventory.
  4. Finite bridge (complex). The complex cover runs over a group of order
     V = 2^(m-1+(n-1)^2) prod_{i<n} (4^i - 1), n = m/2. The literal router charge must stay below
     E = 64(W+m+G+1)^3, and the semantic guard and the external row reserve must hold. The row reserve keeps PR
     #144's conservative bit constants (9909 + 252); see ROW_RESERVE below.
  5. Assembly: all 47 strict constraints and the 7 cost margins of PR #144's semantic assembly, at
     a = min(a_bit, (1 - beta) a_c - 10^-10), with beta = 10^-6 and eta = 10^-8. The minimum margin is
     eps q = (1 - eta) q / (1 + c + q), where q = a(1 - 2 eta) and c = q(1 + eta). kappa is that minimum rounded down
     to the 10^-10 grid, as PR #144 states its own kappa, and every constraint is rechecked at that kappa.

Inventory schema: PR #144's input files (certificates/round9/pr144/paired-cube-{bit,complex}-input.json), plus:
  side: 'bit' or 'complex';
  cert: {saving: 'p/q'} (pinned) or {den: N} (searched), and for the bit {theta, a_old};
  claim: {...} (optional, asserted).

Usage:
  python3 scripts/certificate_round9.py                      # PR #144 reproduction + the frozen round-nine files
  python3 scripts/certificate_round9.py --pr144              # PR #144 reproduction only
  python3 scripts/certificate_round9.py BIT.json CX.json     # evaluate other inventories (nothing written)
  python3 scripts/certificate_round9.py --freeze BIT.json CX.json
        # pin the computed savings as claims; write certificates/round9/{bit,cx}_inventory.json.gz and
        # lean/round9-histograms.json"""
import gzip, json, math, os, sys
from collections import Counter
from fractions import Fraction as Q

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
FROZEN = os.path.join(ROOT, 'certificates', 'round9')
PR144 = os.path.join(FROZEN, 'pr144')

ETA = Q(1, 10**8)
BETA = Q(1, 10**6)              # PR #144's phase stop ("beta" in its assembly)
CX_GAP = Q(1, 10**10)           # the bit saving fed to the assembly stays this far below (1 - beta) a_c
BAD = Q(1, 10**16)              # rare-class fraction of edges that take the fallback
KGRID = 10**10                  # kappa is stated rounded down to this grid
GRID = 1 << 120                 # ln and exp upper bounds are rounded up to this grid
TARGET = Q(1, 2**12)
# External row reserve of the finite bridge. PR #144 keeps two inherited, conservative bit row constants unchanged
# (9909 = 367 * 27, the partial-gauge round's coarse bit family; 252 = 9 * 28, the ordinary leaf), and
# charges the complex halving degree times its wire bits on top. Neither constant is a term of the current bit cover:
# its stock W0 q^(m^2 (w-1)) grows with the atom width, so it is borrowed and restored internally
# (three-stage-cover-rows.tex) at any fixed m, and is never part of this external reserve. This is a proof
# interface: the constants are inherited, not re-derived; see inherited_row_reserve in certificate_round10.py.
ROW_RESERVE = 9909 + 252
ROW_DEGREE = 70000


# ---------------------------------------------------------------- rigorous bounds (PR #144's, reimplemented)
def up(x):
    """round up to the 2^-120 grid"""
    return Q(-((-x.numerator * GRID) // x.denominator), GRID)


def _atanh2(y):
    """2 atanh((y-1)/(y+1)) = ln y from above: 32 series terms plus the geometric tail"""
    z = (y - 1) / (y + 1)
    return 2 * sum((z ** (2 * j + 1) / (2 * j + 1) for j in range(32)), Q(0)) + 2 * z ** 65 / (65 * (1 - z * z))


def log_upper(x):
    """an upper bound of ln x for x >= 1: halve to [1, 2), then ln x = k ln 2 + ln y"""
    x, k = Q(x), 0
    while x >= 2:
        x /= 2
        k += 1
    return up(k * _atanh2(Q(2)) + _atanh2(x))


def exp_upper(u):
    """exp u <= 1 + u + u^2/(2(1 - u/3)) for 0 <= u < 3"""
    assert 0 <= u < 3, 'exponential enclosure'
    return up(1 + u + u * u / (2 * (1 - u / 3)))


# ---------------------------------------------------------------- inventories and cover profiles
def load(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        return json.load(f)


def clean(hist):
    return {int(r): int(n) for r, n in sorted(hist.items(), key=lambda t: int(t[0])) if int(r) and n}


def profile(inv, side):
    """the shared-core cover's child histogram per vertex, rebuilt from the inventory's parts and checked against
    the saved one (PR #144, shared_profile; generalised from its two fixed dimensions to any inventory)"""
    h, v, R, loss = (int(inv[k]) for k in ('h', 'v', 'R', 'loss'))
    m, W, H = 3 * h, 2 * v + R, Counter()
    sel = clean(inv['selected_rank_histogram'])
    assert sum(sel.values()) == inv['selected_roles'], 'selected gauge count'
    for r, n in sel.items():
        assert 0 < r < h and n > 0, 'proper local gauges'
        H[3 * r] += n
    if side == 'complex':
        if 'matched' in inv:
            assert R == inv['c'] + inv['q'] - inv['matched'], 'compatible carrier roles'
        for r, n in enumerate(inv['remaining_internal_histogram']):
            H[r] += 3 * n
        parts = ('source_data_histogram', 'target_data_histogram')
    else:
        parts = ('auxiliary_histogram', 'source_data_histogram', 'target_data_histogram', 'copied_center_histogram')
    for name in parts:
        for r, n in inv[name].items():
            H[int(r)] += 3 * n
    H[2] += 2 * v
    H = clean(H)
    assert all(0 < r < m and n > 0 for r, n in H.items()), 'every child narrower than m'
    s = sum(r * n for r, n in H.items())
    assert H == clean(inv['child_histogram']), 'saved child histogram'
    assert (m, W, s, W * m - s) == (inv['m'], inv['W_per_vertex'], inv['rank_per_vertex'], inv['deficit_per_vertex']), \
        'saved dimensions and rank'
    assert W * m - s == 2 * v - 3 * loss, 'telescoping deficit 2v - 3 loss'
    return dict(m=m, h=h, v=v, R=R, loss=loss, W=W, s=s, D=W * m - s, hist=H, maxchild=max(H), edges=sum(H.values()))


# ---------------------------------------------------------------- moments and savings
def moment_upper(p, a, logs):
    m, W = p['m'], p['W']
    return sum(Q(n * r, W * m) * exp_upper(a * logs[r]) for r, n in p['hist'].items())


def fallback_upper(p, a, logs):
    """the full fallback: BAD of the edges each spawn 32 m^2 children of width 1, nothing displaced"""
    m, W = p['m'], p['W']
    return BAD * Q(32 * m * m * p['edges'], W * m) * exp_upper(a * logs[1])


def gap(p, a, logs, fallback):
    return 1 - moment_upper(p, a, logs) - (fallback_upper(p, a, logs) if fallback else 0)


def _float_root(p, fallback):
    m, W = p['m'], p['W']
    fb = float(BAD) * 32 * m * m * p['edges'] if fallback else 0.0
    F = lambda a: (math.fsum(n * (r / m) ** (1 - a) for r, n in p['hist'].items()) + fb * (1 / m) ** (1 - a)) / W
    lo, hi = 0.0, 0.05
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if F(mid) < 1 else (lo, mid)
    return lo


def saving(p, cert, fallback):
    """pinned: check gap > 0 at cert['saving']. searched: the largest a on the 1/den grid with gap > 0, and the next
    point fails. Returns (a, exact gap at a, next point's gap or None)."""
    logs = {r: log_upper(Q(p['m'], r)) for r in set(p['hist']) | {1}}
    if 'saving' in cert:
        a = Q(cert['saving'])
        g = gap(p, a, logs, fallback)
        assert g > 0, 'strict moment at the pinned saving'
        return a, g, None, logs
    den = int(cert['den'])
    ok = lambda k: gap(p, Q(k, den), logs, fallback) > 0
    lo, hi = 0, int(_float_root(p, fallback) * den) + 2   # the float root brackets the exact one
    assert not ok(hi)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    assert lo > 0
    a = Q(lo, den)
    return a, gap(p, a, logs, fallback), gap(p, a + Q(1, den), logs, fallback), logs


def bit_side(inv):
    p = profile(inv, 'bit')
    cert = inv['cert']
    astar, g, g_next, logs = saving(p, cert, True)
    theta, a_old = Q(cert['theta']), Q(cert['a_old'])
    ab = (1 - theta) * astar + theta * a_old
    m = p['m']
    assert 0 < a_old < astar, 'stopped mix needs a_old < a*'
    assert ab < theta < 1 - ab, 'subordinate adapter and row tolls: a_bit < theta < 1 - a_bit'
    assert Q(2 * m ** 3, 2 ** 80) < BAD, 'fixed prime bad-class allowance'
    assert p['s'] + BAD * 32 * m * m * p['edges'] < p['W'] * m, 'contaminated rank mass'
    return dict(p=p, astar=astar, gap=g, gap_next=g_next, theta=theta, a_old=a_old, a_bit=ab,
                added=fallback_upper(p, astar, logs), moment=moment_upper(p, astar, logs))


def complex_side(inv):
    p = profile(inv, 'complex')
    a, g, g_next, logs = saving(p, inv['cert'], False)
    m = p['m']
    n = m // 2
    V = 2 ** (m - 1 + (n - 1) ** 2) * math.prod(4 ** i - 1 for i in range(1, n))
    return dict(p=p, a_c=a, gap=g, gap_next=g_next, moment=moment_upper(p, a, logs), V=V, N=V * p['v'], W=V * p['W'],
                s=V * p['s'], deficit=V * p['D'])


# ---------------------------------------------------------------- finite bridge and assembly
def halving(m, r):
    """least n with m^n > 2 r^n"""
    assert 0 < r < m
    n = 1
    while m ** n <= 2 * r ** n:
        n += 1
    return n


def finite_bridge(c, inv):
    """PR #144's finite router and semantic guard for the complex cover, and the external row reserve"""
    m, W, s, N, V, r = c['p']['m'], c['W'], c['s'], c['N'], c['V'], c['p']['maxchild']
    h, v, R = c['p']['h'], c['p']['v'], c['p']['R']
    local = 4 * (inv['c'] + v) + 10 * v + 4 * h * v + 4 * h * h + 8 * h + 8 + 2 * h
    local += 8 * R * v * (inv['total_M_operations'] + 16) + 32 * v
    logical = 3 * V * local + 8 * W + 4 * N + 8 * m * R * V
    G = 64 * (m + 1) ** 3 * (logical + 1) * (W + 1) ** 2
    E = 64 * (W + m + G + 1) ** 3
    charge = 2 * G * W * W + 8 * s + 4 * W + 4 + 32 * m
    B = s + E
    C0 = 32 * m * B * B
    assert charge < E and 2 * B * (m - r) >= s + E and 2 * B + 18 < C0, 'finite semantic guard'
    dc, wc = halving(m, r), W.bit_length()
    coefficient = dc * wc + ROW_RESERVE
    degree_gap = Q(ROW_DEGREE) - Q(51 * coefficient, 25)
    assert degree_gap > 0, 'external row reserve'
    return dict(local=local, logical=logical, G=G, E=E, charge=charge, literal_gap=E - charge, B=B, C0=C0,
                halving_degree=dc, wire_bits=wc, coefficient=coefficient, degree=ROW_DEGREE, degree_gap=degree_gap)


def assembly(a, b, bridge, kappa=None, eta=ETA, beta=BETA):
    """PR #144's semantic assembly: 7 cost margins, 47 strict constraints (40 rows + every margin above kappa).
    kappa=None: the minimum margin rounded down to the 1/KGRID grid."""
    tau, sigma = 1 - a, 1 - b
    q = a * (1 - 2 * eta)
    lp = 1 - q
    lam = (tau + lp) / 2
    c = q * (1 + eta)
    eps = (1 - eta) / (1 + c + q)
    minimum = eps * q
    r = (minimum + 1 - eps) / 2
    delta = eta / 8
    internal = tau + (1 - beta) * max(sigma - tau, Q(0))
    leaf = sigma + beta * (1 - sigma)
    if kappa is None:
        kappa = Q(math.floor(minimum * KGRID), KGRID)
    margins = dict(original_prefix=1 - eps * (1 + c), coordinate_movement=a, compact_phase_layer=minimum,
                   bulk_exposure=a, Gaussian_arithmetic=min(1 - eps - delta, r - delta), scalar_work=1 - eps - delta,
                   dimension=eps)
    slacks = dict(
        bit_positive=a, complex_above_bit=b - a, complex_below_one_over32=Q(1, 32) - b, beta_positive=beta,
        beta_below_one=1 - beta, leaf_saving_above_bit=(1 - beta) * b - a, q_positive=q,
        q_below_internal=1 - internal - q, q_below_leaf=1 - leaf - q, c_positive=c, c_below_one=1 - c,
        q_below_reservations=c - q, lambda_above_tau=lam - tau, lambda_above_sigma=lam - sigma,
        lambda_above_internal=lam - internal, lambda_prime_above_lambda=lp - lam, compact_leaf=lp - leaf,
        compact_reservations=lp - (1 - c), lambda_prime_below_one=q, epsilon_positive=eps, epsilon_below_one=1 - eps,
        guard_width=1 - eps, K_geometry=1 - eps * (1 + c), K_dominates_log=eps * c, record_suffix=1 - eps,
        phase_local=1 - eps - delta, phase_boundary=r - delta, gamma_sublinear=1 - eps - r,
        cell_above_band=eps - (1 - r) / 2, prime_interval_packing=1 - eps, alpha_positive=r, alpha_below_one=1 - r,
        alpha_below_one_fourth=Q(1, 4) - r, delta_positive=delta, delta_below_one_eighth=Q(1, 8) - delta,
        short_record_fallback=eps - a, small_field_exposure=1 - eps - minimum,
        artificial_boundary=8 - eps + r - delta - minimum, literal_scalar_guard=Q(bridge['literal_gap']),
        row_product_gap=bridge['degree_gap'])
    slacks.update({name + '_above_kappa': val - kappa for name, val in margins.items()})
    assert len(slacks) == 47 and len(margins) == 7
    bad = sorted(k for k, x in slacks.items() if not x > 0)
    # identities of the assembly: the compact phase layer is the binding margin, eta below the original prefix
    assert min(margins.values()) == minimum and margins['original_prefix'] - minimum == eta
    return dict(a=a, b=b, tau=tau, sigma=sigma, q=q, c=c, eps=eps, lam=lam, lp=lp, r=r, delta=delta,
                internal=internal, leaf=leaf, minimum=minimum, kappa=kappa, absorption_gap=minimum - kappa,
                margins=margins, strict=slacks, bad=bad, ok=not bad)


def certify(bit_inv, cx_inv, kappa=None):
    b, c = bit_side(bit_inv), complex_side(cx_inv)
    bridge = finite_bridge(c, cx_inv)
    a = min(b['a_bit'], (1 - BETA) * c['a_c'] - CX_GAP)
    assert a <= b['a_bit'], 'supported bit interface'
    k = assembly(a, c['a_c'], bridge, kappa)
    for side, d, key in (('bit', b, 'astar'), ('complex', c, 'a_c')):
        inv = bit_inv if side == 'bit' else cx_inv
        if 'claim' in inv and key in inv['claim']:
            assert Q(inv['claim'][key]) == d[key], (side, key, d[key])
    if 'claim' in bit_inv and 'a_bit' in bit_inv['claim']:
        assert Q(bit_inv['claim']['a_bit']) == b['a_bit']
    return dict(bit=b, cx=c, bridge=bridge, k=k, bit_binds=b['a_bit'] < (1 - BETA) * c['a_c'] - CX_GAP)


# ---------------------------------------------------------------- PR #144 reproduction
PR144_BIT = dict(saving='4617656/10000000000', theta='1/1000', a_old='384599/10000000000')
PR144_CX = dict(saving='4856569/10000000000')


def pr144_inputs():
    bit = dict(load(os.path.join(PR144, 'paired-cube-bit-input.json')), side='bit', cert=PR144_BIT)
    cx = dict(load(os.path.join(PR144, 'paired-cube-complex-input.json')), side='complex', cert=PR144_CX)
    return bit, cx


def reproduce_pr144(verbose=True):
    """PR #144's published values from its published inventories (read as data): both moment gaps, the stopped bit
    saving, kappa, the minimum margin, the absorption gap, the row stock, the group order, and every one of its 47
    strict constraints and 7 margins, all exactly."""
    bit, cx = pr144_inputs()
    saved = load(os.path.join(PR144, 'paired-cube-input.json'))
    net = load(os.path.join(PR144, 'paired-cube-network.json'))
    r = certify(bit, cx, kappa=Q(saved['kappa']))
    k, A = r['k'], net['assembly']
    checks = {
        'bit moment gap': r['bit']['gap'] == Q(saved['bit_moment_gap']),
        'complex moment gap': r['cx']['gap'] == Q(saved['complex_moment_gap']),
        'stopped bit saving': r['bit']['a_bit'] == Q(net['bit']['effective_saving']),
        'complex group order': r['cx']['V'] == saved['complex_group_order'],
        'row stock': (r['bridge']['degree'], r['bridge']['coefficient']) == (saved['row_degree'], saved['row_coefficient']),
        'bit children': clean(saved['bit_child_multiplicities']) == r['bit']['p']['hist'],
        'complex children': clean(saved['complex_child_multiplicities']) == r['cx']['p']['hist'],
        'minimum margin': k['minimum'] == Q(saved['assembly_minimum']),
        'absorption gap': k['absorption_gap'] == Q(saved['absorption_gap']),
        'kappa = floor of the minimum on the 1e-10 grid': Q(math.floor(k['minimum'] * KGRID), KGRID) == Q(saved['kappa']),
        '47 strict constraints, exactly': set(k['strict']) == set(A['strict_constraints'])
                                           and all(k['strict'][n] == Q(x) for n, x in A['strict_constraints'].items()),
        '7 margins, exactly': set(k['margins']) == set(A['margins'])
                              and all(k['margins'][n] == Q(x) for n, x in A['margins'].items()),
        'semantic gap': r['bridge']['literal_gap'] == net['finite_bridge']['semantic']['strict_literal_gap'],
        'all constraints positive': k['ok'],
    }
    if verbose:
        print('PR #144 reproduction (published inventories, read as data):')
        for name, ok in checks.items():
            print('  %-50s %s' % (name, 'PASS' if ok else 'FAIL'))
        print('  kappa = %s (%.7e); a_bit = %s, a_c = %s' % (k['kappa'], float(k['kappa']), r['bit']['a_bit'], r['cx']['a_c']))
    assert all(checks.values()), [n for n, ok in checks.items() if not ok]
    return r


# ---------------------------------------------------------------- report and freeze
def report(r, bit_inv, cx_inv):
    b, c, k = r['bit'], r['cx'], r['k']
    for name, d, inv in (('bit', b, bit_inv), ('complex', c, cx_inv)):
        p = d['p']
        print('%s: %s' % (name, inv.get('source', '')))
        print('  h=%d v=%d R=%d loss=%d; m=%d W=%d s=%d D=%d; widths=%d max=%d edges=%d' % (
            p['h'], p['v'], p['R'], p['loss'], p['m'], p['W'], p['s'], p['D'], len(p['hist']), p['maxchild'], p['edges']))
    print('  bit a* = %s (%.12e), gap %.3e, next point %s' % (b['astar'], float(b['astar']), float(b['gap']),
          'pinned' if b['gap_next'] is None else 'fails (gap %.3e)' % float(b['gap_next'])))
    print('  stopped a_bit = (1-%s) a* + %s * %s = %s (%.12e)' % (b['theta'], b['theta'], b['a_old'], b['a_bit'], float(b['a_bit'])))
    print('  complex a_c = %s (%.10e), gap %.3e, next point %s' % (c['a_c'], float(c['a_c']), float(c['gap']),
          'pinned' if c['gap_next'] is None else 'fails (gap %.3e)' % float(c['gap_next'])))
    print('assembly: a = %s (%s binds); minimum margin eps q = %.11e' % (k['a'], 'bit' if r['bit_binds'] else 'complex', float(k['minimum'])))
    print('kappa = %s (%.7e = 2^%.4f), ok=%s, 47 constraints and 7 margins; above 2^-12: %s' % (
        k['kappa'], float(k['kappa']), math.log2(float(k['kappa'])), k['ok'], k['kappa'] > TARGET))
    if k['bad']:
        print('  failing rows:', k['bad'])


def fr(q):
    q = Q(q)
    return [q.numerator, q.denominator]


def _dump_gz(path, d):
    # mtime=0 keeps the gzip bytes reproducible
    with open(path, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as f:
        f.write(json.dumps(d, sort_keys=True).encode())


def lean_input(r, bit_inv, cx_inv):
    """lean/round9-histograms.json: both moments (with the bit fallback), the stopped mix, the finite bridge inputs
    and the assembly"""
    b, c, k, br = r['bit'], r['cx'], r['k'], r['bridge']
    side = lambda d, a: dict(a=fr(a), m=d['p']['m'], W=d['p']['W'], s=d['p']['s'], edges=d['p']['edges'],
                             hist=sorted(d['p']['hist'].items()))
    return dict(kind='cover144',
                bit=dict(side(b, b['astar']), label=bit_inv.get('source', 'bit'), fallback=True,
                         stop=dict(theta=fr(b['theta']), a_old=fr(b['a_old']), ab=fr(b['a_bit']))),
                cx=dict(side(c, c['a_c']), label=cx_inv.get('source', 'complex'), fallback=False,
                        h=c['p']['h'], v=c['p']['v'], R=c['p']['R'], cAdd=cx_inv['c'], mops=cx_inv['total_M_operations'],
                        maxchild=c['p']['maxchild'], halving=br['halving_degree'], wire_bits=br['wire_bits']),
                bad=fr(BAD), beta=fr(BETA), eta=fr(ETA), cx_gap=fr(CX_GAP), row_reserve=ROW_RESERVE,
                row_degree=ROW_DEGREE, kappa=fr(k['kappa']))


def freeze(r, bit_inv, cx_inv):
    os.makedirs(FROZEN, exist_ok=True)
    bit_inv = dict(bit_inv, claim=dict(astar=str(r['bit']['astar']), a_bit=str(r['bit']['a_bit'])))
    cx_inv = dict(cx_inv, claim=dict(a_c=str(r['cx']['a_c'])))
    _dump_gz(os.path.join(FROZEN, 'bit_inventory.json.gz'), bit_inv)
    _dump_gz(os.path.join(FROZEN, 'cx_inventory.json.gz'), cx_inv)
    with open(os.path.join(ROOT, 'lean', 'round9-histograms.json'), 'w') as f:
        json.dump(lean_input(r, bit_inv, cx_inv), f)
    print('froze certificates/round9/{bit,cx}_inventory.json.gz and lean/round9-histograms.json')


def main(argv):
    do_freeze = '--freeze' in argv
    only144 = '--pr144' in argv
    args = [a for a in argv if not a.startswith('--')]
    reproduce_pr144()
    if only144:
        return
    if args:
        bpath, cpath = args
    else:
        bpath, cpath = os.path.join(FROZEN, 'bit_inventory.json.gz'), os.path.join(FROZEN, 'cx_inventory.json.gz')
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
