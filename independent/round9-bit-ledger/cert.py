"""Independent rational moment certificate for a bit child histogram (own code, written separately from omit.py).
Moment: F(a) = (1/W) [ sum_w n_w (w/m)^(1-a) + fb (1/m)^(1-a) ], fb = 1e-16 * 32 m^2 * (#children) (full rare-class
fallback, nothing displaced). a is certified iff an UPPER bound of F(a) is < 1; the next grid point is rejected iff a
LOWER bound of F(a + 1/den) is >= 1 (a genuine rejection, not merely a failed upper bound).
(w/m)^(1-a) = (w/m) exp(a ln(m/w)); ln bounds from ln x = k ln2 + 2 atanh(z) (partial sums below, geometric tail above),
exp bounds from Taylor partial sums with a geometric tail.
Then the stopped saving a_bit = (1 - 1/1000) a* + (1/1000) 384599/10^10 and #144's assembly (paired-cube-network.json):
a = min(a_bit, (1-beta) a_c - 1e-10), q = a(1 - 2 eta), c = q(1 + eta), kappa = (1 - eta) q / (1 + c + q)."""
import json, sys, os, math
if hasattr(sys, 'set_int_max_str_digits'): sys.set_int_max_str_digits(0)   # the published router integers are long
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__))
BAD, THETA, AOLD = Q(1, 10**16), Q(1, 1000), Q(384599, 10**10)
BETA, ETA, AC = Q(1, 10**6), Q(1, 10**8), Q(4856569, 10**10)

def atanh_bounds(z, n=40):
    s = sum(z ** (2 * j + 1) / (2 * j + 1) for j in range(n))
    tail = z ** (2 * n + 1) / ((2 * n + 1) * (1 - z * z))
    return s, s + tail
L2 = tuple(2 * b for b in atanh_bounds(Q(1, 3)))
def ln_bounds(x):
    x = Q(x); assert x >= 1; k = 0
    while x > 2: x /= 2; k += 1
    lo, hi = atanh_bounds((x - 1) / (x + 1)) if x > 1 else (Q(0), Q(0))
    return k * L2[0] + 2 * lo, k * L2[1] + 2 * hi
def exp_bounds(u, n=6):
    assert 0 <= u < 1
    s = sum(u ** k / math.factorial(k) for k in range(n))
    return s, s + u ** n / math.factorial(n) / (1 - u / (n + 1))

def moment(hist, W, m, a, lb, side):
    E = sum(hist.values()); fb = BAD * 32 * m * m * E
    items = list(hist.items()) + [(1, fb)]; tot = Q(0)
    for w, n in items:
        l = lb[w][0 if side == 'lo' else 1]
        e = exp_bounds(a * l)[0 if side == 'lo' else 1]
        tot += Q(n) * Q(w, m) * e
    return tot / W

def certify(hist, W, m, den):
    hist = {int(w): n for w, n in hist.items() if n}
    assert 1 <= min(hist) and max(hist) < m and W * m > sum(w * n for w, n in hist.items())
    lb = {w: ln_bounds(Q(m, w)) for w in set(hist) | {1}}
    fl = lambda a: (sum(n * (w / m) ** (1 - a) for w, n in hist.items()) + float(BAD * 32 * m * m * sum(hist.values())) * (1 / m) ** (1 - a)) / W
    lo, hi = 0.0, 0.01
    for _ in range(100):
        mid = (lo + hi) / 2; lo, hi = (mid, hi) if fl(mid) < 1 else (lo, mid)
    k = int(lo * den) - 3
    while moment(hist, W, m, Q(k + 1, den), lb, 'hi') < 1: k += 1
    while not moment(hist, W, m, Q(k, den), lb, 'hi') < 1: k -= 1
    up = moment(hist, W, m, Q(k, den), lb, 'hi'); nxt = moment(hist, W, m, Q(k + 1, den), lb, 'lo')
    return dict(a=Q(k, den), F_upper_at_a=up, F_lower_at_next=nxt, certified=up < 1, next_rejected=nxt >= 1)

def stopped(a): return (1 - THETA) * a + THETA * AOLD
def assembly(ab, ac=AC):
    a = min(Q(ab), (1 - BETA) * ac - Q(1, 10**10)); q = a * (1 - 2 * ETA); c = q * (1 + ETA)
    k = (1 - ETA) * q / (1 + c + q)
    tau, sigma = 1 - Q(ab), 1 - ac; leaf = sigma + BETA * (1 - sigma); lam, lamp = 1 - a * (1 - ETA), 1 - q
    cons = dict(lambda_above_tau=lam - tau, lambda_above_sigma=lam - sigma, lambda_above_leaf=lam - leaf, lambda_prime_above_lambda=lamp - lam,
                lambda_prime_below_one=1 - lamp, complex_above_bit=ac - Q(ab), bit_positive=Q(ab))
    return dict(a=a, q=q, c=c, kappa=k, kappa_floor_1e10=math.floor(k * 10**10), constraints=cons,
                all_constraints_positive=all(x > 0 for x in cons.values()))

if __name__ == '__main__':
    out = {}
    # 1. reproduce #144 from its published child histogram (read as data)
    P = json.load(open(os.path.join(HERE, '..', '..', 'certificates', 'round9', 'pr144', 'paired-cube-network.json')))
    hb = {int(w): n for w, n in P['bit']['counts']['child_multiplicities'].items()}
    c144 = certify(hb, P['bit']['counts']['W_per_vertex'], 69, 1250000000)
    ab144 = stopped(c144['a']); as144 = assembly(Q(P['assembly']['parameters']['a_bit']))
    out['pr144'] = dict(coarse=str(c144['a']), coarse_matches=c144['a'] == Q(P['bit']['coarse_saving']),
                        next_rejected=c144['next_rejected'], a_bit=str(ab144), a_bit_matches=ab144 == Q(P['assembly']['parameters']['a_bit']),
                        q_matches=as144['q'] == Q(P['assembly']['parameters']['q']), c_matches=as144['c'] == Q(P['assembly']['parameters']['c']),
                        kappa=float(as144['kappa']), kappa_floor_matches=Q(as144['kappa_floor_1e10'], 10**10) == Q(P['kappa']),
                        constraints_match={k: (as144['constraints'][k] == Q(P['assembly']['strict_constraints'][k]) if k in P['assembly']['strict_constraints'] else None) for k in as144['constraints']})
    c12 = certify(hb, 32408, 69, 10**12); out['pr144']['grid_1e12'] = dict(a=str(c12['a']), a_f=float(c12['a']), stopped=float(stopped(c12['a'])),
                                                                          kappa=float(assembly(stopped(c12['a']))['kappa']))
    print('#144 reproduction:', json.dumps(out['pr144'], indent=1), flush=True)
    # 2. our ledgers
    for name in sys.argv[1:]:
        L = json.load(open(os.path.join(HERE, name)))
        c = certify(L['hist'], L['W'], L['m'], 10**12); ab = stopped(c['a']); asb = assembly(ab)
        out[name] = dict(W=L['W'], m=L['m'], rank=L['rank'], children=sum(L['hist'].values()), astar=str(c['a']), astar_f=float(c['a']),
                         certified=c['certified'], next_rejected=c['next_rejected'], F_upper=float(c['F_upper_at_a']), F_lower_next=float(c['F_lower_at_next']),
                         a_bit=str(ab), a_bit_f=float(ab), kappa=str(asb['kappa']), kappa_f=float(asb['kappa']), kappa_floor_1e10=asb['kappa_floor_1e10'],
                         constraints={k: float(x) for k, x in asb['constraints'].items()}, all_constraints_positive=asb['all_constraints_positive'],
                         bit_binds=Q(ab) < (1 - BETA) * AC - Q(1, 10**10))
        print(name, json.dumps(out[name], indent=1), flush=True)
    json.dump(out, open(os.environ.get('CERT_OUT', '/tmp/round9_cert.json'), 'w'), indent=1, default=str)
