"""Independent rational certificate checker. Python 3.9+, standard library only.

This verifies finite moment inequalities and assembly arithmetic, conditional
on the supplied network histograms and upstream algorithmic/analytic arguments.
It does not prove a complete multiplication algorithm or establish priority.
"""
import argparse
import json
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def floor_grid(x, denominator):
    return Q(x.numerator * denominator // x.denominator, denominator)


@lru_cache(None)
def log_bounds(x, terms=72):
    """Atanh series with a positive geometric bound for its omitted tail."""
    x = Q(x)
    require(x >= 1, 'log input below one')
    k = 0
    while x > 2:
        x /= 2
        k += 1

    def series(y):
        z = (y - 1) / (y + 1)
        require(0 <= z <= Q(1, 3), 'log reduction failed')
        low = 2 * sum((z ** (2*j+1) / (2*j+1) for j in range(terms)), Q(0))
        tail = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
        return low, low + tail

    lo, hi = series(x)
    l2, u2 = series(Q(2))
    den = 10**60
    return floor_grid(lo+k*l2, den), floor_grid(hi+k*u2, den)+Q(1, den)


def exp_bounds(x, degree=12):
    """Positive Taylor polynomial; subsequent term ratios are <= x/(d+2)."""
    x = Q(x)
    require(0 <= x < degree+2, 'exp tail not contracting')
    term = Q(1)
    total = term
    for j in range(1, degree+1):
        term *= x / j
        total += term
    omitted = term*x/(degree+1)
    return total, total + omitted/(1-x/(degree+2))


def normalize_network(network):
    m, W = network['m'], network['W']
    hist = {int(w): int(c) for w, c in network['hist'].items()}
    require(m > 1 and W > 0, 'nonpositive network size')
    require(all(0 < w < m and c > 0 for w, c in hist.items()), 'invalid child')
    total = sum(w*c for w, c in hist.items())
    require(total == network['rank_sum'], 'rank sum mismatch')
    require(total < m*W, 'no rank deficit')
    return m, W, hist


def moment_bounds(network, saving):
    """Enclose sum c_w * (w/m)**(1-saving) / W without floating point."""
    saving = Q(saving)
    require(0 <= saving < 1, 'invalid saving')
    m, W, hist = normalize_network(network)
    low, high = Q(0), Q(0)
    for w, count in hist.items():
        l, u = log_bounds(Q(m, w))
        el, _ = exp_bounds(saving*l)
        _, eu = exp_bounds(saving*u)
        weight = Q(count*w, W*m)
        low += weight*el
        high += weight*eu
    return low, high


def approximate_root(network):
    """Discovery only. Every accepted rational is separately certified."""
    m, W, hist = normalize_network(network)
    with localcontext() as ctx:
        ctx.prec = 80
        data = [(Decimal(w*c)/Decimal(m*W), (Decimal(m)/Decimal(w)).ln())
                for w, c in hist.items()]
        lo, hi = Decimal(0), Decimal(1)
        for _ in range(200):
            mid = (lo+hi)/2
            value = sum(weight*(mid*log).exp() for weight, log in data)
            if value < 1:
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2


def certify_root(network, denominator=10**24):
    guess = approximate_root(network)
    with localcontext() as ctx:
        ctx.prec = 80
        tick = int(guess*denominator)
    for _ in range(8):
        candidate = Q(tick, denominator)
        low, high = moment_bounds(network, candidate)
        if high < 1:
            next_low, next_high = moment_bounds(network, candidate+Q(1, denominator))
            if next_low > 1:
                return dict(lower=str(candidate), upper=str(candidate+Q(1, denominator)),
                            upper_moment_gap=str(1-high), next_lower_excess=str(next_low-1),
                            approximate_root=str(guess))
            tick += 1
        else:
            tick -= 1
    raise ValueError('could not bracket root at requested resolution')


def check_jain_assembly(p):
    ab, ac, beta, eps, x, c1, kappa, gap, delta = [Q(p[k]) for k in
        ('a_bit','a_complex','beta','epsilon','x','C1','kappa','lambda_gap','delta')]
    require(0 < ab < ac < 1 and 0 < beta < 1 and gap > 0 and x > 0, 'assembly domains')
    tau, sigma = 1-ab, 1-ac
    chi = tau+(1-beta)*max(sigma-tau, Q(0))
    leaf = sigma+beta*(1-sigma)
    lam = max(tau,sigma,chi)+gap
    lp = max(lam,leaf)+gap
    _, log_upper = log_bounds(Q(p['complex_rank_sum']))
    slacks = {
        'lambda_above_internal': lam-max(tau,sigma,chi),
        'lambda_prime_above_lambda_and_leaf': lp-max(lam,leaf),
        'lambda_prime_below_one': 1-lp,
        'reservations': 1-eps*(2-lp),
        'guard': 1+x-eps*c1,
        'guard_constant': c1-2-p['complex_m']*log_upper,
        'precision': 1+x-3*eps,
        'epsilon_positive': eps,
        'epsilon_below_one': 1-eps,
        'delta_positive': delta,
        'delta_below_eighth': Q(1,8)-delta,
    }
    margins = {
        'prefix_moves': 1-eps,
        'simultaneous_butterflies': eps*(1-lp),
        'fine_bit_exposures': 1-eps-(1+x)*delta,
        'Gaussian_maps': 1-eps-(1+x)*delta,
        'chirps_twists_scalar_products': 1-eps-(1+x)*delta,
        'packed_polynomial_products': eps,
        'reserved_axes': 1-eps-(1+x)*delta,
        'CRT_reversal': max(eps,1-eps)*ab,
    }
    slacks.update({k+'_above_kappa': v-kappa for k,v in margins.items()})
    require(kappa > 0 and all(v > 0 for v in slacks.values()), 'Jain assembly inequality failed')
    return slacks, margins


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', nargs='?', default=str(Path(__file__).with_name('certificate.json')))
    args = parser.parse_args()
    data = json.loads(Path(args.certificate).read_text(encoding='utf-8'))
    for name, network in data['networks'].items():
        cert = data['root_certificates'][name]
        lower, upper = Q(cert['lower']), Q(cert['upper'])
        _, hi = moment_bounds(network,lower)
        lo, _ = moment_bounds(network,upper)
        require(hi < 1 < lo, name+' root enclosure failed')
        print(name, 'root bracket', float(lower), float(upper), 'PASS')
    for name, result in data.get('jain_assemblies',{}).items():
        slacks, margins = check_jain_assembly(result)
        net = data['networks'][result['bit_network']]
        require(moment_bounds(net,Q(result['a_bit']))[1] < 1, 'bit saving invalid')
        net = data['networks'][result['complex_network']]
        require(moment_bounds(net,Q(result['a_complex']))[1] < 1, 'complex saving invalid')
        print(name, 'kappa', result['kappa'], 'all',len(slacks),'strict constraints PASS')
    if data.get('community_assemblies'):
        from community_assembly import assembly
        for name, result in data['community_assemblies'].items():
            ab, ac, kappa, backoff = [Q(result[k]) for k in ('a_bit','a_complex','kappa','h')]
            bit_name=result.get('bit_network','community_bit')
            require(moment_bounds(data['networks'][bit_name],ab)[1] < 1, 'community bit moment')
            require(moment_bounds(data['networks']['community_complex'],ac)[1] < 1, 'community complex moment')
            bridge=result.get('finite_bridge',data['community_finite_bridge'])
            require(bridge['bit']['m']==data['networks'][bit_name]['m'] and
                    bridge['bit']['W']==data['networks'][bit_name]['W'] and
                    bridge['bit']['maxchild']==max(map(int,data['networks'][bit_name]['hist'])),
                    'bridge does not match bit network')
            value = assembly(bridge,ab,kappa,h=backoff,a_complex=ac)
            print(name, 'kappa', str(kappa), 'all',len(value['constraints']),'strict constraints PASS')
    print('Finite arithmetic verified; upstream network and multiplication proofs remain assumptions.')


if __name__ == '__main__':
    main()
