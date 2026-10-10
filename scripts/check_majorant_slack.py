"""How much the exponential majorant costs the round-eleven bit certificate.

The bit moment is certified with E3(x) = 1 + x + x^2/(2(1 - x/3)) >= exp(x) (certificate_round9.exp_upper) on a
1e-12 grid for a*. This script recomputes the same moment with exp and ln evaluated to 80 decimal digits and
bisects for the true root, with and without the fallback term. If the true root is less than one grid step above
the certified a*, no tighter majorant for exp can move the certified saving, and the round-seven refinement of
PR #2 (a degree-eight Taylor majorant) has nothing to give here.

Floating point, named as such: this is an [exp]-status check, not a rational certificate.

    python3 scripts/check_majorant_slack.py            # round eleven bit inventory
    python3 scripts/check_majorant_slack.py certificates/round10/bit_inventory.json.gz
"""
import gzip
import json
import os
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as Q

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certificate_round9 as c9  # noqa: E402
import certificate_round10 as c10  # noqa: E402

getcontext().prec = 80


def D(x):
    x = Q(x)
    return Decimal(x.numerator) / Decimal(x.denominator)


def true_gap(p, a, fallback):
    """1 - M(a) - fallback(a) with exact-to-80-digits exp and ln (no majorant)."""
    m, W = p['m'], p['W']
    a = D(a)
    M = sum(D(n * r) / (D(W) * m) * (D(m) / D(r)).ln().__mul__(a).exp() for r, n in p['hist'].items())
    fb = D(c9.BAD) * D(32 * m * m * p['edges']) / (D(W) * m) * (D(m).ln() * a).exp() if fallback else Decimal(0)
    return 1 - M - fb


def root(p, lo, fallback, span=Q(1, 10 ** 5), steps=90):
    hi = lo + span
    assert true_gap(p, lo, fallback) > 0 > true_gap(p, hi, fallback)
    for _ in range(steps):
        mid = (lo + hi) / 2
        if true_gap(p, mid, fallback) > 0:
            lo = mid
        else:
            hi = mid
    return lo


def main(path):
    inv = json.load(gzip.open(path))
    b = c10.bit_side(inv)
    p, astar = b['p'], b['astar']
    grid = Q(1, inv['cert']['den'])
    print('inventory      ', path)
    print('certified a*   ', astar, '=', float(astar), ' on grid', grid)
    print('certificate gap', float(b['gap']), ' next point', float(b['gap_next']))
    print('E3 moment      ', float(b['moment']), ' fallback term', float(b['added']))
    g = true_gap(p, astar, True)
    print('true gap at a* ', g)
    r = root(p, astar, True)
    r0 = root(p, astar, False)
    slack = r - astar
    print('true root      ', float(r), ' slack over a* =', float(slack), ' grid step', float(grid))
    print('without fallback root', float(r0), ' fallback costs', float(r0 - r), 'of a*')
    moves = slack >= grid
    print('a tighter exp majorant could move a* by at least one grid step:', moves)
    return 0 if not moves else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else 'certificates/round11/bit_inventory_p12.json.gz'))
