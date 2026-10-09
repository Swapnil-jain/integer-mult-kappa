"""Two-sided rational moment check (standing gate; imports nothing but the standard library).
The existing certificates reject the next grid point with the same UPPER bound they use to accept. That shows the
upper bound is tight, but not that the true moment is >= 1 there. This checks both directions rigorously:
  accept a:          F_upper(a) < 1     (ln by atanh series + geometric tail, exp(u) <= 1 + u + u^2 / (2 (1 - u/3)))
  reject a + 1/den:  F_lower(a + 1/den) >= 1   (ln by the truncated atanh series, exp(u) >= 1 + u + ... + u^5/120)
F(a) = (1/W) sum_r n_r (r/m) exp(a ln(m/r)).  Input: JSON with a histogram ('hist' or 'C' or 'child_histogram') and
W ('W' or 'W_per_vertex'), optional 'scale'.  Usage: python3 twosided.py IN.json M A [den=10**12]"""
import sys, json, gzip
from fractions import Fraction as Q


def at(z, N, tail):
    s = 2 * sum(z ** (2 * j + 1) / (2 * j + 1) for j in range(N))
    return s + (2 * z ** (2 * N + 1) / ((2 * N + 1) * (1 - z * z)) if tail else 0)


def ln_bound(x, upper, N=30):
    x, k = Q(x), 0
    while x >= 2: x /= 2; k += 1
    return k * at(Q(1, 3), N, upper) + at((x - 1) / (x + 1), N, upper)


def main():
    p = sys.argv[1]
    d = json.load((gzip.open if p.endswith('.gz') else open)(p, 'rt'))
    m, a = int(sys.argv[2]), Q(sys.argv[3])
    den = int(eval(sys.argv[4].split('=', 1)[1])) if len(sys.argv) > 4 else 10 ** 12
    hist = d.get('hist') or d.get('C') or d.get('child_histogram')
    if isinstance(hist, str): hist = eval(hist)
    W = Q(d.get('W', d.get('W_per_vertex'))); k = int(d.get('scale', 1))
    H = {int(r): Q(n) / k for r, n in hist.items() if Q(n)}; W = W / k
    assert all(0 < r < m for r in H), 'child widths must lie in (0, m)'
    lu = {r: ln_bound(Q(m, r), True) for r in H}; ll = {r: ln_bound(Q(m, r), False) for r in H}
    def Fu(b):
        tot = Q(0)
        for r, n in H.items():
            u = b * lu[r]; assert 0 <= u < 3; tot += n * Q(r, m) * (1 + u + u * u / (2 * (1 - u / 3)))
        return tot / W
    def Fl(b):
        tot = Q(0)
        for r, n in H.items():
            u = b * ll[r]; tot += n * Q(r, m) * sum(u ** i / f for i, f in enumerate((1, 1, 2, 6, 24, 120)))
        return tot / W
    nxt = a + Q(1, den)
    ok1, ok2 = Fu(a) < 1, Fl(nxt) >= 1
    print('%s accept a=%s (F_upper < 1)' % ('PASS' if ok1 else 'FAIL', a))
    print('%s reject next %s by a LOWER bound (F_lower >= 1)' % ('PASS' if ok2 else 'FAIL', nxt))
    print('RESULT twosided: %s' % ('ALL PASS' if ok1 and ok2 else 'FAIL'))
    sys.exit(0 if ok1 and ok2 else 1)


main()
