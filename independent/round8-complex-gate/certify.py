"""Independent certificate from a replay histogram {rank: count}: W, s, exact moment, a_c, kappa.

Moment (round-6 complex condition, PR #10 batched recursion, sigma = 1 - a): M(a) = sum_r (r n_r/(W m)) (m/r)^a < 1.
Rigorous two-sided bounds with a method different from cert.py: ln z for z in [1, 2) by the series
sum_k u^k/k, u = 1 - 1/z <= 1/2 (tail <= u^(K+1)/((K+1)(1-u))), ln 2 = sum_k 1/(k 2^k); exp by its Taylor
polynomial plus the tail bound x^(K+1)/(K+1)! / (1 - x/(K+2)). a_c = largest point of the 1/(2 10^11) grid
with upper(M) < 1; also checks lower(M) > 1 at the next grid point (so the grid point is exactly the largest).
Usage: python3 certify.py hist.pkl H R"""
import sys, pickle, math, os
from fractions import Fraction as Q

P = 10 ** 45
def dn(x): return Q((x.numerator * P) // x.denominator, P)
def upr(x): return Q(-((-x.numerator * P) // x.denominator), P)
K = 70
def ln_z(z):                      # (lo, hi) for ln z, z in [1, 2)
    u = 1 - 1 / z; s = Q(0); p = Q(1)
    for k in range(1, K + 1): p *= u; s += p / k
    return dn(s), upr(s + p * u / ((K + 1) * (1 - u)))
def ln2():
    s = Q(0)
    for k in range(1, K + 1): s += Q(1, k * 2 ** k)
    return dn(s), upr(s + Q(1, (K + 1) * 2 ** K))
L2 = ln2()
def ln_y(y):
    k = 0
    while y >= 2: y /= 2; k += 1
    a, b = ln_z(y); return k * L2[0] + a, k * L2[1] + b
def exp_b(x, hi):
    s = Q(0); t = Q(1)
    for k in range(0, 12): s += t; t = t * x / (k + 1)
    return upr(s + t / (1 - x / 13)) if hi else dn(s)       # t = x^12/12!, x << 1


def certify(H, W, m, grid=2 * 10 ** 11):
    ws = [(Q(r * n, W * m), ln_y(Q(m, r))) for r, n in H.items()]
    M = lambda a, hi: sum(w * exp_b(a * (U[1] if hi else U[0]), hi) for w, U in ws)
    lo, hi = 0, grid // 1000
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if M(Q(mid, grid), True) < 1: lo = mid
        else: hi = mid
    a = Q(lo, grid)
    return a, M(a, True), M(Q(lo + 1, grid), False)


if __name__ == '__main__':
    H = pickle.load(open(sys.argv[1], 'rb')); h = int(sys.argv[2]); R = int(sys.argv[3])
    m = h * h; v = math.comb(h, 3); N = v * v; W = 2 * N + 2 * v * R
    s = sum(r * n for r, n in H.items()); L = s - (W * m - N)
    print(dict(W=W, W_from_hist_wires='2N + 2vR', s=s, L_copies=L, maxrank=max(H), guard_ok=max(H) <= m - 1))
    a, Mup, Mlo_next = certify(H, W, m)
    print(dict(a_c=str(a), a_c_float=float(a), upper_M_at_a=float(Mup), lower_M_next=float(Mlo_next),
               certified=Mup < 1, largest=Mlo_next > 1))
    sys.path.insert(0, os.path.expanduser('~/integer-mult-kappa/scripts'))
    import certificate_round3 as C
    r = C.evaluate(Q(6338633, 5 * 10 ** 10), a, Q(1, 1000), 'crude', m_c=576, s_c=s)
    print(dict(ok=r['ok'], kappa=str(r['kappa']), kappa_float=float(r['kappa']), binding=r['binding'], bad=r['bad'], eps=float(r['eps'])))
