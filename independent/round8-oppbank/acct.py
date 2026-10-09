"""Items 2, 3, 5, 6 (accounting part): independent histogram vs the audited one, independent rigorous moment bound,
stopped composition, kappa, beta sweep, negative controls."""
import sys, os, pickle, math
from fractions import Fraction as Q
from decimal import Decimal, getcontext
HERE = os.path.dirname(os.path.abspath(__file__))
W2 = os.path.normpath(os.path.join(HERE, '..', '..'))
OUT = os.environ.get('GATE_OUT', '/tmp/round8-bitgate'); os.makedirs(OUT, exist_ok=True)
sys.path[:0] = [os.path.join(W2, 'scripts'), os.path.join(W2, 'independent', 'two-stage-bit')]
import moment
from certificate_round3 import evaluate
G = pickle.load(open(os.path.join(OUT, 'hist_gate.pkl'), 'rb'))
import gzip, json
with gzip.open(os.path.join(W2, 'certificates', 'round8', 'bit_hist.json.gz'), 'rt') as _f: AUD = json.load(_f)
AUD['hist'] = {int(w): int(n) for w, n in AUD['hist']}  # the frozen histogram the certificate uses
ok = {}
def rep(k, c): ok[k] = bool(c); print('%-100s %s' % (k, 'PASS' if c else 'FAIL'), flush=True)
h, m, v, N, W, s = G['h'], G['m'], G['v'], G['N'], G['W'], G['s']
H = G['hist']
rep('A1. independent one-run histogram == audited onerun histogram (46 widths)', H == AUD['hist'])
rep('A2. W, s unchanged vs audited (W=%d, s=%d)' % (W, s), (W, s) == (AUD['W'], AUD['s']))
D = W * m - s
rep('A3. D = W m - s = N - L = %d' % D, D == N - G['L'] == 1344189)
rep('A4. rank sum exact, max child %d < m, min child %d >= 1' % (max(H), min(H)), sum(w * c for w, c in H.items()) == s and max(H) < m and min(H) >= 1)
# rebuild the OLD (round-7.1 plain) histogram from the same objects to show only the compile changed
def inner(r): return [1] * (h - r) + [2 * r - h] if 2 * r > h else [1] * r
OLD = {}
def add(d, w, c): d[w] = d.get(w, 0) + c
rk, cr = G['rk'], G['corner_rank']
# per-class rank multisets of one stage (shared by old and new compile)
cls = []
for r, c in cr.items(): cls.append(('aux', m - r, v * c, r))
for r, c in rk.items(): cls.append(('step', r, v * c, None))
NEWchk = {}
# old: aux [m-2r]+inner(r); step inner(r); centre inner(h-1); y/x inner(x); entrance [m-4h+2]+[h-2,1^(h+1)]; copy [1]
# new: aux [m-r]; step [r]; ...  -- rebuilt from H by subtracting the shared y/x/centre parts is circular, so rebuild old from
# the per-step multiset that the new histogram encodes: every new width w (other than aux m-r, entrance, copy) is one step of rank w.
aux_new = {m - r: 2 * v * c for r, c in cr.items()}
rest = dict(H)
for w, c in aux_new.items(): rest[w] -= c
rest[m - 2 * h + 1] -= 2 * N; rest[1] -= N
rest = {w: c for w, c in rest.items() if c}
assert all(c > 0 for c in rest.values())
for w, c in rest.items():
    for x in inner(w): add(OLD, x, c)
for r, c in cr.items():
    add(OLD, m - 2 * r, 2 * v * c)
    for x in inner(r): add(OLD, x, 2 * v * c)
add(OLD, m - 4 * h + 2, 2 * N)
for x in [h - 2] + [1] * (h + 1): add(OLD, x, 2 * N)
add(OLD, 1, N)
OLDc = dict(h=h, m=m, W=W, s=s, hist={w: c for w, c in OLD.items() if c})
a_old, _ = moment.certify(OLDc)
rep('A5. old compile rebuilt from the same steps reproduces round-7.1 plain a_b = %s' % a_old, a_old == Q(12899, 200000000) and sum(w * c for w, c in OLD.items()) == s)
NEWc = dict(h=h, m=m, W=W, s=s, hist=H)
astar, root = moment.certify(NEWc)
rep('A6. moment.certify (their code) on the one-run histogram: a* = %s (%.6e), root %.6e' % (astar, float(astar), root), astar == Q(25367, 200000000))
# independent rigorous bound: Decimal (correctly rounded ln/exp) at 60 digits with explicit slack
getcontext().prec = 60
def F_dec(a):
    a = Decimal(a.numerator) / Decimal(a.denominator); tot = Decimal(0); md = Decimal(m)
    for w, c in H.items():
        x = (a * (md / Decimal(w)).ln())
        tot += Decimal(c) * Decimal(w) / md * x.exp()
    return tot / Decimal(W)
slack = Decimal(10) ** -45                        # >> accumulated rounding of ~100 correctly rounded ops at 60 digits
Fa = F_dec(astar); Fn = F_dec(astar + Q(1, 10**9))
rep('A7. independent Decimal bound: F(a*) = 1 - %.3e < 1 rigorously; F(a* + 1e-9) = 1 + %.3e (a* is the 1e-9 floor)' % (1 - Fa, Fn - 1), Fa + slack < 1 and Fn > 1)
mu1 = Q(s, W * m)
rep('A8. mu_1 = s/(W m) = 1 - %.3e < 1 (linear toll absorbable)' % float(1 - mu1), mu1 < 1)
# ---- stopped composition
TH = Q(1, 1000)
ab = (1 - TH) * astar + TH * a_old
rep('S1. stopped a_b = (1-theta) a* + theta a_old = %s = %.10e' % (ab, float(ab)), ab == Q(6338633, 5 * 10**10))
rep('S2. theta = 1/1000 > a_b (adapter toll e^(1-theta) subordinate to e^(1-a_b))', TH > ab)
rep('S3. a_b > 2^-13 = %.7e' % (2 ** -13), ab > Q(1, 2**13))
# exponent check: cost exponent = max(1-theta, (1-theta)(1-a*) + theta(1-a_old)) = 1 - a_b
e1, e2 = 1 - TH, (1 - TH) * (1 - astar) + TH * (1 - a_old)
rep('S4. exponent identity: (1-theta)(1-a*) + theta(1-a_old) == 1 - a_b, and > 1 - theta', e2 == 1 - ab and e2 > e1)
# best theta just above a*: any theta in (a_b(theta), 1) works; a_b(theta) decreasing in theta
th2 = Q(127, 10**6); ab2 = (1 - th2) * astar + th2 * a_old
rep('S5. (side) smaller theta = 127e-6 > a_b(theta) gives a_b = %.10e (bit not binding, kappa unaffected)' % float(ab2), th2 > ab2)
# ---- kappa
A_C = Q(36926111, 5 * 10**11)
k = evaluate(ab, A_C, Q(1, 1000), 'crude', m_c=576, s_c=119453132304)
print('kappa(beta=1/1000) =', k['kappa'], '= %.10e' % float(k['kappa']), 'binding', k['binding'], 'ok', k['ok'], 'bad', k['bad'], 'log2 %.4f' % math.log2(float(k['kappa'])))
rep('K1. kappa = 184432317329/(2.5e15) exactly, ok', k['ok'] and k['kappa'] == Q(184432317329, 25 * 10**14))
k0 = evaluate(Q(12899, 200000000), A_C, Q(1, 1000), 'crude', m_c=576, s_c=119453132304)
print('  reference: w2 plain a_b -> kappa %.10e (binding %s)' % (float(k0['kappa']), k0['binding']))
kc = evaluate(A_C * 10, A_C, Q(1, 1000), 'crude', m_c=576, s_c=119453132304)
rep('K2. kappa is complex-bound: raising a_b x10 does not change kappa', kc['kappa'] == k['kappa'])
for be in (Q(1, 10**4), Q(1, 10**5), Q(1, 10**6), Q(1, 10**9)):
    kb = evaluate(ab, A_C, be, 'crude', m_c=576, s_c=119453132304)
    kb0 = evaluate(Q(12899, 200000000), A_C, be, 'crude', m_c=576, s_c=119453132304)
    print('  beta=%s: kappa %.10e ok=%s bad=%s binding=%s | same beta with w2 plain a_b: %.10e' % (be, float(kb['kappa']), kb['ok'], kb['bad'], kb['binding'], float(kb0['kappa'])))
# ---- negative controls (accounting)
Hp = dict(H); Hp[1] += 1
rep('N-A1. perturbed histogram breaks the rank sum', sum(w * c for w, c in Hp.items()) != s)
Hm = dict(H); Hm[m] = 1; Hm[m - 1] = Hm.get(m - 1, 0)  # a width-m child is not a recursion
try:
    moment.certify(dict(h=h, m=m, W=W, s=s + m, hist=Hm)); bad = False
except AssertionError:
    bad = True
rep('N-A2. a width-m child (r_u = 0) is rejected by the certifier', bad)
abad = (1 - Q(1, 10**5)) * astar + Q(1, 10**5) * a_old
rep('N-A3. theta = 1e-5 < a_b: toll e^(1-theta) dominates, stopped saving collapses to theta', not (Q(1, 10**5) > abad))
kbad = evaluate(ab, A_C, Q(1), 'crude', m_c=576, s_c=119453132304)
rep('N-A4. beta = 1 is rejected by the assembly constraints', not kbad['ok'])
# ordinary bank: the generic basis in the ORDINARY orientation gives the longest cell: rank-r edge -> r singletons
# (r <= m/2) or (m-r) singletons + block 2r-m: certify that histogram, it must be far below a*
ORD = {}
for w, c in H.items():
    if w <= m // 2: add(ORD, 1, w * c)
    else: add(ORD, 1, (m - w) * c); add(ORD, 2 * w - m, c)
aord, _ = moment.certify(dict(h=h, m=m, W=W, s=s, hist=ORD))
rep('N-A5. same generic basis, ordinary bank order (longest cell) certifies only a = %.3e < a_old' % float(aord), aord < a_old)
print('ALL PASS' if all(ok.values()) else 'FAIL: %s' % [k for k, x in ok.items() if not x])
