"""Finite semantic guard for a frozen complex word (standing gate; standard library only, no construction code).
Reads the word JSON written by cube-prefix's export_word.py. The child histogram and W are the word's claim, which
check_word.py recounts from the role chains and compares; run that first. Applies #144's guard formulas (the same as
cx-gate13/gate/certcheck.py and certificate_round9.finite_bridge) at this word's m, with V = |O(m,2)|, and checks
that every root coefficient denominator divides 6 (the inherited 2^-P 3^-K precision grid).
Usage: python3 cxguard.py WORD.json.gz"""
import sys, json, gzip, math
from fractions import Fraction as Q

J = json.load(gzip.open(sys.argv[1], 'rt'))
h, v, m, R, c = J['h'], J['v'], J['m'], J['R'], J['c']
ops = len(J['ops'])
hist = {int(r): Q(n) for r, n in J['claim']['C'].items() if Q(n)}; W = Q(J['claim']['W'])
s = sum(r * n for r, n in hist.items())
OK = {}
def check(name, ok): OK[name] = bool(ok); print(('PASS ' if ok else 'FAIL ') + name, flush=True)
check('m = 3h = %d, deficit W m - s = %s == 2v - 3 loss = %d' % (m, W * m - s, 2 * v - 3 * J['loss']),
      m == 3 * h and W * m - s == 2 * v - 3 * J['loss'] and W * m - s > 0)
n2 = m // 2
Sp = 2 ** ((n2 - 1) ** 2) * math.prod(4 ** i - 1 for i in range(1, n2))
V = 2 ** (m - 1) * Sp
lcm = 1
for x in list(hist.values()) + [W]: lcm = math.lcm(lcm, x.denominator)
wc, sc, vv = int(W * lcm), int(s * lcm), v * lcm
Wf, sf, Nf = V * wc, V * sc, V * vv
L = 4 * (c + v) + 10 * v + 4 * h * v + 4 * h * h + 8 * h + 8 + 2 * h + 8 * R * v * (ops + 16) + 32 * v
K = 3 * V * L + 8 * Wf + 4 * Nf + 8 * m * R * V; G = 64 * (m + 1) ** 3 * (K + 1) * (Wf + 1) ** 2
E = 64 * (Wf + m + G + 1) ** 3; B = sf + E; C0 = 32 * m * B * B; mx = max(hist)
check('guard at m=%d (V %d bits, lcm %d, c %d, ops %d, R %d, max child %d): 2GW^2+8s+4W+4+32m < E' % (
    m, V.bit_length(), lcm, c, ops, R, mx), 2 * G * Wf * Wf + 8 * sf + 4 * Wf + 4 + 32 * m < E)
check('guard: 2B(m - max child) >= s + E', 2 * B * (m - mx) >= sf + E)
check('guard: 2B + 18 < C0', 2 * B + 18 < C0)
dens = sorted({Q(x).denominator for r in J['roots'] for x in r['coefficients']})
check('root coefficient denominators %s divide 6' % dens, all(6 % d == 0 for d in dens))
print('RESULT cxguard %s p=%d: %s' % (J['name'], J['p'], 'ALL PASS' if all(OK.values()) else 'FAIL'))
sys.exit(0 if all(OK.values()) else 1)
