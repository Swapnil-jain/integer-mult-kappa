"""Decoder identity (H + K) x = x mod 2^61-1 straight from a built graph (no frames, no schedule): the sum over roots of
coefficient x DAG value, plus #144's cube-local K (antipodal +1/2, neighbour -1/2), must equal x on random x.
Usage: python3 ident.py p CFGJSON   (a make_word config as JSON)"""
import sys, json, random
from fractions import Fraction as Q
from make_word import build_cfg
P = (1 << 61) - 1
md = lambda q: Q(q).numerator % P * pow(Q(q).denominator % P, P - 2, P) % P
def ident(p, cfg, seed=1):
    g = build_cfg(p, cfg); args = g['args']; signs = g['signs']; inputs = g['inputs']; v = len(inputs)
    off = 1 if args[0] is None and len(args) > v and args[v] is not None else 0   # zero-based graph (inputs 0..v-1)
    rng = random.Random(seed); x = [rng.randrange(P) for _ in range(v)]; val = [0] * len(args)
    for k in range(len(args)):
        val[k] = x[k] if args[k] is None else (val[args[k][0]] + signs[k] * val[args[k][1]]) % P
    y = [0] * v
    for r in g['roots']:
        for t, c in zip(r['targets'], r['coefficients']): y[t] = (y[t] + md(c) * val[r['node']]) % P
    bad = 0
    for t in range(v):
        k0 = t - t % 8; acc = sum(x[s] if bin(inputs[t] ^ inputs[s]).count('1') == 6 else -x[s] if bin(inputs[t] ^ inputs[s]).count('1') == 2 else 0 for s in range(k0, k0 + 8))
        bad += (y[t] + acc * md(Q(1, 2)) - x[t]) % P != 0
    return bad, v, len(args) - v
if __name__ == '__main__':
    p = int(sys.argv[1]); cfg = json.loads(sys.argv[2])
    for seed in (1, 2):
        bad, v, adds = ident(p, cfg, seed); print('p=%d cfg=%s seed %d: %d wrong of %d (additions %d)' % (p, json.dumps(cfg), seed, bad, v, adds))
