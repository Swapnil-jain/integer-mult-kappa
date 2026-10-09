"""Independent exact check of a cube word's side part: evaluate every node from (args, signs) on random integer x,
accumulate the side roots into the targets with their coefficients, and compare with H computed from the definition
H_TS = (1 - |S cap T|)/2 for S, T in different cubes, 0 inside a cube. Also checks, for every root, that the
F2 span of the node's support is orthogonal to every root target (legality of the read frame).
Usage: python3 replay_alg.py p [variant] where variant in base | pair:<json st> | merge"""
import sys, json, random
from fractions import Fraction as Q
sys.path.insert(0, '.')
import cube_word as cw
from restrict import pair
p = int(sys.argv[1]); var = sys.argv[2] if len(sys.argv) > 2 else 'base'
_orig = cw.mod_pairs
_origt = cw.mod_triples
if var.startswith('cfg:'):   # cfg:<json {"pair": st, "tri": labels, "merge": bool}>
    cfg = json.loads(var[4:])
    if cfg.get('pair'): st = cfg['pair']; cw.mod_pairs = lambda n, w: pair(st['L'], st['s'], st['r'])
    if cfg.get('tri'):
        from restrict import tri; L = cfg['tri']; cw.mod_triples = lambda pp, w: tri(L)
    if cfg.get('pairmod'):
        _M = json.load(open(cfg['pairmod'])); cw.mod_pairs = lambda n, w: _M
    if cfg.get('merge'): cw.MERGE = 0
    var = 'cfg' + ('_ctlcoef' if cfg.get('ctlcoef') else '') + ('merge' if cfg.get('merge') else '')
elif var.startswith('pair:'):
    st = json.loads(var[5:]); cw.mod_pairs = lambda n, w: pair(st['L'], st['s'], st['r'])
if 'merge' in var: cw.MERGE = 0
g = cw.build(p); cw.mod_pairs = _orig; cw.mod_triples = _origt
if 'ctlcoef' in var:   # negative control: flip one merged/side coefficient
    r = next(r for r in g['roots'] if r['kind'] == 'side' and r.get('channel') == ('merge' if 'merge' in var else 'face2'))
    r['coefficients'] = ['-' + c for c in r['coefficients']]
if 'ctlgroup' in var:  # negative control: a face0 node read by its whole cube (one illegal half)
    r = next(r for r in g['roots'] if r['kind'] == 'side' and r.get('channel') == 'face0')
    I = r['targets'][0] // 8; r['targets'] = list(range(8 * I, 8 * I + 8)); r['coefficients'] = [r['coefficients'][0]] * 8
labels = g['labels']; v = len(labels); args = g['args']; signs = g['signs']
rnd = random.Random(7); x = [Q(rnd.randint(-10**6, 10**6)) for _ in range(v)]
val = list(x) + [None] * (len(args) - v); sup = [1 << i for i in range(v)] + [0] * (len(args) - v)
for i in range(v, len(args)):
    a, b = args[i]; val[i] = val[a] + signs[i] * val[b]; sup[i] = sup[a] | sup[b]
y = [Q(0)] * v; bad_legal = 0
cube = [i // 8 for i in range(v)]; vec = g['inputs']
def dot(a, b): return bin(a & b).count('1') & 1
for r in g['roots']:
    if r['kind'] != 'side': continue
    node = r['node']; srcs = [i for i in range(v) if sup[node] >> i & 1]
    for t, c in zip(r['targets'], r['coefficients']):
        y[t] += Q(c) * val[node]
        bad_legal += any(dot(vec[s], vec[t]) for s in srcs)
want = [Q(0)] * v
for t in range(v):
    T = set(labels[t])
    for s in range(v):
        if cube[s] != cube[t]: want[t] += Q(1 - len(T & set(labels[s])), 2) * x[s]
wrong = sum(a != b for a, b in zip(y, want))
print('p', p, 'variant', var[:40], 'targets', v, 'wrong', wrong, 'illegal reads', bad_legal, flush=True)
