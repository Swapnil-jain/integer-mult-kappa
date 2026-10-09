"""Round-eight complex histogram: our two-stage complex word on eumemic's frozen h = 24 addition DAG (PR #117, read as
data, every support recomputed by lib/dagprod.py), with signed dependence reclamation (James Chang's mechanism, PR
#112, reimplemented in reclaim.py) and completed-core sharing (Andrey Mas's mechanism, PR #128, reimplemented in
share.py). Neither upstream's code is used.

The word: c7.build with allE (retained centres), alt (the PR #24 normal form for alternating residuals), links
(PR #104 carrier matching) and dag117; no lifted frames, B-defer or V leaves. reclaim.py then clears dead slots by
exact signed shears (relations dup, diff, sum and deep, depth 6, lazy window 3000) and reuses them. Frame-0
readouts declare all their targets (h_targets=True). c7.walk checks the histogram and every step.

Sharing: the 2024 triples are split into the 87 binary-orthonormal groups of certificates/round8/pr128_partition_mrp24.json
(83 of size 24, 4 of size 8). The cores of a group run consecutively on one bank of auxiliaries with dirty scratch;
each (stage, group, stream) gets ONE separate fix-up child (mode 'own'). The fused variant ('merge') is not used: its
exterior would act before the op's own cleanup. W m - s is unchanged by the sharing (asserted).

independent/round8-complex-gate/ replays the shared word exactly over Q(i) without importing this code.

Usage: python3 independent/round8-coreshare/complex_hist.py [OUT.json] [--gate HIST.pkl]
       no OUT: compare with certificates/round8/cx_hist.json.gz; --gate: also compare with a sreplay.py histogram"""
import gzip, json, os, pickle, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path[:0] = [HERE, os.path.join(HERE, 'lib')]
sys.setrecursionlimit(100000)
import c7
import reclaim as RC
from cert import cert
from share import shared

H_POINTS = 24
CFG = 'allE,alt,links,dag117'
DEEP, LAZY = 6, 3000
PARTITION = os.path.join(ROOT, 'certificates', 'round8', 'pr128_partition_mrp24.json')
FROZEN = os.path.join(ROOT, 'certificates', 'round8', 'cx_hist.json.gz')


def build():
    RC.DEEP = DEEP
    B, Wu = RC.base(H_POINTS, CFG)
    rcx = RC.reclaim(B, Wu, lazy=LAZY, rels=('dup', 'diff', 'sum', 'deep'), verbose=False)
    # frame-0 readouts must declare all their targets: without it some readouts miss targets (50 fail at h = 8)
    B = RC.assemble(B, Wu, rcx, h_targets=True)
    sig = RC.check_signals(B)
    assert sig['terminal_bad'] == 0, sig
    return B, sig


def groups():
    with open(PARTITION) as f:
        return [len(g) for g in json.load(f)['groups']]


def main(argv):
    gate = None
    if '--gate' in argv:
        i = argv.index('--gate'); gate = argv[i + 1]; argv = argv[:i] + argv[i + 2:]
    out_path = argv[0] if argv else None
    t0 = time.time()
    B, sig = build()
    w = c7.walk(B, check=True); m = B['m']
    assert w['sum_ok'] and not w['bad'], (w['sum_ok'], w['bad'])
    sizes = groups(); assert sum(sizes) == B['v'] and len(sizes) == 87
    sh = shared(B, w, groups=sizes, mode='own')
    assert sh['s'] == sh['W'] * m - (w['W'] * m - w['s']), 'deficit changed'
    assert sh['W'] == 2 * B['v'] ** 2 + 2 * len(sizes) * B['R']
    a_cert, _ = cert(sh['H'], sh['W'], m)                       # cert.py: an independent moment bound
    out = dict(side='complex', source='signed reclaim (depth %d, lazy %d) + completed-core sharing (own), cfg %s (h=%d)' % (
        DEEP, LAZY, CFG, H_POINTS), m=m, W=sh['W'], s=sh['s'], hist=sorted([int(r), int(n)] for r, n in sh['H'].items()))
    print('R=%d unshared W=%d s=%d; shared W=%d s=%d D=%d widths=%d max=%d; cert.py a_c=%s (%.9e); %d s' % (
        B['R'], w['W'], w['s'], sh['W'], sh['s'], sh['D'], len(sh['H']), sh['maxrank'], a_cert, float(a_cert), time.time() - t0))
    if gate:
        g = pickle.load(open(gate, 'rb'))
        same = {int(k): int(n) for k, n in g['H'].items() if n} == {int(k): int(n) for k, n in sh['H'].items()} and g['W'] == sh['W']
        print('matches the gate histogram %s: %s' % (os.path.basename(gate), same)); assert same
    if out_path:
        json.dump(out, open(out_path, 'w')); print('wrote', out_path); return
    with gzip.open(FROZEN, 'rt') as f:
        fz = json.load(f)
    same = all(fz[k] == out[k] for k in ('m', 'W', 's', 'hist'))
    print('matches certificates/round8/cx_hist.json.gz:', same)
    assert same


if __name__ == '__main__':
    main(sys.argv[1:])
