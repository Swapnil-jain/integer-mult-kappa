"""Round-eight complex histogram: our two-stage complex word with completed-core scratch sharing (PR #128's mechanism,
reimplemented in share.py from its description; PR #128's code is not used).

The word is c7.build at h = 24 with: allE (retained centres), alt (PR #24 normal form), links (PR #104 carrier
matching), dag117 (PR #117's frozen addition DAG, certificates/round8/pr117_dag.json.gz, read as data and every
support recomputed by dagprod.py), lift (lifted frames), defer (B-defer), vleaf (V leaves), clos (phase-1 closure).
c7.walk checks the histogram and every step; share.shared regroups the cores into PR #128's 87 orthogonal triple
groups (83 of size 24, 4 of size 8) and returns the shared histogram (mode 'merge': each group's fix-up child merged
into the core's last step where legal). W m - s is unchanged by the sharing (asserted).
e2e_share.py replays the sharing exactly over Q(i); partition_check.py checks PR #128's partition.

Usage: python3 independent/round8-coreshare/complex_hist.py [OUT.json] [cfg] [mode]
       (no OUT: compare with certificates/round8/cx_hist.json.gz)"""
import gzip, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
sys.setrecursionlimit(100000)
import c7
from cert import cert
from share import shared

CFG = 'allE,alt,links,dag117,lift,defer,vleaf,clos'
DAG = os.path.join(ROOT, 'certificates', 'round8', 'pr117_dag.json.gz')


def build(cfg=CFG, h=24):
    fl = set(cfg.split(','))
    kw = dict(allE='allE' in fl, lift='lift' in fl, defer='defer' in fl, vleaf='vleaf' in fl, alt='alt' in fl,
              links='links' in fl, pasm='pasm' in fl, ivec='ivec' in fl, clos='clos' in fl,
              pstar=(4 if 'dual+' in fl else 0), dag=(DAG if 'dag117' in fl else None))
    return c7.build(h, verbose=False, **kw)


def main(argv):
    out_path = argv[0] if argv else None
    cfg = argv[1] if len(argv) > 1 else CFG
    mode = argv[2] if len(argv) > 2 else 'merge'
    t0 = time.time()
    B = build(cfg)
    w = c7.walk(B, check=True); m = B['m']
    assert w['sum_ok'] and not w['bad'], (w['sum_ok'], w['bad'])
    sh = shared(B, w, mode=mode)
    assert sh['s'] == sh['W'] * m - (w['W'] * m - w['s']), 'deficit changed'
    a_cert, _ = cert(sh['H'], sh['W'], m)                       # cert.py: a second, different moment bound
    out = dict(side='complex', source='completed-core sharing (%s), cfg %s (h=24)' % (mode, cfg),
               m=m, W=sh['W'], s=sh['s'], hist=sorted([int(r), int(n)] for r, n in sh['H'].items()))
    print('R=%d unshared W=%d s=%d; shared W=%d s=%d D=%d widths=%d max=%d; cert.py a_c=%s (%.7e); %d s' % (
        B['R'], w['W'], w['s'], sh['W'], sh['s'], sh['D'], len(sh['H']), sh['maxrank'], a_cert, float(a_cert), time.time() - t0))
    if out_path:
        json.dump(out, open(out_path, 'w')); print('wrote', out_path); return
    with gzip.open(os.path.join(ROOT, 'certificates', 'round8', 'cx_hist.json.gz'), 'rt') as f:
        fz = json.load(f)
    same = all(fz[k] == out[k] for k in ('m', 'W', 's', 'hist'))
    print('matches certificates/round8/cx_hist.json.gz:', same)
    assert same
    if 'claim' in fz: assert str(a_cert) == fz['claim']['a_c'], 'cert.py and the certificate disagree'


if __name__ == '__main__':
    main(sys.argv[1:])
