"""Build the reclaimed complex word exactly as independent/round8-coreshare/complex_hist.py does, and pickle
plain data for the independent gate: ops with explicit updates, producer supports for the walk, c7.walk histogram,
and the construction's claimed shared histograms (share.shared, modes merge / own) for the given group sizes.
This file is the only place construction code is imported; the checkers read the pickle only.
Usage: python3 build_word.py H CFG GROUPSPEC OUT.pkl      (GROUPSPEC: 'pr128' or a comma list of sizes)"""
import sys, os, pickle, time, json, resource
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'round8-coreshare')); sys.path.insert(1, os.path.join(HERE, '..', 'round8-coreshare', 'lib'))
sys.setrecursionlimit(100000)
import c7, reclaim as RC, share
h = int(sys.argv[1]); cfg = sys.argv[2]; gs = sys.argv[3]; out = sys.argv[4]
fl = set(cfg.split(','))
t0 = time.time()
RC.DEEP = int(os.environ.get('DEEP', '12'))
B, Wu = RC.base(h, cfg)
rcx = RC.reclaim(B, Wu, lazy=int(os.environ.get('LAZY', '3000')), rels=('dup', 'diff', 'sum', 'deep'), verbose=False)
B = RC.assemble(B, Wu, rcx, h_targets=os.environ.get('HT', '1') == '1')   # frame-0 readouts declare all targets (as e2e_share.py)
sig = RC.check_signals(B)
if 'defer' in fl: RC.gfix(B, verbose=False)
w = c7.walk(B, check=True)
print('built', dict(R=B['R'], W=w['W'], s=w['s'], sum_ok=w['sum_ok'], bad=w['bad'], signals=sig), round(time.time() - t0), 's', flush=True)
if gs == 'pr128':
    G = json.load(open(os.path.join(HERE, '..', '..', 'certificates', 'round8', 'pr128_partition_mrp24.json')))['groups']
    sizes = [len(g) for g in G]
else: sizes = [int(x) for x in gs.split(',')]
shared = {}
for mode in ('own',):   # merge mode is not part of the release (its exterior would sit before the op's cleanup)
    try: sh = share.shared(B, w, groups=sizes, mode=mode); shared[mode] = sh
    except AssertionError as e: shared[mode] = dict(error=str(e))
c = B['c']
nodes = set(n for _, n in c.retained) | set(p[1] for p in c.pieces)
keep = dict(h=h, cfg=cfg, ops=[list(o) for o in B['ops']], T=B['T'], R=B['R'], m=B['m'], v=B['v'],
            retslot=B['retslot'], pout={i: B['rmap'][q] for i, q in B['kk']['pout'].items()},
            pieces=list(c.pieces), retained=list(c.retained), sup={n: c.sup[n] for n in nodes}, tid=dict(c.tid),
            sigma=dict(B.get('sigma', {})), walkH=w['H'], walkW=w['W'], walks=w['s'], sizes=sizes, shared=shared,
            construction_signals=sig)
pickle.dump(keep, open(out, 'wb'), protocol=5)
for mode in shared: print(mode, {k: v for k, v in shared[mode].items() if k != 'H'}, flush=True)
print('rss_MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // (1024 if sys.platform != 'darwin' else 1 << 20), 'secs', round(time.time() - t0))
