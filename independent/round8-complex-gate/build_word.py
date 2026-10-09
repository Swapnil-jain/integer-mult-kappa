"""Build the round-eight complex word (before sharing; independent/round8-coreshare/complex_hist.py builds the same
word) and pickle it for the independent replay tools in this directory (replay.py, walkcheck.py), which do not
import c7. Usage: python3 independent/round8-complex-gate/build_word.py OUT.pkl"""
import os, sys, pickle, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'round8-coreshare'))
sys.setrecursionlimit(100000)
import c7
from complex_hist import build
t0 = time.time()
B = build()
w = c7.walk(B, check=True)
assert w['sum_ok'] and not w['bad']
print({k: w[k] for k in ('W', 's', 'sum_ok', 'bad', 'L', 'maxrank')}, round(time.time() - t0), 's', flush=True)
c = B['c']
keep = dict(h=B['h'], ops=B['ops'], sigma=B['sigma'], M=B['M'], m=B['m'], v=B['v'], R=B['R'], T=B['T'], retslot=B['retslot'],
            cone=B['cone'], tg=B['tg'], vls=B['vls'], vorder=B['vorder'], kk=B['kk'], walkH=w['H'], walkW=w['W'], walks=w['s'],
            args=c.args, pieces=c.pieces, retained=c.retained, tid=c.tid, triples=c.triples, allE=c.allE,
            active=sorted(c.active), sup=None)
pickle.dump(keep, open(sys.argv[1], 'wb'), protocol=5)
print('wrote', sys.argv[1])
