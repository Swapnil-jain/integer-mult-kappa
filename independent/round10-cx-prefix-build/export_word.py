"""The ONLY release file that imports construction code. Rebuilds a final word (make_word config + frozen plan pickle
from final.py / final2.py) and writes everything the independent checker needs as one JSON:
  signed DAG (node 0 dummy, inputs 1..v), roots, ops (dest, ctl, node), schedule + cut, sources, root roles,
  gauges (role, annihilator A, targets), late read clocks tau, birth pairs (donor, recipient), the frame of every op,
  and the CLAIMED child histogram, W and a_c (the checker recomputes all three and compares).
Usage: python3 export_word.py P NAME PLAN.pkl OUT.json.gz"""
import sys, os, json, gzip, pickle
from resel_br import *
import make_word

p, name, plan, out = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
W = make_word.make(p, name); PL = pickle.load(open(plan, 'rb'))
sel, pairs, tau = PL['sel'], PL['pairs'], PL['tau']
W2 = with_sel(W, sel); F = gfirst_frames(W2)
assert late_tau(W2, F) == tau, 'plan clocks differ from the rebuilt schedule'
if 'frames' in PL: F = dict(F, frames=PL['frames'])
C, Wv, inf = recount(W2, F, pairs, tau, True); assert not inf['bad'], inf['bad']
a, _ = acert(C, Wv, W['m'])
G = W2['wit']['g']; s = W2['s']
signs = list(G['signs'])
if len(signs) == len(G['args']) - 1: signs = [1] + signs   # compile_frames shifts args, not signs
assert len(signs) == len(G['args'])
J = dict(p=p, h=W['h'], v=W['pr']['v'], m=W['m'], name=name, loss=W['pr']['loss'], R=W['pr']['R'], c=W['pr']['c'],
         inputs=G['inputs'], args=G['args'], signs=signs,
         roots=[dict(node=r['node'], targets=r['targets'], coefficients=r['coefficients'], kind=r['kind'],
                     coordinate=r.get('coordinate')) for r in G['roots']],
         ops=[list(o) for o in s['ops']], sched=F['sched'], cut=F['cut'],
         sources={str(x): sr for x, sr in s['sources'].items()}, rootroles=s['rootroles'],
         gauges=[dict(role=z['role'], A=list(z['A']), targets=z['targets']) for z in sel],
         tau={str(b): t for b, t in tau.items()}, pairs=[list(q) for q in pairs],
         frames=[list(f) for f in F['frames']],
         claim=dict(C={str(r): str(n) for r, n in sorted(C.items())}, W=str(Wv), a_c=str(a)))
with gzip.open(out, 'wt') as f: json.dump(J, f, separators=(',', ':'))
log('wrote %s: a_c %s = %.10e, W %s, pairs %d, ops %d' % (out, a, float(a), Wv, len(pairs), len(s['ops'])))
