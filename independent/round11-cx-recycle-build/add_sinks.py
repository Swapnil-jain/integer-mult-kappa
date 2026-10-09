"""Adapter step: put the terminal-output deletions chosen by t176.py (audit-176's selection: disjoint, individually
frame-legal cube-root sinks) into a frozen word, and replace its claim by the recount after deletion, certified on
the 1e-12 grid. Usage: python3 add_sinks.py WORD.json.gz T176_OUT.json OUT.json.gz"""
import sys, json, gzip
from fractions import Fraction as Q
from evalc import acert
J = json.load(gzip.open(sys.argv[1], 'rt')); T = json.load(open(sys.argv[2]))
assert T['word'] == J['name'] and T['p'] == J['p'] and all(T['ok'].values()), 't176 run did not pass'
J['sinks'] = [dict(role=z['role'], targets=z['targets'], pivot=z['pivot']) for z in T['sinks']]
C = {int(k): Q(c) for k, c in T['C1'].items() if int(c)}; W = Q(T['W1'])
a, _ = acert(C, W, J['m'])
J['claim'] = dict(C={str(k): str(c) for k, c in sorted(C.items())}, W=str(W), a_c=str(a))
with gzip.open(sys.argv[3], 'wt') as f: json.dump(J, f, separators=(',', ':'))
print('wrote %s: %d sinks, W %s, a_c %s = %.10e' % (sys.argv[3], len(J['sinks']), W, a, float(a)))
