"""Dump the round-seven bit histograms of both witnesses (staircase and round-six entrance corners) with the
round-six complex histogram, and every assembly, for the Lean check. The bit histograms are rebuilt from the frozen
schedules by scripts/certificate_round7.py; the complex side is round six's (lean/round6-histograms.json).
'bit' / 'kappa' is the headline (the last witness of certificate_round7.WITNESSES, staircase corner).
Usage (from the repository root): python3 lean/dump7.py"""
import sys, json, os
from fractions import Fraction as Q
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'scripts'))
import certificate_round7 as c7

fr = lambda q: [Q(q).numerator, Q(q).denominator]
cx = json.load(open(os.path.join(HERE, 'round6-histograms.json')))['cx']
out = dict(sides=[], kappas=[], cx=cx)
n = len(c7.WITNESSES)
for idx, wit in enumerate(c7.WITNESSES):
    _, res = c7.certify(wit)
    tag = '' if idx == n - 1 else str(idx + 1)                 # the headline witness gets the plain names
    for key, corner_label, suffix in (('staircase', 'staircase entrance corner', ''), ('plain', 'round-six entrance corner', 'plain')):
        a, k, c = res[key]
        side = 'bit' + tag + suffix; kk = 'kappa' + tag + suffix
        out['sides'].append([side, 'bit interchange, witness %d (%s), %s' % (idx + 1, wit['name'], corner_label)])
        out['kappas'].append([kk, side])
        out[side] = dict(a=fr(a), m=c['m'], W=c['W'], s=c['s'], hist=sorted(c['hist'].items()))
        out[kk] = dict(beta=[1, 1000], eps=fr(k['eps']), x=fr(k['x']), c1=fr(k['c1']), kappa=fr(k['kappa']), mc=cx['m'])
        assert sum(w * m for w, m in out[side]['hist']) == c['s'] and cx['m'] == c7.M_C and cx['s'] == c7.S_C
        print(side, 'a_b', a, 'kappa', k['kappa'], k['ok'])
out['sides'].append(['cx', 'complex interchange (round six)'])
json.dump(out, open(os.path.join(HERE, 'round7-histograms.json'), 'w'))
