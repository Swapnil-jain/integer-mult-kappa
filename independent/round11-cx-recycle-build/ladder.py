"""Lever ladder for one word config (certified a_c at each rung):
  caps      : our relaxed arcs, every node at its cap, #144 gauges, no reuse
  desc      : + node descent
  reuse144  : + birth reuse with #144's gauges (exact matching)
  pairaware : + pair-aware gauge reselection (this is final2's 'before pdescent')
(pdescent on top is final2's result.) Usage: python3 ladder.py p name"""
import sys, os
from resel_br import *
import make_word, cube_word as cw
p, name = int(sys.argv[1]), sys.argv[2]; a0 = 5.6e-4
W = make_word.make(p, name)
def rung(Wx, pairs_on, tag):
    F = gfirst_frames(Wx); tau = late_tau(Wx, F)
    pairs = match(candidates(Wx, F, tau, a0, True, True)) if pairs_on else []
    C, Wv, inf = recount(Wx, F, pairs, tau, False)
    a, _ = acert(C, Wv, Wx['m']); log('%s %s: gauges %d pairs %d R %d W %s bad %s a_c %s = %.10e' % (name, tag, len(Wx['s']['sel']), len(pairs), Wx['pr']['R'], Wv, inf['bad'], a, float(a)))
pr0, s0 = cw.select_gauges(W['prof0'], W['wit0'])
Wc = dict(W, prof=W['prof0'], wit=W['wit0'], pr=pr0, s=s0)
log('%s arcs matched %d, R %d, c %d, q %d' % (name, W['prof0']['matched'], W['prof0']['R'], W['prof0']['c'], W['prof0']['q']))
rung(Wc, False, 'caps')
rung(W, False, 'desc')
rung(W, True, 'reuse144')
F0 = gfirst_frames(W); sel, _, _ = reselect(W, F0, a0, set(), True); rung(with_sel(W, sel), True, 'pairaware')
