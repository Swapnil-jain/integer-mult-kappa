"""Mutation controls for dead-copy recycling words: each output word must be REJECTED by check_word2_recycle.py.
  unequal   one host's control replaced by a retiring row holding a DIFFERENT value (the pivot keeps signal);
  early     one host's birth replaced by a fresh birth that starts BEFORE the pair has retired;
  claim     the claimed child histogram is changed by moving one child between widths.
Usage: python3 mutate_hosts.py WORD.json.gz OUTDIR"""
import sys, gzip, json, copy, os
J = json.load(gzip.open(sys.argv[1], 'rt')); out = sys.argv[2]
ops = J['ops']; sched = J['sched']; pos = {i: k for k, i in enumerate(sched)}
rops = {}
for i in sched:
    for r in ops[i][:2]: rops.setdefault(r, []).append(i)
hold = {}
for xx, sr in J['sources'].items(): hold[sr] = int(xx)
for i in sched: hold[ops[i][0]] = ops[i][2]
H = J['hosts']; z0 = H[0]; used = {z[k] for z in H for k in ('pivot', 'control', 'birth')}
def dump(name, W):
    with gzip.open(os.path.join(out, 'mut_%s.json.gz' % name), 'wt') as f: json.dump(W, f, separators=(',', ':'))
# unequal: any other non-host role whose final value differs and that retired before the birth
tb = pos[rops[z0['birth']][0]]
alt = next(r for r, L in rops.items() if r not in used and hold.get(r) != hold.get(z0['pivot']) and pos[L[-1]] < tb
           and r not in J['rootroles'] and r not in {int(b) for b in J['tau']})
W = copy.deepcopy(J); W['hosts'][0]['control'] = alt; dump('unequal', W)
# early: a fresh birth whose first op precedes the pair's retirement
t0 = max(pos[rops[z0['pivot']][-1]], pos[rops[z0['control']][-1]])
srcs = set(J['sources'].values()); gauged = {int(b) for b in J['tau']}; sinks = {z['role'] for z in J.get('sinks', [])}
eb = next(r for r, L in rops.items() if r not in used and r not in srcs and r not in gauged and r not in sinks
          and ops[L[0]][0] == r and pos[L[0]] < t0)
W = copy.deepcopy(J); W['hosts'][0]['birth'] = eb; dump('early', W)
# claim
W = copy.deepcopy(J); C = W['claim_hosts']['C']; k1 = '1'; k2 = '2'
C[k1] = str(int(C[k1]) - 2); C[k2] = str(int(C[k2]) + 1); dump('claim', W)
print('wrote mutations unequal (control %d), early (birth %d), claim' % (alt, eb))
