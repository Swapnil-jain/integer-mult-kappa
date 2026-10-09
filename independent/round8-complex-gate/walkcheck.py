"""Time-order value walk of a c7 word on its actual op sequence (stage one; stage two is its exact inverse).

Slots hold sums of inputs x_t (supports as bitsets over triples, no cancellation allowed); the scratch is zero here
(the dirty part is the replay's job). Checked, in word order:
  - every gate: its two input slots hold exactly the two argument values (links: the linked slot still holds the
    donor's argument), every output slot then holds the node value; one gate per active addition, roles == ins+outs;
  - every inject reads a slot holding its piece node; the centre read sees E_0..E_h-1 in the retained slots;
  - every readout is pure (its slot unwritten so far) and deferred readouts come after the centre read;
    phase-one gates (before the centre read) are exactly the retained cones' gates; no deferred slot is touched
    before the centre read;
  - after 'high' and the closing V gates every slot is zero again;
  - y_S = sum of pieces (coef/2) + scatter of the centres equals x_S exactly (rationals), for every S.
Usage: python3 walkcheck.py word_H.pkl"""
import sys, pickle
from fractions import Fraction as Q
from collections import Counter
import numpy as np


def main(path, ctl=None):
    d = pickle.load(open(path, 'rb'))
    if ctl == 'reorder':                       # move one phase-two gate before the gate of one of its arguments
        ops = d['ops']; pos = {op[3][0]: i for i, op in enumerate(ops) if op[2] == 'gate'}
        i = next(pos[n] for n in sorted(pos) if d['args'][n] and n not in d['cone'] and any(a in pos for a in d['args'][n]))
        n = ops[i][3][0]; j = max(pos[a] for a in d['args'][n] if a in pos); ops.insert(j, ops.pop(i))
    if ctl == 'linkslot':                      # a linked gate reads the donor's pivot slot instead of the linked slot
        kk = d['kk']; u, (x, pos_) = next((u, xp) for u, xp in kk['links'].items() if u[0] == 'g')
        gl = [list(g) for g in kk['gates']]; gx = next(g for g in gl if g[0] == x); gu = next(g for g in gl if g[0] == u[1])
        linked = gx[1][1]; gu[1] = tuple(gx[1][0] if s_ == linked else s_ for s_ in gu[1]); kk['gates'] = [tuple(g) for g in gl]
    h = d['h']; T = d['T']; v = len(T); R = d['R']; args = d['args']; kk = d['kk']; tid = {t: i for i, t in enumerate(T)}
    sup = {}
    for n in d['active']:
        if args[n] is None: sup[n] = 1 << (n - 1)
        else:
            a, b = args[n]; assert not sup[a] & sup[b], 'cancellation in producer'; sup[n] = sup[a] | sup[b]
    gate_of = {n: (ins, outs) for n, ins, outs in kk['gates']}
    rout = kk['rout']; pout = kk['pout']
    cone = set(d['cone']); sigma = d['sigma']
    val = [0] * R; bad = Counter(); seen_gate = Counter(); written = set(); after_cread = False; after_high = False
    p2 = []; Y = np.zeros((v, v), np.int64); L42 = 2 * (h - 3)
    scat = [[Q(1, h - 3) - (Q(1, 2) if i in S else 0) for i in range(h)] for S in T]
    def addv(dst, src_val):
        if val[dst] & src_val: bad['overlap'] += 1
        val[dst] |= src_val
    def bits(x):
        out = []
        while x: lo = x & -x; out.append(lo.bit_length() - 1); x ^= lo
        return out
    for op in d['ops']:
        roles, F, kind = op[:3]
        if kind == 'read':
            q = roles[0][1]
            if q in written: bad['impure read'] += 1
            if (q in sigma) != bool(F): bad['read frame vs deferral'] += 1
            if q in sigma and not after_cread: bad['deferred read before centre read'] += 1
            continue
        if kind == 'gate':
            n, excl = op[3]; ins, outs = gate_of[n]; seen_gate[n] += 1
            if set(r[1] for r in roles) != set(ins + outs) - set(excl): bad['gate roles'] += 1
            if (n in cone) == after_cread: bad['gate phase'] += 1
            if args[n]:
                got = sorted([val[ins[0]], val[ins[1]]]); want = sorted([sup[args[n][0]], sup[args[n][1]]])
                if got != want: bad['gate inputs'] += 1
                addv(ins[0], val[ins[1]]); p2.append((ins[0], ins[1])); written.add(ins[0])
            for o in outs[1:]:
                if o in excl: continue
                addv(o, val[outs[0]]); p2.append((o, outs[0])); written.add(o)
            if any(val[o] != sup[n] for o in outs if o not in excl): bad['gate outputs'] += 1
            continue
        if kind == 'vgate':
            t = tid[roles[0][1]]; q = roles[1][1]
            if not after_high:
                if not after_cread and q in sigma: bad['deferred slot touched in phase one'] += 1
                addv(q, 1 << t); written.add(q)
            else:
                if not val[q] >> t & 1: bad['vgate undo'] += 1
                val[q] ^= 1 << t
            continue
        if kind == 'inject':
            if len(op) > 3 and op[3] is not None:
                S, node, cf = d['pieces'][op[3]]; q = roles[1][1]
                if q != pout[op[3]] or val[q] != sup[node]: bad['inject value'] += 1
                for t in bits(sup[node]): Y[tid[S], t] += cf * (L42 // 2)
            continue
        if kind == 'cread':
            after_cread = True
            for i in range(h):
                q = rout[('E', i)]; node = dict(d['retained'])[('E', i)]
                if val[q] != sup[node]: bad['centre value'] += 1
                b = bits(sup[node])
                for S in T: Y[tid[S], b] += int(scat[tid[S]][i] * L42)
            continue
        if kind == 'high':
            after_high = True
            for dst, src in reversed(p2):
                if val[src] & ~val[dst]: bad['high undo'] += 1
                val[dst] ^= val[src]
            continue
        bad['unknown op ' + kind] += 1
    nadd = sum(1 for n in d['active'] if args[n])
    res = dict(h=h, gates=sum(seen_gate.values()), additions=nadd, each_gate_once=all(seen_gate[n] == 1 for n in d['active'] if args[n]),
               slots_zero_at_end=all(x == 0 for x in val), y_equals_x=bool((Y == np.eye(v, dtype=np.int64) * L42).all()),
               links=len(kk.get('links', {})), deferred=len(sigma), bad=dict(bad))
    res['PASS'] = res['each_gate_once'] and res['gates'] == nadd and res['slots_zero_at_end'] and res['y_equals_x'] and not bad
    print(res, flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
