"""Time-order value walk of a reclaimed complex word (explicit signed updates), stage one (stage two is its literal
inverse, replayed exactly by sreplay.py). Scratch is zero here; slots hold exact integer combinations of the inputs
x_t (dict t -> coefficient), so no cancellation pattern is assumed. Checked in word order:
  - every update is dst += c * src between distinct roles, c in {+1, -1} on slots (signed shears), c = piece/2 on y;
  - a +1 update adds a support disjoint from dst's, a -1 update removes a sub-sum dst actually holds (no cancellation
    through unrelated terms); counts of each;
  - at the centre read every retained slot holds exactly its node (0/1 indicator of the node's support);
  - every inject reads a slot holding exactly its piece node;
  - at 'high' every piece / retained slot holds its terminal node; after 'high' and the tail every slot is zero;
  - y_S = sum of pieces (coef/2) + scatter of the centres equals x_S exactly (rationals) for every S;
  - readouts: in stage one, every readout slot is unwritten up to its read (pure dirty value).
Usage: python3 walk2.py word.pkl [control]"""
import sys, pickle
from fractions import Fraction as Q
from collections import Counter, defaultdict


def main(path, ctl=None):
    d = pickle.load(open(path, 'rb'))
    h = d['h']; T = d['T']; v = len(T); R = d['R']; tid = {t: i for i, t in enumerate(T)}
    ops = [list(o) for o in d['ops']]
    if ctl == 'dropshear':                      # remove one signed clear: the walk must notice the leftover signal
        k = next(i for i, o in enumerate(ops) if o[2] == 'shear'); ops.pop(k)
    if ctl == 'swapgate':                       # move one gate before the previous gate writing one of its sources
        gs = [i for i, o in enumerate(ops) if o[2] == 'gate' and o[3]]
        for i in gs[1:]:
            srcs = {s for _, s, _ in ops[i][3]}
            j = max((k for k in gs if k < i and any(dd in srcs for dd, _, _ in ops[k][3])), default=None)
            if j is not None: ops.insert(j, ops.pop(i)); break
    sup = d['sup']; ind = lambda n: {b: 1 for b in range(v) if sup[n] >> b & 1}
    retslot = d['retslot']; pout = d['pout']; pieces = d['pieces']
    retnode = dict(d['retained'])
    val = defaultdict(dict); bad = Counter(); cnt = Counter(); written = set()
    Y = [defaultdict(Q) for _ in range(v)]; L = Q(1, h - 3)
    after_high = False
    def get(r):
        if r[0] == 'x': return {d['tid'][r[1]]: 1}
        return val[r]
    for op in ops:
        roles, F, kind, ups = op[:4]
        if kind == 'read':
            q = roles[0]
            if q in written: bad['impure read'] += 1
            continue
        if kind == 'cread':
            for q, nm in retslot.items():
                if val[('s', q)] != ind(retnode[nm]): bad['centre value'] += 1
            for dst, src, cf in ups:
                for b, c in val[tuple(src)].items(): Y[tid[dst[1]]][b] += Q(cf) * c
            cnt['cread'] += 1; continue
        if kind == 'high':
            for q, nm in retslot.items():
                if val[('s', q)] != ind(retnode[nm]): bad['retained terminal'] += 1
            for i, q in pout.items():
                if val[('s', q)] != ind(pieces[i][1]): bad['piece terminal'] += 1
            after_high = True
        if kind == 'inject':
            for dst, src, cf in ups or ():
                q = tuple(src)[1]; S = tuple(dst)[1]
                ok = [i for i, pq in pout.items() if pq == q and pieces[i][0] == S and Q(pieces[i][2], 2) == Q(cf)]
                if len(ok) != 1 or val[('s', q)] != ind(pieces[ok[0]][1]): bad['inject value'] += 1
                cnt['inject'] += 1
        for dst, src, cf in ups or ():
            dst = tuple(dst); src = tuple(src); cf = Q(cf)
            if dst == src: bad['self'] += 1; continue
            s_ = get(src)
            if dst[0] == 'y':
                if kind != 'inject': bad['y written outside inject/cread'] += 1
                for b, c in s_.items(): Y[tid[dst[1]]][b] += cf * c
                cnt['inject_terms'] += 1; continue
            if dst[0] != 's': bad['write to data role'] += 1; continue
            dv = val[dst]
            if cf == 1:
                if set(dv) & set(s_): bad['+1 overlap'] += 1
                for b, c in s_.items(): dv[b] = dv.get(b, 0) + c
                cnt['add'] += 1
            elif cf == -1:
                if any(dv.get(b, 0) != c for b, c in s_.items()): bad['-1 not a held sub-sum'] += 1
                for b, c in s_.items():
                    dv[b] = dv.get(b, 0) - c
                    if dv[b] == 0: del dv[b]
                cnt['sub'] += 1
            else: bad['slot coefficient'] += 1
            if not after_high: written.add(dst)
    # inject checks: every piece slot read by an inject held its node at that time is implied by the terminal check
    # (pieces are never written after their terminal value); verify directly by re-walking injects
    nz = sum(1 for r, dv in val.items() if r[0] == 's' and dv)
    yok = all({b: c for b, c in Y[tid[S]].items() if c} == {tid[S]: 1} for S in T)
    cnt['injects_expected'] = len(pieces)
    res = dict(h=h, R=R, ops=len(ops), counts=dict(cnt), slots_nonzero_at_end=nz, y_equals_x=yok, bad=dict(bad), control=ctl)
    res['PASS'] = nz == 0 and yok and not bad
    print(res, flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
