"""Independent exact Q(i) dirty-scratch replay of a complete c7 word (pickled by build_word.py), both stages.

Model (the two-stage complex network): at a dual point eta in F2^m (m = h^2) every role carries a frame (a subspace
of F2^m); a frame change U -> V (nested) is one whole-residual child multiplying the role by i^{wt(p_R eta)} (going
up) or its inverse (going down), R the residual, computed from an explicit orthonormal basis of R; it is asserted
equal to the reference ratio i^{wt(p_V eta) - wt(p_U eta)} (Gram solve), except for alternating residuals, whose
phase is the reference (PR #24 normal form; counted).  Stage one, for every fixed second triple a2, runs the word on
frames L (x) t_a2 over the pairs (t, a2); stage two, for every fixed first triple a1, runs the inverse word (ops
reversed, updates inverted, labels complemented) on frames t_a1^perp (x) F + t_a1 (x) L over the pairs (a1, t).
Each (stage, fixed) instance has its own auxiliaries, started dirty (random) at their first frame, restored by the
wrap child (residual last^perp + start).  All v fixed instances run at once as numpy columns.

Exactness: values are Gaussian integers in int64 (inputs are multiples of 42^2, so every division by 2 or 42 in the
word is exact; each division is asserted exact, and every operand is asserted below 2^50 so no int64 overflow).
Lifted quantities are computed from h-level data through the identities (b (x) t).eta = b.(E t), popc(b (x) t) =
3 popc(b), Gram(b (x) t) = Gram(b) for weight-3 t; these are checked literally at 576 bits by lift_sample.py.

Target: (x, y) -> X_out = F y, Y_out = F x on every pair, with the endpoint (x' = Z_w x, B' = B + P^-1 A,
Y_out = i^9 Z_w B', X_out = -A); every auxiliary of every instance restored to F r0; histogram == c7.walk's.
Usage: python3 replay.py word_H.pkl SEED [control]"""
import sys, os, time, pickle, random
from collections import Counter, defaultdict
from fractions import Fraction as Q
import numpy as np
import gf2

COS = np.array([1, 0, -1, 0], np.int64); SIN = np.array([0, 1, 0, -1], np.int64)
LIM = 1 << 50


class Fail(Exception): pass


def rot(re, im, ph):
    k = np.asarray(ph) % 4; c = COS[k]; s = SIN[k]
    return c * re - s * im, s * re + c * im


def chk(*arrs):
    for a in arrs:
        if a.size and int(np.abs(a).max()) >= LIM: raise Fail('magnitude bound exceeded')


def exdiv(a, d):
    if (a % d).any(): raise Fail('inexact division by %d' % d)
    return a // d


class Lift:
    """h-level frames with per-column lifted phases for one stage."""
    def __init__(s, h, stage, eta, T, ctl, stats):
        s.h = h; s.stage = stage; s.ctl = ctl; s.st = stats
        n = h; E = [[(eta >> (i * h + j)) & 1 for j in range(h)] for i in range(h)]
        tv = [sum(1 << p for p in t) for t in T]
        # column c is the fixed triple T[c]; stage 1 contracts eta with it on the right, stage 2 on the left
        if stage == 1: col = [sum((sum(E[i][j] for j in range(h) if t >> j & 1) & 1) << i for i in range(h)) for t in tv]
        else: col = [sum((sum(E[i][j] for i in range(h) if t >> i & 1) & 1) << j for j in range(h)) for t in tv]
        s.eta = np.array(col, np.uint32)
        s.ids = {}; s.basis = []; s.wtc = {}; s.mv = {}
    def fid(s, vs):
        k = gf2.key(vs)
        if k not in s.ids: s.ids[k] = len(s.basis); s.basis.append(list(k))
        return s.ids[k]
    def par(s, b): return (np.bitwise_count(np.uint32(b) & s.eta) & 1).astype(np.int64)
    def wt(s, f):
        """lifted wt(p_{L-part} eta) per column, L = basis[f] (stage 2: the t (x) L part; D0 handled separately)."""
        if f in s.wtc: return s.wtc[f]
        B = s.basis[f]
        if not B: r = np.zeros(len(s.eta), np.int64)
        else:
            if not gf2.nondegenerate(B): raise Fail('degenerate frame')
            Gi = gf2.gram_inverse(B); Y = [s.par(b) for b in B]; P = np.zeros(len(s.eta), np.uint32)
            for i, b in enumerate(B):
                c = np.zeros(len(s.eta), np.int64)
                for j in range(len(B)):
                    if Gi[i] >> j & 1: c ^= Y[j]
                P ^= np.where(c == 1, np.uint32(b), np.uint32(0))
            r = (3 * (np.bitwise_count(P).astype(np.int64) % 4)) % 4
        r = r.astype(np.int8); s.wtc[f] = r; return r
    def space_phase(s, B):
        """lifted wt(p_{R (x) t} eta) for a nondegenerate h-level span R, from an orthonormal basis when R has one."""
        O = None if s.ctl == 'badbasis' else gf2.orthonormal(B)
        if s.ctl == 'badbasis': O = list(gf2.reduce_basis(B).values())
        if O is None:                                    # alternating residual: PR #24 normal form, phase = reference
            s.st['alt_distinct'] += 1; return None
        ph = np.zeros(len(s.eta), np.int64)
        for o in O: ph += s.par(o) * ((3 * (gf2.pc(o) % 4)) % 4)
        return ph % 4
    def move(s, u, w):
        """(phase per column, rank, alternating?) of the child u -> w."""
        k = (u, w)
        if k in s.mv: return s.mv[k]
        U, W = s.basis[u], s.basis[w]
        if gf2.contains(W, U): lo, hi, up = U, W, True
        elif gf2.contains(U, W): lo, hi, up = W, U, False
        else: raise Fail('frames not nested')
        R = gf2.perp_within(hi, lo)
        if len(R) != len(hi) - len(lo) or not gf2.nondegenerate(R): raise Fail('residual degenerate')
        ref = (s.wt(w) - s.wt(u)) % 4
        ph = s.space_phase(R); alt = ph is None
        if alt: ph = ref
        else:
            if not up: ph = (-ph) % 4
            if s.ctl != 'badbasis' and (ph != ref).any(): raise Fail('child phase != reference ratio')
        s.st['moves_distinct'] += 1
        s.mv[k] = (ph.astype(np.int8), len(hi) - len(lo), alt); return s.mv[k]


def run(path, seed, ctl=None):
    t0 = time.time()
    d = pickle.load(open(path, 'rb'))
    h = d['h']; T = d['T']; v = len(T); m = h * h; N = v * v; R = d['R']; args = d['args']; kk = d['kk']
    tid = {t: i for i, t in enumerate(T)}; Fh = [1 << p for p in range(h)]
    tvec = [sum(1 << p for p in t) for t in T]
    rnd = random.Random(seed); eta = rnd.getrandbits(m)
    stats = Counter(); hist = Counter()
    gate_of = {n: (ins, outs) for n, ins, outs in kk['gates']}
    rout = kk['rout']; names = sorted(rout, key=lambda nm: nm[1])
    assert d['allE'] and [nm[0] for nm in names] == ['E'] * h
    ret_slots = [rout[nm] for nm in names]
    # ---------- forward instruction list (word semantics) ----------
    fwd = []; p2 = []; after_high = False
    for op in d['ops']:
        roles, F, kind = op[:3]
        ups = []
        if kind == 'gate':
            n, excl = op[3]; ins, outs = gate_of[n]
            ups = ([(('s', ins[0]), ('s', ins[1]), Q(1))] if args[n] else []) + \
                  [(('s', o), ('s', outs[0]), Q(1)) for o in outs[1:] if o not in excl]
            p2 += ups
        elif kind == 'latecopy': ups = [(roles[1], roles[0], Q(1))]; p2 += ups
        elif kind in ('xcopy', 'vgate'): ups = [(roles[1], roles[0], Q(-1 if after_high else 1))]
        elif kind == 'inject':
            if len(op) > 3 and op[3] is not None: ups = [(roles[0], roles[1], Q(d['pieces'][op[3]][2], 2))]
        elif kind == 'high': ups = [(a, b, -c) for a, b, c in reversed(p2)]; after_high = True
        elif kind in ('read', 'cread'): ups = None
        else: raise Fail('unknown op ' + kind)
        fwd.append([list(roles), F, kind, ups])
    stats['ops'] = len(fwd)
    # scatter coefficients x42: y_S += sum_i (1/(h-3) - [i in S]/2) E_i
    den = h - 3 if ctl != 'scat22' else h - 2
    scat = np.array([[Q(1, den) - (Q(1, 2) if i in S else 0) for i in range(h)] for S in T], dtype=object)
    L42 = 2 * (h - 3) if ctl != 'scat22' else 2 * (h - 2) * (h - 3)
    scat42 = np.array([[int(x * L42) for x in row] for row in scat], np.int64)
    assert all(Q(int(x * L42), 1) == x * L42 for row in scat for x in row)
    # ---------- garbage map G (pass two at x = 0), by the adjoint, scaled by L42 ----------
    ih = next(i for i, w in enumerate(fwd) if w[2] == 'high')
    adj = np.zeros((R, v), np.int64); adjy = np.eye(v, dtype=np.int64) * L42
    for i in range(ih - 1, -1, -1):
        roles, F, kind, ups = fwd[i]
        if kind == 'read': continue
        if kind == 'cread':
            # y_S += scat[S,i] E_i  =>  adj(E_i) += sum_S scat[S,i] adj(y_S)
            for j, q in enumerate(ret_slots): adj[q] += exdiv(scat42[:, j] @ adjy, L42)
            continue
        for dst, src, c in reversed(ups):
            if src[0] == 'x': continue
            a_d = adj[dst[1]] if dst[0] == 's' else adjy[tid[dst[1]]]
            inc = a_d * c.numerator
            if c.denominator != 1: inc = exdiv(inc, c.denominator)
            if src[0] == 's': adj[src[1]] += inc
            else: adjy[tid[src[1]]] += inc
    assert int(np.abs(adj).max()) < 2 ** 15
    G42 = adj.T.astype(np.int16); del adj                          # G42[T, q] = L42 * coefficient of z_q in y_T
    stats['G_nnz'] = int((G42 != 0).sum())
    # readout coefficients with an odd (non-dyadic) denominator: these need #104's odd-denominator grid (allE)
    odd = (G42 % (L42 // 2)) != 0
    stats['G_nondyadic_entries'] = int(odd.sum()); stats['G_nondyadic_slots'] = int(odd.any(axis=0).sum()); stats['G_absmax_x%d' % L42] = int(np.abs(G42).max())
    # readouts: support and purity
    written = set(); bad = Counter()
    for w in fwd:
        roles, F, kind, ups = w
        if kind == 'read':
            q = roles[0][1]; tg = np.array(sorted(tid[r[1]] for r in roles[1:]), np.int64)
            sup = np.nonzero(G42[:, q])[0]
            if q in written: bad['impure'] += 1
            if not np.isin(sup, tg).all(): bad['support'] += 1
            w[3] = (q, tg)
        elif kind == 'high': break
        elif ups:
            for dst, _, _ in ups:
                if dst[0] == 's': written.add(dst[1])
    if bad: raise Fail('readout checks %s' % dict(bad))
    stats['reads_pure_and_supported'] = 1
    if ctl == 'skipread':
        k = next(i for i, w in enumerate(fwd) if w[2] == 'read' and w[1] and G42[:, w[3][0]].any()); fwd[k][2] = 'noop'
    if ctl == 'lateread':
        k = next(i for i, w in enumerate(fwd) if w[2] == 'read' and w[1] and G42[:, w[3][0]].any()); q = fwd[k][3][0]
        j = next(i for i in range(k + 1, len(fwd)) if fwd[i][2] not in ('read', 'cread') and fwd[i][3] and any(u[0] == ('s', q) for u in fwd[i][3]))
        fwd.insert(j + 1, fwd.pop(k))
    if ctl == 'relink':                                  # one linked gate adds a wrong slot instead of its linked argument
        lk = kk['links']; g = next(u[1] for u in lk if u[0] == 'g' and args[u[1]])
        ins_g = gate_of[g][0]; other = next(q for q in range(R) if q not in ins_g)
        for w in fwd:
            if w[2] == 'gate' and w[3] and w[3][0][0] == ('s', ins_g[0]) and w[3][0][1] == ('s', ins_g[1]):
                w[3][0] = (('s', ins_g[0]), ('s', other), Q(1)); w[0] = list(w[0]) + [('s', other)]; break
        else: raise Fail('relink target not found')
    if ctl == 'pivswap':                                 # a link-reoriented gate (pivot = second argument) adds into the wrong slot
        piv = kk['piv']; g = next(n for n in sorted(piv) if piv[n] == 1 and args[n])
        ins_g = gate_of[g][0]
        for w in fwd:
            if w[2] == 'gate' and w[3] and w[3][0][:2] == (('s', ins_g[0]), ('s', ins_g[1])):
                w[3][0] = (('s', ins_g[1]), ('s', ins_g[0]), Q(1)); break
        else: raise Fail('pivswap target not found')
    print('instructions ready', round(time.time() - t0), 's', dict(stats), flush=True)
    # ---------- inputs ----------
    rng = np.random.default_rng(seed)
    S0 = (h - 3) ** 2 * 4 if ctl != 'scat22' else L42 * L42
    gi = lambda *sh: (rng.integers(-1000, 1001, sh, dtype=np.int64) * S0, rng.integers(-1000, 1001, sh, dtype=np.int64) * S0)
    xin = gi(v, v); yin = gi(v, v)                         # [a1, a2]
    Fe = gf2.pc(eta) % 4
    cpl = {}
    def comp(F):
        if F not in cpl: cpl[F] = tuple(gf2.perp_within(Fh, list(F)))
        return cpl[F]
    # D0 phases: stage 1 D0' = F (x) t^perp, stage 2 D0 = t^perp (x) F, from orthonormal bases of t^perp
    Erow = [(eta >> (i * h)) & ((1 << h) - 1) for i in range(h)]
    Ecol = [sum(((eta >> (i * h + j)) & 1) << i for i in range(h)) for j in range(h)]
    Ot = [gf2.orthonormal(gf2.perp_within(Fh, [tv])) for tv in tvec]
    def phiD(stage):
        rows = Erow if stage == 1 else Ecol
        return np.array([sum(gf2.dot(r, o) * (gf2.pc(o) % 4) for r in rows for o in O) % 4 for O in Ot], np.int64)
    # ---------- one stage ----------
    ctx = {}
    def stage_run(stage, ins, data, lf):
        """ins: instruction list; data: dict kind -> (re, im, frame ids per row, pending entrance flags)."""
        lab = (lambda F: F) if stage == 1 else comp
        PHI = phiD(stage)
        # identity check: Fe = wt(D0 part) + lifted wt(F)
        fF = lf.fid(Fh)
        if ((PHI + lf.wt(fF)) % 4 != Fe).any(): raise Fail('orthogonal split of the full space')
        are, aim = gi(R, v); r0 = (are.copy(), aim.copy())
        afr = np.full(R, -1, np.int64)
        first = {}; rset = set(ret_slots)
        for roles, F, kind, ups in ins:
            if kind == 'noop': continue
            for r in roles:
                if r[0] == 's' and r[1] not in first and not (kind == 'cread' and r[1] in rset): first[r[1]] = F
        for q in range(R): afr[q] = lf.fid(lab(first[q]))
        start = afr.copy()
        retset = set(ret_slots)
        ci = next(i for i, w in enumerate(ins) if w[2] == 'cread')
        nxt = {}
        for w in ins[ci + 1:]:
            for r in w[0]:
                if r[0] == 's' and r[1] in retset and r[1] not in nxt: nxt[r[1]] = w[1]
        def mv_rows(arr, fr, idx, f, rowcount_hist=True, pend=None, kindname=None):
            """move rows idx of arr (re, im) to frame f, grouped by current frame."""
            idx = np.asarray(idx, np.int64)
            if pend is not None and pend[idx].any():
                entr = idx[pend[idx]]; entrance(kindname, entr, f); idx = idx[~pend[idx]]
            if not len(idx): return
            cur = fr[idx]
            for u in np.unique(cur):
                if u == f: continue
                sel = idx[cur == u]; ph, rk, alt = lf.move(int(u), f)
                arr[0][sel], arr[1][sel] = rot(arr[0][sel], arr[1][sel], ph[None, :])
                hist[rk] += v * len(sel); stats['alt_children'] += v * len(sel) * alt
                fr[sel] = f
        def entrance(kindname, rows, f):
            """first stage-two move of data rows (row = a2, column = a1) from their stage-one frames."""
            re_, im_, fr_, pend_ = data[kindname]
            want = (lambda a2: [tvec[a2]]) if kindname == 'y' else (lambda a2: [])
            for a2 in rows:
                if gf2.key(lf.basis[f]) != gf2.key(want(a2)): raise Fail('stage-two entrance frame')
            ph = ctx['PSI'][:, rows].T                          # [row a2, column a1]
            # reference: wt(new) - wt(old) = PHI(a1) + wt2(L) - wt1(old)
            ref = (PHI[None, :] + lf.wt(f)[None, :] - ctx['OLDWT'][kindname][:, rows].T) % 4
            if (ph % 4 != ref).any(): raise Fail('entrance phase != reference')
            re_[rows], im_[rows] = rot(re_[rows], im_[rows], ph)
            hist[(h - 1) ** 2] += v * len(rows); fr_[rows] = f; pend_[rows] = False
        def rowsof(r):
            return ('s', r[1]) if r[0] == 's' else (r[0], tid[r[1]])
        for oi, (roles, F, kind, ups) in enumerate(ins):
            f = lf.fid(lab(F))
            if kind == 'noop': continue
            if kind == 'cread':
                copies = np.zeros((h, v), np.int64), np.zeros((h, v), np.int64)
                for j, q in enumerate(ret_slots):
                    if stage == 2:
                        mv_rows((are, aim), afr, [q], lf.fid(lab(nxt[q])))
                    u = int(afr[q])
                    if u == f: copies[0][j], copies[1][j] = are[q], aim[q]; continue
                    ph, rk, alt = lf.move(u, f); hist[rk] += v; stats['alt_children'] += v * alt
                    copies[0][j], copies[1][j] = rot(are[q], aim[q], ph)
                ys = [tid[r[1]] for r in roles if r[0] == 'y']
                Yre, Yim, Yfr, Ypd = data['y']
                mv_rows((Yre, Yim), Yfr, ys, f, pend=Ypd, kindname='y')
                sg = 1 if stage == 1 else -1
                for part, cp in ((Yre, copies[0]), (Yim, copies[1])):
                    chk(cp); part[ys] += sg * exdiv(scat42[ys] @ cp, L42)
                chk(Yre, Yim); continue
            if kind == 'read':
                q, tg = ups
                mv_rows((are, aim), afr, [q], f)
                Yre, Yim, Yfr, Ypd = data['y']
                mv_rows((Yre, Yim), Yfr, tg, f, pend=Ypd, kindname='y')
                sg = -1 if stage == 1 else 1
                pend_reads.append((q, tg, sg))
                # flush when the next instruction is not a read at the same frame (batched, commuting adds)
                nx = ins[oi + 1] if oi + 1 < len(ins) else None
                if nx is None or nx[2] != 'read' or lf.fid(lab(nx[1])) != f: flush_reads(Yre, Yim, are, aim)
                continue
            # generic op: move every role, then apply updates in order
            for r in roles:
                if r[0] == 's': mv_rows((are, aim), afr, [r[1]], f)
                else:
                    D = data[r[0]]; mv_rows((D[0], D[1]), D[2], [tid[r[1]]], f, pend=D[3], kindname=r[0])
            for dst, src, c in ups:
                def ref_(r):
                    if r[0] == 's': return are, aim, r[1]
                    D = data[r[0]]; return D[0], D[1], tid[r[1]]
                dre, dim_, di = ref_(dst); sre, sim, si = ref_(src)
                for dp, sp in ((dre, sre), (dim_, sim)):
                    inc = sp[si] * c.numerator
                    if c.denominator != 1: inc = exdiv(inc, c.denominator)
                    dp[di] = dp[di] + inc
                chk(dre[di], dim_[di])
        # wrap children
        Fid = fF
        okw = 0
        for f_last in np.unique(afr):
            for f_st in np.unique(start[afr == f_last]):
                sel = np.nonzero((afr == f_last) & (start == f_st))[0]
                Ll = lf.basis[int(f_last)]; Ls = lf.basis[int(f_st)]
                Wh = gf2.perp_within(Fh, Ll) + list(Ls)
                if gf2.dim(Wh) != h - len(Ll) + len(Ls) or not gf2.nondegenerate(Wh): raise Fail('wrap residual')
                if ctl == 'nowrapsig': Wh = gf2.perp_within(Fh, Ll)
                ph = lf.space_phase(Wh) if Wh else np.zeros(v, np.int64)
                if ph is None: ph = lf.wt(lf.fid(Wh))
                ph = (PHI + ph) % 4
                ref = (Fe - lf.wt(int(f_last)) + lf.wt(int(f_st))) % 4
                if ctl not in ('badbasis', 'nowrapsig') and (ph != ref).any(): raise Fail('wrap phase != reference')
                are[sel], aim[sel] = rot(are[sel], aim[sel], ph[None, :])
                hist[m - len(Ll) + len(Ls)] += v * len(sel)
        e_re, e_im = rot(r0[0], r0[1], Fe)
        notrest = int(((are != e_re) | (aim != e_im)).any(axis=1).sum())
        stats['aux_not_restored_s%d' % stage] = notrest
        return lf
    pend_reads = []
    def flush_reads(Yre, Yim, are, aim):
        if not pend_reads: return
        qs = np.array([q for q, _, _ in pend_reads]); sg = pend_reads[0][2]
        rows = np.unique(np.concatenate([tg for _, tg, _ in pend_reads]))
        Gs = G42[np.ix_(rows, qs)].astype(np.float64)
        for part, src in ((Yre, are), (Yim, aim)):
            Z = src[qs]; chk(Z)
            bound = float(np.abs(Gs).sum(axis=1).max()) * float(np.abs(Z).max() + 1)
            if bound >= 2.0 ** 52: raise Fail('float matmul bound')
            prod = np.rint(Gs @ Z.astype(np.float64)).astype(np.int64)
            part[rows] += sg * exdiv(prod, L42)
        chk(Yre[rows], Yim[rows]); pend_reads.clear()
    # stage one: rows = first triple a1, columns = fixed a2; X rows start at their line, Y rows at 0
    LF1 = Lift(h, 1, eta, T, ctl, stats); LF2 = Lift(h, 2, eta, T, ctl, stats)
    xfr = np.array([LF1.fid([tvec[t]]) for t in range(v)], np.int64); yfr = np.full(v, LF1.fid([]), np.int64)
    data1 = {'x': (xin[0].copy(), xin[1].copy(), xfr, np.zeros(v, bool)), 'y': (yin[0].copy(), yin[1].copy(), yfr, np.zeros(v, bool))}
    if True:
        stage_run(1, fwd, data1, LF1)
        print('stage 1 done', round(time.time() - t0), 's', dict(stats), flush=True)
        # stage-one end frames must be X: F (x) t_a2, Y: t_a1^perp (x) t_a2
        if any(gf2.key(LF1.basis[int(f)]) != gf2.key(Fh) for f in data1['x'][2]): raise Fail('stage-one X end frame')
        if any(gf2.key(LF1.basis[int(f)]) != gf2.key(gf2.perp_within(Fh, [tvec[t]])) for t, f in enumerate(data1['y'][2])):
            raise Fail('stage-one Y end frame')
        # entrance phases PSI[a1, a2] on t_a1^perp (x) t_a2^perp, and the stage-one wt of the old frames
        Em = np.array([[(eta >> (i * h + j)) & 1 for j in range(h)] for i in range(h)], np.float64)
        Om = [np.array([[(o >> p) & 1 for p in range(h)] for o in O], np.float64) for O in Ot]
        Wv = [np.array([gf2.pc(o) % 4 for o in O], np.int64) for O in Ot]
        Qall = np.concatenate([Em @ Om[a2].T for a2 in range(v)], axis=1)          # (h, v*(h-1))
        W2 = np.stack(Wv)                                                           # (v, h-1)
        PSI = np.zeros((v, v), np.int64)
        for a1 in range(v):
            M = (np.rint(Om[a1] @ Qall).astype(np.int64) % 2).reshape(h - 1, v, h - 1)
            PSI[a1] = (Wv[a1] @ (M * W2[None, :, :]).sum(axis=2)) % 4
        # old-frame wt: X: wt1(F) at column a2; Y: wt1(t_a1^perp) at column a2
        wtF1 = LF1.wt(LF1.fid(Fh))
        ctx['PSI'] = PSI
        ctx['OLDWT'] = {'y': np.tile(wtF1, (v, 1)), 'x': np.stack([LF1.wt(LF1.fid(gf2.perp_within(Fh, [tvec[a1]]))) for a1 in range(v)])}
        # stage two: rows = second triple a2, columns = fixed a1.  word role 'y' <- X data, 'x' <- Y data
        dummy = LF2.fid([])
        data2 = {'y': (data1['x'][0].T.copy(), data1['x'][1].T.copy(), np.full(v, dummy, np.int64), np.ones(v, bool)),
                 'x': (data1['y'][0].T.copy(), data1['y'][1].T.copy(), np.full(v, dummy, np.int64), np.ones(v, bool))}
        del data1
        inv = [[w[0], w[1], w[2], ([(a, b, -c) for a, b, c in reversed(w[3])] if w[2] not in ('read', 'cread', 'noop') else w[3])]
               for w in reversed(fwd)]
        stage_run(2, inv, data2, LF2)
        print('stage 2 done', round(time.time() - t0), 's', flush=True)
        # final moves: X to the full space, Y to (t_a1 (x) t_a2)^perp; stage-two frames t_a1^perp (x) F + t_a1 (x) L
        Xre, Xim, Xfr, Xpd = data2['y']; Yre, Yim, Yfr, Ypd = data2['x']
        if Xpd.any() or Ypd.any(): raise Fail('data role untouched in stage two')
        fF2 = LF2.fid(Fh)
        for a2 in range(v):
            for (re_, im_, fr_), tgt in (((Xre, Xim, Xfr), fF2), ((Yre, Yim, Yfr), LF2.fid(gf2.perp_within(Fh, [tvec[a2]])))):
                u = int(fr_[a2])
                if u != tgt:
                    ph, rk, alt = LF2.move(u, tgt); re_[a2], im_[a2] = rot(re_[a2], im_[a2], ph); hist[rk] += v; fr_[a2] = tgt
        hist[1] += N
        A = (Xre.T, Xim.T); Bv = (Yre.T, Yim.T)                                     # back to [a1, a2]
        Tm = np.array([[(tv >> p) & 1 for p in range(h)] for tv in tvec], np.float64)
        cw = (np.rint(Tm @ Em @ Tm.T).astype(np.int64)) % 2                        # w.eta, w = t_a1 (x) t_a2
        Pe = (9 * cw) % 4
        res = {}
        eq = lambda a, b: int(((a[0] != b[0]) | (a[1] != b[1])).sum())
        add = lambda a, b: (a[0] + b[0], a[1] + b[1])
        res['A'] = eq(A, rot(*yin, Fe + 2))
        res['B'] = eq(Bv, add(rot(*xin, Fe - 2 * Pe), rot(*yin, Fe - Pe)))
        pinv = -cw if ctl != 'noPinv' else 0 * cw
        Bc = add(Bv, rot(*A, pinv)) if ctl != 'noPinvchild' else Bv
        pw = (Fe + 9 - 2 * cw) % 4                                                  # popc(eta + w) mod 4
        # literal spot check of the popcount formula
        for _ in range(500):
            a1, a2 = rnd.randrange(v), rnd.randrange(v)
            if gf2.pc(eta ^ gf2.kron(tvec[a1], tvec[a2], h)) % 4 != pw[a1, a2]: raise Fail('popc(eta+w) formula')
        res['Yout'] = eq(rot(*Bc, 9), rot(*xin, pw))
        res['Xout'] = eq(rot(*A, 2), rot(*yin, Fe))
        H = {k: n for k, n in sorted(hist.items()) if n}
        Wn = 2 * N + 2 * v * R
        s = sum(k * n for k, n in H.items())
        ok = not any(res.values()) and stats['aux_not_restored_s1'] == 0 and stats['aux_not_restored_s2'] == 0
        out = dict(h=h, seed=seed, control=ctl, PASS=ok and H == d['walkH'], wrong=res,
                   aux_not_restored=(stats['aux_not_restored_s1'], stats['aux_not_restored_s2']),
                   hist_equals_walk=(H == d['walkH']), W=Wn, s=s, s_equals_walk=(s == d['walks']),
                   moves_distinct=stats['moves_distinct'], alt_distinct=stats['alt_distinct'], alt_children=stats['alt_children'],
                   G_nnz=stats['G_nnz'], G_nondyadic_slots=stats['G_nondyadic_slots'], G_nondyadic_entries=stats['G_nondyadic_entries'], secs=round(time.time() - t0))
        pickle.dump(H, open('hist_replay_%d_%d%s.pkl' % (h, seed, '_' + ctl if ctl else ''), 'wb'))
        return out


if __name__ == '__main__':
    path = sys.argv[1]; seed = int(sys.argv[2]); ctl = sys.argv[3] if len(sys.argv) > 3 else None
    try: print(run(path, seed, ctl), flush=True)
    except Fail as e: print(dict(path=path, seed=seed, control=ctl, PASS=False, FAIL=str(e)), flush=True)
