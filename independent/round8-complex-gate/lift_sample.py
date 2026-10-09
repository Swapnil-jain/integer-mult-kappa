"""Literal 576-bit check of the lift reductions the replay uses, on children sampled from the actual word.
For sampled (stage, h-level move U -> V, column t): builds the lifted frames literally (stage 1: L (x) t;
stage 2: t^perp (x) F + t (x) L), checks nesting, residual dimension and nondegeneracy, computes the child phase
from a literal orthonormal basis of the lifted residual and the literal reference wt(V) - wt(U) (Gram solve at 576
bits), and compares both with the replay's reduced (h-level) values.  Also: lifted frame wt, wrap children
(residual last^perp + start), and the stage-two entrance t_a1^perp (x) t_a2^perp.
Usage: python3 lift_sample.py word_H.pkl SEED NMOVE1 NMOVE2 NWRAP NENT"""
import sys, random, pickle, time
import numpy as np
import gf2
from replay import Lift

def main(path, seed, n1, n2, nw, ne):
    t0 = time.time(); d = pickle.load(open(path, 'rb'))
    h = d['h']; T = d['T']; v = len(T); m = h * h; R = d['R']; ret = set(d['kk']['rout'].values())
    Fh = [1 << p for p in range(h)]; tvec = [sum(1 << p for p in t) for t in T]
    ops = [(op[0], op[1], op[2]) for op in d['ops']]; del d
    rnd = random.Random(seed); eta = rnd.getrandbits(m)
    comp = lambda F: tuple(gf2.perp_within(Fh, list(F)))
    moves = {1: set(), 2: set()}; wraps = {1: set(), 2: set()}
    for stage in (1, 2):
        seq = ops if stage == 1 else [(r, F, k) for r, F, k in reversed(ops)]
        cur = {}; start = {}
        for roles, F, k in seq:
            for r in roles:
                if r[0] == 's' and k == 'cread' and r[1] in ret: continue
                if r not in cur:
                    if r[0] == 's': cur[r] = F; start[r] = F
                    elif stage == 1: cur[r] = (tvec[T.index(r[1])],) if r[0] == 'x' else ()
                    else: cur[r] = F                       # stage-two data first frame: the entrance (checked apart)
                if cur[r] != F: moves[stage].add((cur[r], F)); cur[r] = F
        for r in start: wraps[stage].add((cur[r], start[r]))
    print('distinct h-level moves', {k: len(x) for k, x in moves.items()}, 'wraps', {k: len(x) for k, x in wraps.items()}, flush=True)
    LF = {s: Lift(h, s, eta, T, None, __import__('collections').Counter()) for s in (1, 2)}
    lab = {1: lambda F: list(F), 2: lambda F: list(comp(F))}
    Erow = [(eta >> (i * h)) & ((1 << h) - 1) for i in range(h)]
    Ecol = [sum(((eta >> (i * h + j)) & 1) << i for i in range(h)) for j in range(h)]
    Ot = {}
    def O(t):
        if t not in Ot: Ot[t] = gf2.orthonormal(gf2.perp_within(Fh, [tvec[t]]))
        return Ot[t]
    def phi(stage, t):
        rows = Erow if stage == 1 else Ecol
        return sum(gf2.dot(r, o) * (gf2.pc(o) % 4) for r in rows for o in O(t)) % 4
    def lifted(stage, t, L):
        if stage == 1: return [gf2.kron(l, tvec[t], h) for l in L]
        return [gf2.kron(u, e, h) for u in gf2.perp_within(Fh, [tvec[t]]) for e in Fh] + [gf2.kron(tvec[t], l, h) for l in L]
    def lit_phase(B):
        Ob = gf2.orthonormal(B)
        if Ob is None: return None
        return sum(gf2.dot(o, eta) * (gf2.pc(o) % 4) for o in Ob) % 4
    stats = dict(moves=0, wt=0, alt=0, wraps=0, entr=0, bad=0)
    for stage, n in ((1, n1), (2, n2)):
        lf = LF[stage]
        for (U, V) in rnd.sample(sorted(moves[stage]), min(n, len(moves[stage]))):
            c = rnd.randrange(v); Uh, Vh = lab[stage](U), lab[stage](V)
            Ul, Vl = lifted(stage, c, Uh), lifted(stage, c, Vh)
            up = gf2.contains(Vl, Ul)
            if not (up or gf2.contains(Ul, Vl)): stats['bad'] += 1; print('not nested', stage); continue
            lo, hi = (Ul, Vl) if up else (Vl, Ul)
            Rl = gf2.perp_within(hi, lo)
            if gf2.dim(Rl) != abs(len(Vh) - len(Uh)) or not gf2.nondegenerate(Rl): stats['bad'] += 1; print('residual', stage); continue
            ref_l = (gf2.wt_point(Vl, eta) - gf2.wt_point(Ul, eta)) % 4
            pl = lit_phase(Rl)
            red, rk, alt = lf.move(lf.fid(Uh), lf.fid(Vh))
            if pl is None:
                stats['alt'] += 1
                if not alt: stats['bad'] += 1; print('alt mismatch')
            else:
                pl = pl if up else (-pl) % 4
                if pl != ref_l: stats['bad'] += 1; print('literal child != literal reference', stage)
            if int(red[c]) % 4 != ref_l: stats['bad'] += 1; print('reduced != literal', stage)
            # lifted frame wt
            w_red = (int(lf.wt(lf.fid(Vh))[c]) + (phi(2, c) if stage == 2 else 0)) % 4
            if w_red != gf2.wt_point(Vl, eta): stats['bad'] += 1; print('wt reduced != literal', stage)
            stats['moves'] += 1; stats['wt'] += 1
        print('stage', stage, 'moves done', round(time.time() - t0), 's', stats, flush=True)
    full = [1 << p for p in range(m)]; Fe = gf2.pc(eta) % 4
    for stage in (1, 2):
        lf = LF[stage]
        for (Lst, Sst) in rnd.sample(sorted(wraps[stage]), min(nw, len(wraps[stage]))):
            c = rnd.randrange(v); Ll, Sl = lifted(stage, c, lab[stage](Lst)), lifted(stage, c, lab[stage](Sst))
            Rw = gf2.perp_within(full, Ll) + Sl
            if gf2.dim(Rw) != m - gf2.dim(Ll) + gf2.dim(Sl) or not gf2.nondegenerate(Rw): stats['bad'] += 1; print('wrap residual'); continue
            pl = lit_phase(Rw); ref_l = (Fe - gf2.wt_point(Ll, eta) + gf2.wt_point(Sl, eta)) % 4
            ref_red = (Fe - int(lf.wt(lf.fid(lab[stage](Lst)))[c]) + int(lf.wt(lf.fid(lab[stage](Sst)))[c])) % 4
            if pl != ref_l or ref_l != ref_red: stats['bad'] += 1; print('wrap phase', pl, ref_l, ref_red)
            stats['wraps'] += 1
            print('wrap', stage, round(time.time() - t0), 's', flush=True)
    for _ in range(ne):
        a1, a2 = rnd.randrange(v), rnd.randrange(v)
        tp1, tp2 = gf2.perp_within(Fh, [tvec[a1]]), gf2.perp_within(Fh, [tvec[a2]])
        want = gf2.key([gf2.kron(u, w, h) for u in tp1 for w in tp2])
        o1, o2 = O(a1), O(a2)
        red = sum(gf2.dot(gf2.kron(x, y, h), eta) * ((gf2.pc(x) * gf2.pc(y)) % 4) for x in o1 for y in o2) % 4
        for old, new in (([gf2.kron(f, tvec[a2], h) for f in Fh], lifted(2, a1, [tvec[a2]])), ([gf2.kron(u, tvec[a2], h) for u in tp1], lifted(2, a1, []))):
            if not gf2.contains(new, old): stats['bad'] += 1; print('entrance not nested'); continue
            Rl = gf2.perp_within(new, old)
            if gf2.key(Rl) != want: stats['bad'] += 1; print('entrance residual'); continue
            pl = lit_phase(Rl); ref_l = (gf2.wt_point(new, eta) - gf2.wt_point(old, eta)) % 4
            if pl != ref_l or pl != red: stats['bad'] += 1; print('entrance phase', pl, ref_l, red)
            stats['entr'] += 1
        print('entrance', round(time.time() - t0), 's', flush=True)
    stats['PASS'] = stats['bad'] == 0; print(stats, flush=True)

if __name__ == '__main__':
    a = sys.argv; main(a[1], int(a[2]), int(a[3]), int(a[4]), int(a[5]), int(a[6]))
