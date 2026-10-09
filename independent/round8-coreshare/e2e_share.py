"""Exact Q(i) two-stage replay of completed-core scratch sharing on our complex words (e2e7's model, unchanged per
core): the triples are partitioned into binary-orthonormal groups (Gram I, even pairwise overlaps); per stage, the
cores of a group run CONSECUTIVELY on ONE bank of R auxiliaries holding arbitrary Gaussian-rational scratch. Each
core relabels its incoming scratch at its own first frames (no child: a frame on pure scratch is a label) and leaves
its residual. After the group, ONE fix-up child on residual Rfix = (sum_b R_b)^perp, R_b = last_b minus start_b, is
applied (sum must be an orthogonal direct sum; child phase asserted equal to the reference ratio). PASS = every
data pair (x, y) -> (F y, F x) (endpoint correction as e2e7), every auxiliary of every bank ends at C_full r0, and the
replay histogram equals share.shared(..., mode='own') for the same group sizes.
Merge check: for every mergeable auxiliary (its last op is the step into F) the merged child on (Rfix + last-step
residual) has phase = sum of the two (so mode 'merge' has the same action).
Controls: 'nofix' (skip the fix-up), 'badgroup' (a group with a non-orthogonal pair), 'wrapfix' (fix-up residual
missing the start parts, i.e. the old per-core wrap without sigma).
Usage: python3 e2e_share.py H cfg NPTS [control] [reclaim]"""
import sys, os, random
from fractions import Fraction as Q
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
RCM = 'reclaim' in sys.argv[5:] or 'reclaim' in sys.argv[4:]
if RCM: sys.path.insert(0, os.path.join(HERE, 'rc')); sys.path.insert(1, os.path.join(HERE, 'rc', 'lib'))
else: sys.path.insert(0, HERE)
import c7, e2e7
from e2e7 import Frames, canon, kron, rot, gadd, popc
from frames import perp_in, dot
import share

pc = lambda x: bin(x).count('1')


def groups_of(T, h, control=None, seed=3):
    """greedy partition of the triples into orthonormal groups (pairwise intersections of size 0 or 2)."""
    rnd = random.Random(seed); order = list(T); rnd.shuffle(order)
    if control == 'single': return [[t] for t in order]
    tv = lambda t: sum(1 << p for p in t)
    G = []
    for t in order:
        for g in G:
            if len(g) < h and all(pc(tv(t) & tv(u)) in (0, 2) for u in g):
                # rank check: the group must stay linearly independent (orthonormal => automatic)
                g.append(t); break
        else: G.append([t])
    if control == 'badgroup':                        # merge two groups with a size-1 overlap pair
        for i in range(len(G)):
            for j in range(i + 1, len(G)):
                if any(pc(tv(a) & tv(b)) == 1 for a in G[i] for b in G[j]) and len(G[i]) + len(G[j]) <= h:
                    G[i] += G.pop(j); return G
    return G


def build(h, cfg):
    fl = set(cfg.split(','))
    if RCM:
        import copy, reclaim as RC
        B, W0 = RC.base(h, cfg)
        rcx = RC.reclaim(B, W0, lazy=300, verbose=False)
        B = RC.assemble(B, W0, rcx, h_targets=True)
        if 'defer' in fl: RC.gfix(B, verbose=False)
        W = copy.deepcopy([list(o) for o in B['ops']]); badG = e2e7.garbage(B, W)
        print('reclaim', dict(rcx['stat']), 'R', B['R'], flush=True)
        return B, W, badG
    B = c7.build(h, allE='allE' in fl, lift='lift' in fl, defer='defer' in fl, vleaf='vleaf' in fl, alt='alt' in fl,
                 links='links' in fl, pasm='pasm' in fl, ivec='ivec' in fl, clos='clos' in fl,
                 pstar=(4 if 'dual+' in fl else 0), verbose=False,
                 dag=(os.path.join(HERE, '..', '..', 'certificates', 'round8', 'pr117_dag.json.gz') if 'dag117' in fl else None))
    W = e2e7.updates(B); badG = e2e7.garbage(B, W)
    return B, W, badG


def run(h, B, W, badG, eta, seed, control=None):
    T = B['T']; m = h * h; v = len(T); N = v * v; R = B['R']; Fs = [1 << p for p in range(h)]
    rnd = random.Random(seed)
    g = lambda: (Q(rnd.randrange(-10**6, 10**6), rnd.randrange(1, 1000)), Q(rnd.randrange(-10**6, 10**6), rnd.randrange(1, 1000)))
    fr = Frames(eta); full = canon([1 << p for p in range(m)]); Fe = fr.wt(full)
    tvec = {t: sum(1 << p for p in t) for t in T}
    tperp = {t: perp_in(Fs, [tvec[t]]) for t in T}
    D0 = {t: [kron(u, e, h) for u in tperp[t] for e in Fs] for t in T}
    lc = {}
    def lift(stage, t, lab):
        key = (stage, t, lab)
        if key not in lc:
            lc[key] = canon([kron(l, tvec[t], h) for l in lab] if stage == 1 else D0[t] + [kron(tvec[t], l, h) for l in lab])
        return lc[key]
    comp = {}
    def cpl(lab):
        if lab not in comp: comp[lab] = canon(perp_in(Fs, list(lab)))
        return comp[lab]
    fwd = W
    inv = [[roles, cpl(lab), kind, [(d, s_, -cf) for d, s_, cf in reversed(ups)] if ups else ups]
           for roles, lab, kind, ups in reversed(fwd)]
    pairs = [(a1, a2) for a1 in T for a2 in T]
    wv = {p: kron(tvec[p[0]], tvec[p[1]], h) for p in pairs}
    U = {p: canon([wv[p]]) for p in pairs}
    xin = {p: g() for p in pairs}; yin = {p: g() for p in pairs}
    val = {('X', p): xin[p] for p in pairs}; val.update({('Y', p): yin[p] for p in pairs})
    lab = {('X', p): U[p] for p in pairs}; lab.update({('Y', p): () for p in pairs})
    stage_of = {r: 1 for r in val}
    hist = Counter(); ent = dict(n=0, ok=True); prev = {}
    def move(r, L, vd, ld, data_key=None, stage=None):
        old = ld[r]
        if old == L: return
        if data_key is not None and stage == 2 and stage_of[data_key] == 1:
            a1, a2 = data_key[1]
            want = canon([kron(u, e, h) for u in tperp[a1] for e in tperp[a2]])
            if canon(perp_in(list(L), list(old))) != want: ent['ok'] = False
            ent['n'] += 1; stage_of[data_key] = 2
        ph, rk = fr.edge(old, L); hist[abs(rk)] += 1
        vd[r] = rot(vd[r], ph); ld[r] = L
        if data_key is None: prev[r] = old
    aux = [('s', q) for q in range(R)]; retset = {('s', q) for q in B['retslot']}
    aux_bad = 0; merge_bad = 0; nmerge = 0; orth_bad = 0
    G = groups_of(T, h, control)
    sizes = sorted(len(x) for x in G)
    for stage in (1, 2):
        ops = fwd if stage == 1 else inv
        first = {}; nxt = {}
        for oi, (roles, lb, kind, ups) in enumerate(ops):
            for rr in roles:
                if rr[0] == 's' and rr not in first and not (kind == 'cread' and rr in retset): first[rr] = lb
        ci = next(i for i, o in enumerate(ops) if o[2] == 'cread')
        for oi in range(ci + 1, len(ops)):
            for rr in ops[oi][0]:
                if rr in retset and rr not in nxt: nxt[rr] = ops[oi][1]
        for grp in G:
            r0 = {a: g() for a in aux}; rv = dict(r0)          # ONE bank for the whole group
            resid = {a: [] for a in aux}; lastprev = {}
            for fixed in grp:
                rl = {a: lift(stage, fixed, first[a]) for a in aux}     # relabel the incoming scratch
                start = dict(rl); prev.clear()
                if stage == 1: xmap = {t: ('X', (t, fixed)) for t in T}; ymap = {t: ('Y', (t, fixed)) for t in T}
                else: xmap = {t: ('Y', (fixed, t)) for t in T}; ymap = {t: ('X', (fixed, t)) for t in T}
                km = lambda rr: (xmap if rr[0] == 'x' else ymap)[rr[1]] if rr[0] in ('x', 'y') else rr
                lasttouch = {}
                for oi, (roles, lb, kind, ups) in enumerate(ops):
                    L = lift(stage, fixed, lb)
                    if kind == 'cread':
                        copies = {}
                        for rr in sorted(retset):
                            if stage == 2: move(rr, lift(stage, fixed, nxt[rr]), rv, rl)
                            if rl[rr] == L: copies[rr] = rv[rr]; continue
                            ph, rk = fr.edge(rl[rr], L); hist[abs(rk)] += 1; copies[rr] = rot(rv[rr], ph)
                        for rr in roles:
                            if rr[0] == 'y': move(km(rr), L, val, lab, data_key=km(rr), stage=stage)
                        continue_ups = ups
                        for d, s_, cf in continue_ups:
                            dk = km(d); val[dk] = gadd(val[dk], copies[s_], cf)
                        continue
                    for rr in roles:
                        if rr[0] in ('x', 'y'): move(km(rr), L, val, lab, data_key=km(rr), stage=stage)
                        else:
                            before = rl[rr]; move(rr, L, rv, rl)
                            lasttouch[rr] = (before != L)          # was this op a frame step?
                    for d, s_, cf in ups or ():
                        dk, sk = km(d), km(s_)
                        dv = val if dk[0] in ('X', 'Y') else rv; sv = val if sk[0] in ('X', 'Y') else rv
                        dv[dk] = gadd(dv[dk], sv[sk], cf)
                for a in aux:
                    Rb = canon(perp_in(list(rl[a]), list(start[a])))
                    assert len(Rb) == len(rl[a]) - len(start[a])
                    resid[a].append(Rb)
                    lastprev[a] = (prev.get(a), rl[a], lasttouch.get(a, False))
            for a in aux:
                Gs = [x for Rb in resid[a] for x in Rb]
                if control == 'wrapfix':
                    Gs = [x for Rb in resid[a][:1] for x in Rb]
                Gc = canon(Gs)
                if len(Gc) != sum(len(Rb) for Rb in resid[a]) or any(dot(x, y) for i, Rb in enumerate(resid[a]) for Rc in resid[a][i + 1:] for x in Rb for y in Rc):
                    orth_bad += 1
                if control != 'nofix':
                    Rfix = canon(perp_in([1 << p for p in range(m)], list(Gc)))
                    if Rfix:
                        ref = (Fe - sum(0 for _ in ()) - fr.wt(Gc)) % 4
                        ph = fr.child(Rfix); hist[len(Rfix)] += 1
                        if ph is None: ph = (Fe - fr.wt(Gc)) % 4
                        assert ph == ref, 'fix-up child != reference'
                        rv[a] = rot(rv[a], ph)
                        # merge check: the fix-up merged with the last core's last step (when it is the last op)
                        pv, lastF, stepped = lastprev[a]
                        if stepped and pv is not None and len(lastF) > len(pv):
                            Rl = canon(perp_in(list(lastF), list(pv)))
                            Rm = canon(list(Rl) + list(Rfix))
                            if len(Rm) == len(Rl) + len(Rfix):
                                nmerge += 1; phm = fr.child(Rm)
                                phl = fr.child(Rl)
                                if phm is None: phm = (fr.wt(canon(list(Rm))) ) % 4
                                if phl is None: phl = fr.wt(Rl) % 4
                                if phm is not None and phl is not None and (phm - phl - ph) % 4: merge_bad += 1
                if rv[a] != rot(r0[a], Fe): aux_bad += 1
    for p in pairs:
        move(('X', p), full, val, lab)
        move(('Y', p), canon(perp_in([1 << q for q in range(m)], list(U[p]))), val, lab)
    hist[1] += N
    bad = dict(A=0, B=0, Xout=0, Yout=0)
    for p in pairs:
        w = wv[p]; Pe = (9 * dot(w, eta)) % 4; Ee = (Fe - Pe) % 4
        A = val[('X', p)]; Bv = val[('Y', p)]
        if A != rot(yin[p], Fe + 2): bad['A'] += 1
        if Bv != gadd(rot(xin[p], Fe - 2 * Pe), rot(yin[p], Ee), 1): bad['B'] += 1
        Bc = gadd(Bv, rot(A, -Pe), 1)
        if rot(Bc, 9) != rot(xin[p], popc(eta ^ w) % 4): bad['Yout'] += 1
        if rot(A, 2) != rot(yin[p], Fe): bad['Xout'] += 1
    wk = c7.walk(B)
    pred = share.shared(B, wk, groups=[len(x) for x in G], mode='own')
    hm = {k: n for k, n in hist.items() if n} == pred['H']
    ok = aux_bad == 0 and not any(bad.values()) and ent['ok'] and not badG and orth_bad == 0 and merge_bad == 0
    return dict(h=h, control=control, PASS=ok, aux_not_restored=aux_bad, **bad, orth_bad=orth_bad, merge_checked=nmerge,
                merge_bad=merge_bad, garbage_checks=dict(badG), entrances_ok=ent['ok'], hist_matches_shared=hm,
                groups=len(G), sizes=dict(Counter(sizes)), R=R, deferred=len(B.get('sigma', {})))


if __name__ == '__main__':
    h = int(sys.argv[1]); cfg = sys.argv[2]; npts = int(sys.argv[3])
    control = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] not in ('reclaim', '-') else None
    B, W, badG = build(h, cfg)
    rr = random.Random(12345)
    for i in range(npts):
        eta = rr.getrandbits(h * h)
        try: print(run(h, B, W, badG, eta, 1000 + i, control), flush=True)
        except AssertionError as e: print(dict(h=h, control=control, FAILED_ASSERT=str(e)), flush=True)
