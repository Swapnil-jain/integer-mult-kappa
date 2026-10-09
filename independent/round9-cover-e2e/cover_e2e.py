"""Exact Q(i) replay of PR #130's three-stage Cayley cover, built on our own local complex word (c7.build) and
e2e7's dual-point model (C_U = H diag(i^wt(p_U x)) H; at a fixed dual point eta in F2^m the network is a scalar
network over Q(i); a role carries a frame; before a gate every participant moves to the gate's frame; every move is
one whole-residual child whose phase is checked against the reference ratio).

Geometry: E = A ⊥ B ⊥ C, dim A = h, dim B = dim C = h-1, m = 3h-2. Vertices are orthogonal maps k; a local label U
(subset of A) at invocation (stage s, vertex k) is the global label k(U + O_s), O = 0, B, B+C. R12_S / R23_S are the
orthogonal involutions of the note (tau_{q+e_d} basis of A cap q-perp exchanged with B or C). Stage 2 runs the
inverse word at complemented local labels with register roles exchanged (source Y, target X).

Tracked positions (g = I, S) for S in TS: stage 1 at R12_S, stage 2 at I, stage 3 at R23_S. The other ports of those
invocations get fresh data at their stage's entrance labels. Every auxiliary starts with arbitrary Q(i) scratch at
its first frame and gets its tail child (residual (last)^perp + start) at the end; it must end at F times its start.
Checks for tracked positions: zero data connectors (the entrance label of the next stage equals the exit label of the
previous one), X_out = -F y, Y_out = F T^-2 x (T = C_q), aux restored, all child phases equal reference ratios.
Controls (each must FAIL): fwdmid (stage 2 runs the forward word), badoffset (stage 3 offset B instead of B+C),
skipread (one deferred readout dropped), nonorth (R12 replaced by a non-orthogonal map with the same image of A).
Usage: python3 cover_e2e.py H flags NPTS NTRACK [control]"""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from fractions import Fraction as Q
from collections import Counter
import c7, e2e7
from e2e7 import Frames, rot, gadd, canon, updates, garbage
from frames import perp_in, dot

popc = lambda x: bin(x).count('1')


def geometry(h, control=None):
    m = 3 * h - 2
    A = [1 << i for i in range(h)]; Bb = [1 << (h + i) for i in range(h - 1)]; Cb = [1 << (2 * h - 1 + i) for i in range(h - 1)]
    def apply(M, x):
        y = 0; i = 0
        while x:
            if x & 1: y ^= M[i]
            x >>= 1; i += 1
        return y
    def solve(basis, imgs):
        # column images of the standard basis for the linear map basis[t] -> imgs[t] (basis spans F2^m)
        piv = {}
        for t, b in enumerate(basis):
            x, c = b, 1 << t
            for p, (y, cy) in piv.items():
                if x & p: x ^= y; c ^= cy
            assert x, 'basis dependent'
            p = x & -x
            for pp in list(piv):
                y, cy = piv[pp]
                if y & p: piv[pp] = (y ^ x, cy ^ c)
            piv[p] = (x, c)
        M = []
        for j in range(m):
            c = piv[1 << j][1]; y = 0; t = 0
            while c:
                if c & 1: y ^= imgs[t]
                c >>= 1; t += 1
            M.append(y)
        return M
    R = {}
    for S in c7_triples(h):
        q = sum(1 << i for i in S); d = next(i for i in range(h) if i not in S); w = q ^ (1 << d)
        tau = lambda x: x ^ (w if dot(x, w) else 0)
        kb = [tau(1 << i) for i in range(h) if i != d]
        basis = [q] + kb + Bb + Cb
        for nm, tgt in (('12', Bb), ('23', Cb)):
            if nm == '12': imgs = [q] + Bb + kb + Cb
            else: imgs = [q] + Cb + Bb + kb
            M = solve(basis, imgs)
            if control == 'nonorth' and nm == '12':
                # same image of A and of A cap q-perp, but compose with a non-orthogonal shear inside B
                if len(Bb) >= 2:
                    sh = [x for x in range(m)]
                    G = [1 << j for j in range(m)]; G[h] = (1 << h) | (1 << (h + 1))   # e_h -> e_h + e_{h+1}
                    M = [apply(G, y) for y in M]
            R[(nm, S)] = M
    return dict(m=m, A=A, B=Bb, C=Cb, R=R, apply=apply)


def c7_triples(h):
    from itertools import combinations
    return [tuple(c) for c in combinations(range(h), 3)]


def run(h, flags, eta, seed, ntrack, control=None):
    B = c7.build(h, allE='allE' in flags, lift='lift' in flags, defer='defer' in flags, late='late' in flags,
                 vleaf='vleaf' in flags, climb='climb' in flags, pstar=(4 if 'dual+' in flags else 0),
                 links='links' in flags, pasm='pasm' in flags, ivec='ivec' in flags, prio='prio' in flags,
                 clos='clos' in flags, dag=None, alt='alt' in flags, verbose=False)
    W = updates(B); badG = garbage(B, W)
    T = B['T']
    if control == 'skipread':           # drop a deferred readout that reaches a tracked target
        k = next(i for i, w in enumerate(W) if w[2] == 'read' and w[1] and w[3] and any(d[1] in T[:ntrack] for d, _, _ in w[3]))
        W[k][3] = []
    R = B['R']; Fs = [1 << p for p in range(h)]
    Gm = geometry(h, control); m = Gm['m']; ap = Gm['apply']
    assert sorted(T) == sorted(c7_triples(h)) or set(T) <= set(c7_triples(h))
    rnd = random.Random(seed)
    g = lambda: (Q(rnd.randrange(-10**6, 10**6), rnd.randrange(1, 1000)), Q(rnd.randrange(-10**6, 10**6), rnd.randrange(1, 1000)))
    fr = Frames(eta)
    full = canon([1 << p for p in range(m)])
    comp = {}
    def cpl(lab):
        if lab not in comp: comp[lab] = canon(perp_in(Fs, list(lab)))
        return comp[lab]
    inv = [[roles, cpl(lab), kind, [(d, s_, -cf) for d, s_, cf in reversed(ups)] if ups else ups]
           for roles, lab, kind, ups in reversed(W)]
    offs = {1: [], 2: Gm['B'], 3: (Gm['B'] if control == 'badoffset' else Gm['B'] + Gm['C'])}
    gl = {}
    def glob(M, stage, lab):
        key = (id(M), stage, lab)
        if key not in gl: gl[key] = canon([ap(M, x) for x in list(lab) + offs[stage]]) if (lab or offs[stage]) else ()
        return gl[key]
    I = [1 << j for j in range(m)]
    kap = [1 << (j + h) if j < h else (1 << (j - h) if j < 2 * h else 1 << j) for j in range(m)]
    if control == 'badkappa':                    # partner whose active space meets A (shift by h-1): not orthogonal
        kap = [1 << ((j + h - 1) % m) for j in range(m)]
    def glob_nooff(M, lab):
        return [ap(M, x) for x in lab]
    TS = T[:ntrack]
    qv = {S: sum(1 << i for i in S) for S in T}
    hist = Counter(); stats = Counter()
    def move(vd, ld, r, L):
        old = ld[r]
        if old == L: return
        ph, rk = fr.edge(old, L); hist[abs(rk)] += 1
        vd[r] = rot(vd[r], ph); ld[r] = L
    retset = {('s', q) for q in B['retslot']}
    # tracked data registers: values and global labels
    xin = {S: g() for S in TS}; yin = {S: g() for S in TS}
    val = {}; lab = {}
    for S in TS:
        val[('X', S)] = xin[S]; lab[('X', S)] = canon([qv[S]])     # entrance gauge T = C_q on x (declared, no op)
        val[('Y', S)] = yin[S]; lab[('Y', S)] = ()
    aux_bad = 0
    def invocation(stage, M, carried, shared=None):
        """run one local word at vertex M; carried = {port S: tracked key} for ports whose data arrive from the
        previous stage; other ports get fresh data at the stage's entrance labels."""
        nonlocal aux_bad
        fwd = (stage != 2) or control == 'fwdmid'
        ops = W if fwd else inv
        vd = {}; ld = {}
        for S in T:
            if S in carried:
                for rr in ('X', 'Y'):
                    k = (rr, carried[S]); vd[(rr, S)] = val[k]; ld[(rr, S)] = lab[k]
            else:
                if stage == 2 and not fwd:    # inverse word: register X plays local 'y' (starts at q), Y plays 'x' (starts at 0)
                    vd[('X', S)] = g(); ld[('X', S)] = glob(M, stage, canon([qv[S]]))
                    vd[('Y', S)] = g(); ld[('Y', S)] = glob(M, stage, ())
                else:
                    vd[('X', S)] = g(); ld[('X', S)] = glob(M, stage, canon([qv[S]]))
                    vd[('Y', S)] = g(); ld[('Y', S)] = glob(M, stage, ())
        # entrance check (no connectors): the carried labels must equal this stage's entrance labels
        for S in carried:
            ex = glob(M, stage, canon([qv[S]])), glob(M, stage, ())
            if (ld[('X', S)], ld[('Y', S)]) != ex: stats['connector_mismatch'] += 1
        if fwd: km = lambda rr: ('X', rr[1]) if rr[0] == 'x' else ('Y', rr[1]) if rr[0] == 'y' else rr
        else: km = lambda rr: ('Y', rr[1]) if rr[0] == 'x' else ('X', rr[1]) if rr[0] == 'y' else rr
        aux = [('s', q) for q in range(R)]
        first = {}; nxt = {}
        for roles, lb, kind, ups in ops:
            for rr in roles:
                if rr[0] == 's' and rr not in first and not (kind == 'cread' and rr in retset): first[rr] = lb
        ci = next(i for i, o in enumerate(ops) if o[2] == 'cread')
        for oi in range(ci + 1, len(ops)):
            for rr in ops[oi][0]:
                if rr in retset and rr not in nxt: nxt[rr] = ops[oi][1]
        if shared is None: r0 = {a: g() for a in aux}; rv = dict(r0)
        else: rv = shared['rv']                      # consecutive sharing: scratch left by the previous core, relabelled
        rl = {a: glob(M, stage, first[a]) for a in aux}; start = dict(rl)
        for roles, lb, kind, ups in ops:
            L = glob(M, stage, lb)
            if kind == 'cread':
                copies = {}
                for rr in sorted(retset):
                    if not fwd: move(rv, rl, rr, glob(M, stage, nxt[rr]))
                    if rl[rr] == L: copies[rr] = rv[rr]; continue
                    ph, rk = fr.edge(rl[rr], L); hist[abs(rk)] += 1; copies[rr] = rot(rv[rr], ph)
                for rr in roles:
                    if rr[0] == 'y': move(vd, ld, km(rr), L)
                for d, s_, cf in ups:
                    dk = km(d); vd[dk] = gadd(vd[dk], copies[s_], cf)
                continue
            for rr in roles:
                if rr[0] in ('x', 'y'): move(vd, ld, km(rr), L)
                else: move(rv, rl, rr, L)
            for d, s_, cf in ups or ():
                dk, sk = km(d), km(s_)
                dv = vd if dk[0] in ('X', 'Y') else rv; sv = vd if sk[0] in ('X', 'Y') else rv
                dv[dk] = gadd(dv[dk], sv[sk], cf)
        Fe = fr.wt(full)
        if shared is not None:
            # record this core's net action (start -> last, offsets cancel) and its local start/last without offsets
            for a in aux:
                shared['net'][a] = (shared['net'].get(a, 0) + fr.wt(rl[a]) - fr.wt(start[a])) % 4
                Moff = glob_nooff(M, offs[stage])
                shared['last'].setdefault(a, []).extend(perp_in(list(rl[a]), Moff) if Moff else list(rl[a]))
                shared['start'].setdefault(a, []).extend(perp_in(list(start[a]), Moff) if Moff else list(start[a]))
            aux = []
        for a in aux:                               # tail child: residual (last)^perp + start
            last = rl[a]; st = start[a]
            Rw = canon(perp_in([1 << p for p in range(m)], list(last)) + list(st))
            assert len(Rw) == m - len(last) + len(st)
            ph = fr.child(Rw); hist[len(Rw)] += 1
            if ph is None: ph = (Fe - fr.wt(last) + fr.wt(st)) % 4
            assert ph == (Fe - fr.wt(last) + fr.wt(st)) % 4, 'tail child != reference'
            rv[a] = rot(rv[a], ph)
            if rv[a] != rot(r0[a], Fe): aux_bad += 1
        # exit: every data register must have reached its local end label (X: A, Y: A - q)
        for S in T:
            want = (glob(M, stage, canon(Fs)), glob(M, stage, cpl(canon([qv[S]]))))
            if (ld[('X', S)], ld[('Y', S)]) != want: stats['exit_label_mismatch'] += 1
        for S, key in carried.items():
            for rr in ('X', 'Y'): val[(rr, key)] = vd[(rr, S)]; lab[(rr, key)] = ld[(rr, S)]
        stats['invocations'] += 1
    def kmul(M):                                 # the vertex M*kappa as column images
        return [ap(M, kap[j]) for j in range(m)]
    def paired(stage, M, carried):
        nonlocal aux_bad
        r0 = {('s', q): g() for q in range(R)}
        sh = dict(rv=dict(r0), net={}, last={}, start={})
        invocation(stage, M, carried, sh); invocation(stage, kmul(M), {}, sh)
        Fe = fr.wt(full)
        for a in r0:
            last = canon(sh['last'][a]); st = canon(sh['start'][a])
            Rw = canon(perp_in([1 << p for p in range(m)], list(last)) + list(st))
            assert len(Rw) == m - len(last) + len(st), 'pair tail rank'
            stats['pair_tail_' + str(len(Rw))] += 0
            ph = fr.child(Rw); hist[len(Rw)] += 1
            ref = (Fe - sh['net'][a]) % 4
            if ph is None: ph = ref
            assert ph == ref, 'pair tail child != reference'
            sh['rv'][a] = rot(sh['rv'][a], ph)
            if sh['rv'][a] != rot(r0[a], Fe): aux_bad += 1
    run_inv = paired if control in ('pair', 'badkappa') or 'pair' in flags else (lambda st, M, c: invocation(st, M, c))
    for S in TS: run_inv(1, Gm['R'][('12', S)], {S: S})
    run_inv(2, I, {S: S for S in TS})
    for S in TS: run_inv(3, Gm['R'][('23', S)], {S: S})
    Fe = fr.wt(full); bad = Counter()
    for S in TS:
        if lab[('X', S)] != full: bad['X_label'] += 1
        if lab[('Y', S)] != canon(perp_in([1 << p for p in range(m)], [qv[S]])): bad['Y_label'] += 1
        Pe = fr.wt(canon([qv[S]]))
        if val[('X', S)] != rot(yin[S], Fe + 2): bad['X_out'] += 1
        if val[('Y', S)] != rot(xin[S], Fe - 2 * Pe): bad['Y_out'] += 1
        # alternative endpoint convention (-F T^-1 y, F x): must fail whenever T is not a scalar here (Pe != 0)
        if Pe % 4:
            stats['Pe_nonzero'] += 1
            if val[('X', S)] == rot(yin[S], Fe + 2 - Pe) and val[('Y', S)] == rot(xin[S], Fe): stats['alt_convention_matches'] += 1
    ok = aux_bad == 0 and not bad and not stats['connector_mismatch'] and not stats['exit_label_mismatch'] and not badG
    return dict(h=h, m=m, control=control, PASS=ok, aux_not_restored=aux_bad, bad=dict(bad), connectors=stats['connector_mismatch'],
                exit_mismatch=stats['exit_label_mismatch'], invocations=stats['invocations'], tracked=len(TS), R=R,
                garbage=dict(badG), alt_children=fr.nalt, Pe_nonzero=stats['Pe_nonzero'], alt_convention_matches=stats['alt_convention_matches'])


if __name__ == '__main__':
    h = int(sys.argv[1]); flags = set(sys.argv[2].split(',')); npts = int(sys.argv[3]); ntrack = int(sys.argv[4])
    control = sys.argv[5] if len(sys.argv) > 5 else None
    rr = random.Random(4242)
    for i in range(npts):
        eta = rr.getrandbits(3 * h - 2)
        try: print(run(h, flags, eta, 77 + i, ntrack, control), flush=True)
        except AssertionError as e: print(dict(h=h, control=control, FAILED_ASSERT=str(e)), flush=True)
