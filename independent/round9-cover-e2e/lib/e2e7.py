"""Exact end-to-end replay of the round-7 complex word (c7.build) over Q(i), both stages, arbitrary scratch, every
frame change compiled as ONE whole-residual child (PR #10 signed identity on an explicit orthonormal residual basis,
asserted equal to the reference phase ratio), copied retained centres, the wrap child of every auxiliary
(F_last -> start + full, residual F_last^perp + start), and negative controls.

The model is gate15c/e2e_copy.py's (C_U = H diag(i^wt(p_U x)) H; at a fixed dual point eta the network is a scalar
network over Q(i); a role carries a frame; before a gate every participant moves to the gate's frame). Our word adds:
direct garbage readouts y_T -= G_{T,u} z_u (G = the garbage map of pass two, computed symbolically here and checked
to be supported on the readout gate's roles, each readout reading slot u before anything writes it), deferred
readouts at sigma_u (the slot STARTS at sigma_u: no child from 0, its wrap child carries sigma_u), lifted frames,
late copies, V leaves. Stage two is the inverse word (ops reversed, updates inverted, labels complemented).
Target: (x, y) -> (F y, F x) on every data pair, every auxiliary restored (to F r0), histogram == c7.walk's.
Controls: 'nowrapsig' (wrap child without the start part), 'skipread' (drop one deferred readout),
'lateread' (a deferred readout after its slot's first write), 'badbasis' (non-orthonormal child basis).
Usage: python3 e2e7.py H flags(comma) NPTS [control]"""
import sys, random
from fractions import Fraction as Q
from collections import Counter, defaultdict
import c7
from frames import echelon, red, perp_in, ok_res, dot

HALF = Q(1, 2)
popc = lambda x: bin(x).count('1')
def gadd(a, b, cf): return (a[0] + cf * b[0], a[1] + cf * b[1])
def rot(a, k):
    k %= 4
    return a if k == 0 else (-a[1], a[0]) if k == 1 else (-a[0], -a[1]) if k == 2 else (a[1], -a[0])
def kron(a, b, h):
    x = 0; i = 0
    while a:
        if a & 1: x |= b << (i * h)
        a >>= 1; i += 1
    return x
def canon(B): return tuple(sorted(echelon(list(B)).values()))

def orthonormal_basis(V):
    out = []; W = list(V)
    while W:
        e = next((x for x in W if popc(x) & 1), None)
        if e is None:
            f = out.pop(); u = W[0]; w = next(x for x in W if dot(u, x))
            out += [f ^ u, f ^ w, f ^ u ^ w]; W = perp_in(W, [u, w]); continue
        out.append(e); W = perp_in(W, [e])
    assert all(dot(a, b) == (i == j) for i, a in enumerate(out) for j, b in enumerate(out))
    return out

class Frames:
    def __init__(s, eta, badbasis=False): s.eta = eta; s.bad = badbasis; s.c = {}; s.w = {}; s.rc = {}; s.nalt = 0
    def wt(s, U):
        if U in s.w: return s.w[U]
        B = list(U); n = len(B)
        if n == 0: s.w[U] = 0; return 0
        rows = [(sum(dot(B[i], B[j]) << j for j in range(n)), dot(B[i], s.eta)) for i in range(n)]
        sol = []
        for col in range(n):
            bit = 1 << col
            k = next(q for q, r in enumerate(rows) if r[0] & bit)
            pr = rows.pop(k)
            rows = [(r[0] ^ pr[0], r[1] ^ pr[1]) if r[0] & bit else r for r in rows]
            sol = [(r[0] ^ pr[0], r[1] ^ pr[1]) if r[0] & bit else r for r in sol]; sol.append(pr)
        p = 0
        for r in sol:
            assert popc(r[0]) == 1
            if r[1]: p ^= B[r[0].bit_length() - 1]
        s.w[U] = popc(p) % 4; return s.w[U]
    def child(s, R):
        """phase of the whole-residual child on residual R (orthonormal basis, PR #10 signed identity)."""
        if R in s.rc: return s.rc[R]
        assert ok_res(list(R)), 'residual degenerate or alternating'
        if not any(popc(x) & 1 for x in R):
            s.rc[R] = None; return None              # alternating: PR #24 Gauss normal form, phase = reference
        basis = list(R) if s.bad else orthonormal_basis(list(R))
        y = [dot(b, s.eta) for b in basis]; Bs = [v for v, b in enumerate(basis) if popc(b) % 4 == 3]
        ph = (sum(yv ^ (v in Bs) for v, yv in enumerate(y)) + 3 * len(Bs)) % 4
        s.rc[R] = ph; return ph
    def edge(s, U, V):
        if (U, V) in s.c: return s.c[U, V]
        up = len(V) >= len(U); lo, hi = (U, V) if up else (V, U)
        bh = echelon(list(hi)); assert all(red(bh, x) == 0 for x in lo), 'frames not nested'
        R = canon(perp_in(list(hi), list(lo))); assert len(R) == len(hi) - len(lo)
        ph = s.child(R)
        if ph is None: ph = (s.wt(hi) - s.wt(lo)) % 4; s.nalt += 1
        if not up: ph = (-ph) % 4
        if not s.bad: assert ph == (s.wt(V) - s.wt(U)) % 4, 'child != reference ratio'
        s.c[U, V] = (ph, len(V) - len(U)); return s.c[U, V]


def updates(B):
    """stage-one word as a list of (roles, label, kind, ups); ups = [(dst, src, coef)] applied in order.
    'cread' ups read copies of the retained slots. 'high' undoes pass two's gate / late-copy updates."""
    c = B['c']; kk = B['kk']; h = B['h']; T = B['T']; args = c.args
    gate_of = {n: (ins, outs) for n, ins, outs in kk['gates']}
    rout = kk['rout']; allE = c.allE
    def scat(S):
        if allE: return [(('s', rout[('E', i)]), Q(1, h - 3) - (HALF if i in S else 0)) for i in range(h)]
        last = h - 1
        if last not in S: terms = [(('*',), Q(1))] + [(('E', i), -HALF) for i in S]
        else: terms = [(('*',), (5 - h) * HALF)] + [(('E', i), HALF) for i in range(h - 1) if i not in S]
        return [(('s', rout[nm]), cf) for nm, cf in terms]
    out = []; p2 = []; after_high = False
    for op in B['ops']:
        roles, lab, kind = op[:3]
        if kind == 'gate':
            n, excl = op[3]; ins, outs = gate_of[n]
            ups = ([(('s', ins[0]), ('s', ins[1]), 1)] if args[n] else []) + \
                  [(('s', o), ('s', outs[0]), 1) for o in outs[1:] if o not in excl]
            p2 += ups
        elif kind == 'latecopy':
            ups = [(roles[1], roles[0], 1)]; p2 += ups
        elif kind in ('xcopy', 'vgate'):
            ups = [(roles[1], roles[0], -1 if after_high else 1)]
        elif kind == 'inject':
            if len(op) > 3 and op[3] is not None:
                i = op[3]; ups = [(roles[0], roles[1], c.pieces[i][2] * HALF)]
            else: ups = []
        elif kind == 'cread':
            ups = [(('y', S), sl, cf) for S in T for sl, cf in scat(S)]
        elif kind == 'high':
            ups = [(d, s_, -cf) for d, s_, cf in reversed(p2)]; after_high = True
        elif kind == 'read':
            ups = None                                   # filled from G
        out.append([roles, lab, kind, ups])
    return out


def garbage(B, W):
    """G[T][q]: coefficient of the initial garbage z_q in y_T after pass two (x = 0), symbolically; then fill the
    readouts (y_T -= G z) and check support and purity."""
    val = defaultdict(dict)
    for q in range(B['R']): val[('s', q)] = {q: Q(1)}
    for roles, lab, kind, ups in W:
        if kind in ('read',): continue
        if kind == 'high': break
        for d, s_, cf in ups:
            if s_[0] == 'x': continue                    # x = 0
            dv = val[d]
            for k, v in val[s_].items():
                dv[k] = dv.get(k, 0) + cf * v
                if dv[k] == 0: del dv[k]
    G = defaultdict(dict)                                # q -> {T: coef}
    for S in B['T']:
        for q, v in val[('y', S)].items(): G[q][S] = v
    touched = set(); bad = Counter()
    for w in W:
        roles, lab, kind, ups = w
        if kind == 'read':
            q = roles[0][1]
            if ('s', q) in touched: bad['impure'] += 1
            tg = set(r[1] for r in roles[1:])
            if not set(G[q]) <= tg: bad['support'] += 1
            w[3] = [(('y', S), ('s', q), -cf) for S, cf in G[q].items()]
        elif ups:
            for d, s_, cf in ups: touched.add(d)
    return bad


def run(h, flags, eta, seed, control=None):
    B = c7.build(h, allE='allE' in flags, lift='lift' in flags, defer='defer' in flags, late='late' in flags,
                 vleaf='vleaf' in flags, climb='climb' in flags, pstar=(6 if 'ival+' in flags else 5 if 'ival' in flags else 4 if 'dual+' in flags else 3 if 'dual' in flags else 2 if 'pstar2' in flags else 1 if 'pstar' in flags else 0), links='links' in flags, pasm='pasm' in flags, ivec='ivec' in flags, prio='prio' in flags, clos='clos' in flags, dag=None, alt='alt' in flags, verbose=False)
    W = updates(B); badG = garbage(B, W)
    if control == 'skipread':
        k = next(i for i, w in enumerate(W) if w[2] == 'read' and w[1] and w[3]); W[k][3] = []
    if control == 'lateread':                           # a deferred readout moved after its slot's first write
        k = next(i for i, w in enumerate(W) if w[2] == 'read' and w[1] and w[3]); q = W[k][0][0]
        j = next(i for i, w in enumerate(W) if i > k and w[2] != 'read' and w[3] and any(d == q for d, _, _ in w[3]))
        W.insert(j + 1, W.pop(k))
    c = B['c']; T = B['T']; m = h * h; v = len(T); N = v * v; R = B['R']; Fs = [1 << p for p in range(h)]
    rnd = random.Random(seed)
    g = lambda: (Q(rnd.randrange(-10**6, 10**6), rnd.randrange(1, 1000)), Q(rnd.randrange(-10**6, 10**6), rnd.randrange(1, 1000)))
    fr = Frames(eta, control == 'badbasis'); full = canon([1 << p for p in range(m)])
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
    hist = Counter(); ent = dict(n=0, ok=True)
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
    aux = [('s', q) for q in range(R)]; retset = {('s', q) for q in B['retslot']}
    aux_bad = 0
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
        for fixed in T:
            r0 = {a: g() for a in aux}; rv = dict(r0)
            rl = {a: lift(stage, fixed, first[a]) for a in aux}     # every auxiliary starts at its first frame
            start = dict(rl)
            if stage == 1: xmap = {t: ('X', (t, fixed)) for t in T}; ymap = {t: ('Y', (t, fixed)) for t in T}
            else: xmap = {t: ('Y', (fixed, t)) for t in T}; ymap = {t: ('X', (fixed, t)) for t in T}
            km = lambda rr: (xmap if rr[0] == 'x' else ymap)[rr[1]] if rr[0] in ('x', 'y') else rr
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
                    for d, s_, cf in ups:
                        dk = km(d); val[dk] = gadd(val[dk], copies[s_], cf)
                    continue
                for rr in roles:
                    if rr[0] in ('x', 'y'): move(km(rr), L, val, lab, data_key=km(rr), stage=stage)
                    else: move(rr, L, rv, rl)
                for d, s_, cf in ups or ():
                    dk, sk = km(d), km(s_)
                    dv = val if dk[0] in ('X', 'Y') else rv; sv = val if sk[0] in ('X', 'Y') else rv
                    dv[dk] = gadd(dv[dk], sv[sk], cf)
            Fe = fr.wt(full)
            for a in aux:                                   # wrap child: residual (last)^perp + start
                last = rl[a]; st = start[a] if control != 'nowrapsig' else ()
                Rw = canon(perp_in([1 << p for p in range(m)], list(last)) + list(st))
                assert len(Rw) == m - len(last) + len(st)
                ph = fr.child(Rw); hist[len(Rw)] += 1
                if ph is None: ph = (Fe - fr.wt(last) + fr.wt(st)) % 4
                if control != 'badbasis': assert ph == (Fe - fr.wt(last) + fr.wt(st)) % 4, 'wrap child != reference'
                rv[a] = rot(rv[a], ph)
                if rv[a] != rot(r0[a], Fe): aux_bad += 1
    for p in pairs:
        move(('X', p), full, val, lab)
        move(('Y', p), canon(perp_in([1 << q for q in range(m)], list(U[p]))), val, lab)
    hist[1] += N
    Fe = fr.wt(full); bad = dict(A=0, B=0, Xout=0, Yout=0)
    for p in pairs:
        w = wv[p]; Pe = (9 * dot(w, eta)) % 4; Ee = (Fe - Pe) % 4
        A = val[('X', p)]; Bv = val[('Y', p)]
        if A != rot(yin[p], Fe + 2): bad['A'] += 1
        if Bv != gadd(rot(xin[p], Fe - 2 * Pe), rot(yin[p], Ee), 1): bad['B'] += 1
        Bc = gadd(Bv, rot(A, -Pe), 1)
        if rot(Bc, 9) != rot(xin[p], popc(eta ^ w) % 4): bad['Yout'] += 1
        if rot(A, 2) != rot(yin[p], Fe): bad['Xout'] += 1
    wk = c7.walk(B)
    ok = aux_bad == 0 and not any(bad.values()) and ent['ok'] and not badG
    return dict(h=h, flags=flags, control=control, PASS=ok, aux_not_restored=aux_bad, **bad, garbage_checks=dict(badG),
                entrances_ok=ent['ok'], hist_matches_walk=({k: n for k, n in hist.items() if n} == wk['H']),
                R=R, deferred=len(B['sigma']), alternating_children=fr.nalt)


if __name__ == '__main__':
    h = int(sys.argv[1]); flags = set(sys.argv[2].split(',')); npts = int(sys.argv[3])
    control = sys.argv[4] if len(sys.argv) > 4 else None
    rr = random.Random(12345)
    for i in range(npts):
        eta = rr.getrandbits(h * h)
        try: print(run(h, flags, eta, 1000 + i, control), flush=True)
        except AssertionError as e: print(dict(h=h, flags=sorted(flags), control=control, FAILED_ASSERT=str(e)), flush=True)
