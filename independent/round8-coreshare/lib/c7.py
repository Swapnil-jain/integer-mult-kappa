"""Round-7 levers on the two-stage complex network (NStar3 producer, retained totals, copied centres, full batching).

The stage-one word is built explicitly as a list of ops (roles, frame), every frame a subspace of F2^h (stage-one
lift L (x) t_a2). Stage two is its complement time reversal (frames L -> L^perp, op order reversed), so its role
chains have the same rank multiset; the histogram counts stage one twice for auxiliaries and copies, and pairs the
stage-one data chains with their mirrors for the data wires (exactly as hist.py does).

Options
  allE   retained basis E_0..E_{h-1} instead of (E_0..E_{h-2}, T)   (sum_j E_j = (h-3) T)
  lift   lifted frames: M_n = Z_n = (cap of the frames of n's users) cap (cap of t_S^perp over n's pieces), when every
         step it creates is a valid residual; else the label (top-down, users first)
  defer  B-defer: garbage readouts are direct (y_T -= C_{T,u} z_u). Slots untouched by phase 1 (the gates of the
         retained totals' cones) read out after the centre reads, at sigma_u <= F0(u) (first pass-two frame); the
         slot's chain starts at sigma_u and the wrap child F -> sigma_u + full has rank m - h + dim sigma_u.
         On every y_T the readout frames (0, deferred sigma's sorted by dim, t_T^perp) must be a chain (hosting).
  late   late copies: a fan-out copy of a cone node for a non-cone user is made in phase two (after its readout),
         from the pivot, at the pivot's frame M_n (so the copy slot is free).
Usage: python3 c7.py H [allE] [lift] [defer] [late]"""
import sys, math, pickle
from collections import Counter, defaultdict
from producer import NStar3, compile_roles
from frames import echelon, red, perp_in, ok_res, nondeg, Checker


NTRY = 6
CACHE = 400000
DEPTH = 3


class Sub:
    """subspaces of F2^h as canonical reduced-echelon tuples, with cached meet / containment / step validity."""
    def __init__(s, h):
        s.h = h; s.F = tuple(sorted(echelon([1 << p for p in range(h)]).values())); s.Fl = list(s.F)
        s.cp = {}; s.cm = {}; s.cs = {}; s.cn = {}
        import random; s.rnd = random.Random(11)
    def canon(s, vs): return tuple(sorted(echelon(list(vs)).values()))
    def perp(s, U):
        if U not in s.cp: s.cp[U] = s.canon(perp_in(s.Fl, list(U)))
        return s.cp[U]
    def meet(s, U, V):
        k = (U, V) if U <= V else (V, U)
        if len(s.cm) > CACHE: s.cm.clear()
        if k not in s.cm: s.cm[k] = s.perp(s.canon(list(s.perp(U)) + list(s.perp(V))))
        return s.cm[k]
    def le(s, U, V):
        """U <= V"""
        if len(U) > len(V): return False
        b = dict((x & -x, x) for x in V)            # V is reduced echelon keyed by lowest bit
        return all(red(b, u) == 0 for u in U)
    def step_ok(s, U, V):
        """U -> V is one valid whole-residual child (nested either way, residual nondegenerate nonalternating)."""
        if U == V: return True
        k = (U, V)
        if len(s.cs) > CACHE: s.cs.clear()
        if k not in s.cs:
            lo, hi = (U, V) if len(U) <= len(V) else (V, U)
            s.cs[k] = s.le(lo, hi) and ok_res(perp_in(list(hi), list(lo)))
        return s.cs[k]
    def fit(s, lo, cap, ok, tries=None, depth=None, rnd=None):
        """a large subspace Y with lo <= Y <= cap and ok(Y): cap itself, then nondegenerate complements of cap's
        radical through lo, then random hyperplane descents through lo (depth levels, a few per level)."""
        tries = tries or NTRY; depth = depth if depth is not None else DEPTH
        rnd = rnd or s.rnd
        if ok(cap): return cap
        lo_l = list(lo); seen = {cap}
        frontier = [cap]
        rad = s.canon(perp_in(list(cap), list(cap)))
        if rad and len(s.canon(lo_l + list(rad))) == len(lo) + len(rad):
            bb = echelon(lo_l + list(rad)); ext = []
            for x in cap:
                if red(bb, x): ext.append(x); bb = echelon(list(bb.values()) + [x])
            rv = list(rad)
            for k in range(tries):
                e2 = []
                for x in ext:
                    y = x
                    if k:
                        for r_ in rv:
                            if rnd.random() < 0.5: y ^= r_
                    e2.append(y)
                W = s.canon(lo_l + e2)
                if W in seen or len(W) != len(cap) - len(rad): continue
                seen.add(W)
                if ok(W): return W
                frontier.append(W)
        perp_lo = list(s.perp(lo))
        for d in range(depth):
            new = []
            for Wsp in frontier:
                if len(Wsp) <= len(lo) + 1: continue
                for k in range(tries):
                    f = 0
                    for b in perp_lo:
                        if rnd.random() < 0.5: f ^= b
                    if not f: continue
                    H = s.canon(perp_in(list(Wsp), [f]))
                    if len(H) != len(Wsp) - 1 or H in seen: continue
                    seen.add(H)
                    if ok(H): return H
                    new.append(H)
            rnd.shuffle(new); frontier = new[:tries]
            if not frontier: break
        return lo if lo != cap and ok(lo) else None
    def nd(s, U):
        if U not in s.cn: s.cn[U] = (not U) or nondeg(list(U))
        return s.cn[U]


def compile_roles2(c, cone, climb=False, defer=False):
    """compile_roles with one change: a cone node's pivot goes to a non-cone user when it has one (so after phase 1
    the pivot still holds the node's value at its frame, and copies for non-cone users can be made late)."""
    users = {x: [] for x in c.active}
    for n in sorted(c.active):
        if c.args[n]:
            for pos, x in enumerate(c.args[n]): users[x].append(('g', n, pos))
    for i, (_, n, _) in enumerate(c.pieces): users[n].append(('p', i))
    for name, n in c.retained: users[n].append(('r', name))
    if climb:
        # users in time order (phase-1 gates, retained read, phase-2 gates, pieces); the pivot goes to the last
        def tkey(u):
            if u[0] == 'g': return ((1 if (defer and u[1] in cone) else 3), u[1])
            if u[0] == 'r': return (2, 0)
            return (4, u[1])
        for x in c.active:
            if not c.args[x]: continue
            us = sorted(users[x], key=tkey)
            users[x] = [us[-1]] + us[:-1]
    else:
        for x in cone:
            us = users[x]
            k = next((i for i, u in enumerate(us) if (u[0] == 'g' and u[1] not in cone) or u[0] == 'p'), None)
            if k: us.insert(0, us.pop(k))
    edge = {}; src = {}; pout = {}; rout = {}; gates = []; size = 0
    for n in sorted(c.active):
        if c.args[n]: ins = (edge[n, 0], edge[n, 1]); piv = ins[0]
        else: piv = size; size += 1; ins = (piv,); src[c.triples[n - 1]] = piv
        outs = (piv,) + tuple(range(size, size + len(users[n]) - 1)); size += len(users[n]) - 1
        gates.append((n, ins, outs))
        for u, sl in zip(users[n], outs):
            if u[0] == 'g': edge[u[1], u[2]] = sl
            elif u[0] == 'p': pout[u[1]] = sl
            else: rout[u[1]] = sl
    assert size == c.roles
    return dict(size=size, gates=gates, src=src, pout=pout, rout=rout)


def build(h, allE=False, lift=False, defer=False, late=False, vleaf=False, dp=False, climb=False, pstar=False, alt=False, links=False, pasm=False, ivec=False, prio=False, dag=None, clos=False, verbose=True):
    import frames as _fr; _fr.ALT[0] = alt
    if dag:
        from dagprod import DagProducer
        c = DagProducer(dag); assert c.h == h and allE
    else:
        c = NStar3(h, allE=allE, pstar=pstar, nosplit=alt, pasm=pasm, ivec=ivec)
    ch = Checker(c); X = Sub(h)
    cone0 = set(); st = [n for _, n in c.retained]
    while st:
        x = st.pop()
        if x in cone0: continue
        cone0.add(x)
        if c.args[x]: st.extend(c.args[x])
    T = c.triples; v = len(T); m = h * h; args = c.args
    tv = lambda S: sum(1 << p for p in S)
    tperp = {S: X.perp(X.canon([tv(S)])) for S in T}
    line = {S: X.canon([tv(S)]) for S in T}
    lab = {n: X.canon(ch.label(n)) for n in c.active}
    tkey = lambda n: ((1 if (defer and n in cone0 and not clos) else 2), len(lab[n]), n)
    if links:
        from links import compile_links
        kk = compile_links(c, lab, X, tperp, tkey, verbose)
    else:
        kk = compile_roles2(c, cone0 if (late or climb) else set(), climb, defer)
    R = kk['size']
    link_tg = defaultdict(list)                          # donor -> frames its linked slot must reach
    for u, (x, pos) in kk.get('links', {}).items():
        link_tg[x].append(u)
    gates = kk['gates']; gate_of = {n: (ins, outs) for n, ins, outs in gates}
    # users (gates) and pieces of every node
    gusers = defaultdict(list); npieces = defaultdict(list); retnode = {}
    for n in c.active:
        if args[n]:
            for a in args[n]: gusers[a].append(n)
    for i, (S, n, cf) in enumerate(c.pieces): npieces[n].append(S)
    for nm, n in c.retained: retnode[n] = nm
    # ---------------- frames of the nodes ----------------
    M = dict(lab)
    nlift = Counter()
    if lift:
        for n in sorted(c.active, key=tkey, reverse=True):   # users and link targets come later in time
            if n in retnode or not args[n]: continue      # retained totals and leaves keep their labels
            Z = X.F
            for g in gusers[n]: Z = X.meet(Z, M[g])
            for u in link_tg.get(n, ()):
                Z = X.meet(Z, M[u[1]] if u[0] == 'g' else tperp[c.pieces[u[1]][0]])
            for S in npieces[n]: Z = X.meet(Z, tperp[S])
            assert X.le(lab[n], Z)
            def okY(Y):
                if not (X.nd(Y) and ok_res(list(Y)) and X.step_ok(Y, X.F)): return False
                if not all(X.step_ok(Y, M[g]) for g in gusers[n]): return False
                if not all(X.step_ok(Y, M[u[1]] if u[0] == 'g' else tperp[c.pieces[u[1]][0]]) for u in link_tg.get(n, ())): return False
                if not all(X.step_ok(Y, tperp[S]) for S in npieces[n]): return False
                return all(X.step_ok(M[a], Y) for a in args[n])      # args still at their labels here
            ok = X.fit(lab[n], Z, okY) if Z != lab[n] else None
            if ok is not None and ok != lab[n]: M[n] = ok; nlift['lifted'] += 1
            else: nlift['label'] += 1
        # args were checked against labels; re-check every DAG step with final frames, demote on failure
        bad = 1
        while bad:
            bad = 0
            for n in sorted(c.active):
                if not args[n] or M[n] == lab[n]: continue
                if not all(X.step_ok(M[a], M[n]) for a in args[n]):
                    M[n] = lab[n]; bad += 1; nlift['demoted'] += 1
            # demoting can break a user step: re-check users' steps from demoted frames
            for n in sorted(c.active):
                if not args[n]: continue
                for a in args[n]:
                    if not X.step_ok(M[a], M[n]) and M[n] != lab[n]: M[n] = lab[n]; bad += 1; nlift['demoted'] += 1
        if verbose: print('lift', dict(nlift), flush=True)
    # ---------------- slot structure ----------------
    touch = defaultdict(list)                            # slot -> gates in order
    reader = {}                                          # slot -> gate reading it as an input
    readers_all = defaultdict(list)
    writer = {}
    for n, ins, outs in gates:
        for q in set(ins + outs): touch[q].append(n)
        if args[n]:
            for q in ins:
                readers_all[q].append(n)
                if q not in reader: reader[q] = n          # first gate reading q (its garbage flows from there)
        for q in outs:
            if q not in writer: writer[q] = n
    pout = kk['pout']; rout = kk['rout']; src = kk['src']; srcslot = {q: S for S, q in src.items()}
    pslot = {q: c.pieces[i][0] for i, q in pout.items()}
    retslot = {q: nm for nm, q in rout.items()}
    # retained scatter: which y_S read which retained names
    def sret_names(S):
        if allE: return [('E', i) for i in range(h)]
        return [('*',)] + ([('E', i) for i in S] if h - 1 not in S else [('E', i) for i in range(h - 1) if i not in S])
    readers_of = defaultdict(list)
    for S in T:
        for nm in sret_names(S): readers_of[nm].append(S)
    # roots reachable from every node: pieces S and retained names
    reach = {}; vleafs_box = [set()]
    for n in sorted(c.active, reverse=True):
        r = set(('p', S) for S in npieces[n])
        if n in retnode: r.add(('r', retnode[n]))
        for g in gusers[n]: r |= reach[g]
        reach[n] = frozenset(r)
    def targets(q):
        if q in srcslot and not (vleaf and (c.tid[srcslot[q]] + 1) in vleafs_box[0]):
            rr = reach[c.tid[srcslot[q]] + 1]                # the leaf gate fans the source slot out
        elif q in reader:
            rr = set()
            for g in readers_all[q]: rr |= reach[g]      # a linked slot is read by several gates
        elif q in pslot: rr = {('p', pslot[q])}
        elif q in retslot: rr = {('r', retslot[q])}
        else: rr = set()
        out = set()
        for kind, x in rr:
            if kind == 'p': out.add(x)
            else: out.update(readers_of[x])
        return out
    # cones of the retained totals
    cone = set(); st = [n for _, n in c.retained]
    while st:
        x = st.pop()
        if x in cone: continue
        cone.add(x)
        if args[x]: st.extend(args[x])
    if clos:
        # phase 1 = closure of the cone gates under slot-sequence precedence in the time order (links may carry a
        # slot from a later-phase donor into a cone gate; that donor then joins phase 1)
        pos = {n: i for i, (n, _, _) in enumerate(gates)}
        st = list(cone); P = set(cone)
        while st:
            g = st.pop()
            ins, outs = gate_of[g]
            for q in set(ins + outs):
                for g2 in touch[q]:
                    if pos[g2] < pos[g] and g2 not in P: P.add(g2); st.append(g2)
        if verbose: print('phase-1 closure', len(cone), '->', len(P), flush=True)
        cone = P
    # late copies: the fan-out copy of a cone node n for a non-cone user is made in phase two from the pivot.
    # (the pivot outs[0] goes to users[n][0]; compile_roles orders users gate-users first, so if n has a cone
    #  gate user the pivot is consumed in phase 1; the copy is then made from a dedicated holder: we instead
    #  require the pivot to stay at M_n until the copy, which holds since the pivot's next touch is its consumer.)
    latecp = set()
    if late:
        for n, ins, outs in gates:
            if n not in cone or (vleaf and not args[n]): continue
            p0 = outs[0]
            if not ((reader.get(p0) is not None and reader[p0] not in cone) or p0 in pslot): continue
            for q in outs[1:]:
                g = reader.get(q)
                if (g is not None and g not in cone) or (q in pslot):
                    latecp.add(q)
    # V leaves: every use slot of a leaf gets x_t by its own gate on the data wire X_t (no fan-out gate), at
    # s_i = N_i cap N_{i+1} cap ... (late copies on X_t; uses ordered phase-1 first, then by dim N), so the
    # slots of phase-two users are untouched by phase 1.
    vls = {}; vorder = {}
    vleafs_box[0] = set()
    if vleaf:
        for n, ins, outs in gates:
            if args[n]: continue
            t = T[n - 1]
            us = []
            for q in outs:
                Nq = M[reader[q]] if q in reader else tperp[pslot[q]]
                ph = 1 if (q in reader and reader[q] in cone) else 2
                us.append((ph, len(Nq), q, Nq))
            us.sort()
            cur = line[t]; caps = [None] * len(us); acc = X.F
            for i in range(len(us) - 1, -1, -1):
                acc = X.meet(acc, us[i][3]); caps[i] = acc
            seq = []; got = {}
            for i, (ph, _, q, Nq) in enumerate(us):
                cap = caps[i]
                assert X.le(cur, cap)
                okf = lambda Y, Nq=Nq, cur=cur: X.nd(Y) and ok_res(list(Y)) and X.step_ok(Y, Nq) and X.step_ok(cur, Y)
                Y = X.fit(cur, cap, okf)
                if Y is None: break
                got[q] = Y; cur = Y; seq.append(q)
            else:
                vls.update(got); vorder[n] = seq      # every use placed: this leaf uses V gates
    vleafs = set(vorder); vleafs_box[0] = vleafs
    leafout = set(vls)
    # climbing late copies: the pivot (holder) of a node goes to its last user and climbs through
    # C_i = N_i cap N_{i+1} cap ... ; the copy for user i is made at C_i just before that user's gate.
    climbcp = {}; cbefore = defaultdict(list); cpieces = []; nclimb = Counter()
    if climb:
        for n, ins, outs in gates:
            if not args[n] or len(outs) < 2: continue
            if any(q in retslot for q in outs): nclimb['retained'] += 1; continue
            seq = list(outs[1:]) + [outs[0]]                 # time order of the users (holder's user last)
            def user_frame(q):
                rl_ = [g for g in readers_all.get(q, ()) if g != n]
                if rl_: return M[rl_[0]]
                return tperp[pslot[q]]
            Ns = [user_frame(q) for q in seq]
            caps = [None] * len(seq); acc = X.F
            for i in range(len(seq) - 1, -1, -1):
                acc = X.meet(acc, Ns[i]); caps[i] = acc
            cur = M[n]; got = {}
            for i, q in enumerate(seq[:-1]):
                Nq = Ns[i]
                okf = lambda Y, Nq=Nq, cur=cur: X.nd(Y) and ok_res(list(Y)) and X.step_ok(Y, Nq) and X.step_ok(cur, Y)
                Y = X.fit(cur, caps[i], okf)
                if Y is None: break
                got[q] = Y; cur = Y
            else:
                if X.step_ok(cur, Ns[-1]):
                    for q, Y in got.items():
                        climbcp[q] = (outs[0], Y)
                        if q in reader: cbefore[reader[q]].append(q)
                        else: cpieces.append(q)
                    nclimb['climbed'] += 1; continue
            nclimb['failed'] += 1
        latecp = set(climbcp)
        if verbose: print('climb', dict(nclimb), flush=True)
    def phase1_touch(q):
        if q in climbcp: return reader.get(q) in cone
        if q in leafout: return q in reader and reader[q] in cone
        return any(g in cone and not (q in latecp and g == writer.get(q)) for g in touch[q])
    free = [q for q in range(R) if q not in retslot and not phase1_touch(q)] if defer else []
    freeset = set(free)
    # first pass-two frame of every slot
    def F0(q):
        if q in vls: return vls[q]
        if q in climbcp: return climbcp[q][1]
        if q in srcslot: return line[srcslot[q]]
        return M[touch[q][0]]
    # ---------------- hosting: choose sigma_q for free slots ----------------
    sigma = {}
    chain = {S: [] for S in T}                           # deferred frames on y_S (sorted by dim), excluding 0, t^perp
    tg = {q: targets(q) for q in free}
    lamf = lambda r: r * math.log(m / r) if r else 0.0
    def slot_gain(f0, s_):
        if s_ == 0: return 0.0
        return lamf(f0) + lamf(m - h) - (lamf(f0 - s_) + lamf(m - h + s_))
    singles = [q for q in free if len(tg[q]) == 1] if dp else []
    sset = set(singles)
    if prio: order = sorted([q for q in free if q not in sset], key=lambda q: -len(F0(q)) ** 2 / max(1, len(tg[q])))
    else: order = sorted([q for q in free if q not in sset], key=lambda q: -len(F0(q)))
    nd_fail = Counter()
    import random; hrnd = random.Random(7)
    for q in order:
        X0 = F0(q); Xc = X0
        for _ in range(50):
            changed = False
            for S in tg[q]:
                ch_ = chain[S]
                # largest chain element <= Xc, then cap with the next one
                nxt = tperp[S]
                for C_ in ch_:
                    if not X.le(C_, Xc):
                        nxt = C_; break
                Y = X.meet(Xc, nxt) if not X.le(Xc, nxt) else Xc
                if Y != Xc: Xc = Y; changed = True
            if not changed: break
        # sigma must also contain every chain element below it (comparable): ensured by construction when stable
        s_ = Xc
        def valid(s_):
            if not s_: return True
            if not X.nd(s_): nd_fail['r:nd'] += 1; return False
            if not X.step_ok(s_, X0): nd_fail['r:slotstep'] += 1; return False
            for S in tg[q]:
                ch_ = chain[S]
                below = [C_ for C_ in ch_ if len(C_) <= len(s_)]
                above = [C_ for C_ in ch_ if len(C_) > len(s_)]
                if below and below[-1] == s_: continue
                if not all(X.le(C_, s_) for C_ in below) or not all(X.le(s_, C_) for C_ in above):
                    nd_fail['r:chain'] += 1; return False
                lo = below[-1] if below else ()
                hi = above[0] if above else tperp[S]
                if not X.step_ok(lo, s_): nd_fail['r:ylo'] += 1; return False
                if not X.step_ok(s_, hi): nd_fail['r:yhi'] += 1; return False
            return True
        got = ()
        if s_:
            Bl = []
            for S in tg[q]:
                bl = [C_ for C_ in chain[S] if X.le(C_, s_)]
                if bl: Bl += list(bl[-1])
            Bl = X.canon(Bl)
            got = X.fit(Bl, s_, lambda Y: bool(Y) and valid(Y)) or ()
            if not got:
                pool = set()
                for S in tg[q]:
                    for C_ in chain[S]:
                        if X.le(C_, s_): pool.add(C_)
                for C_ in sorted(pool, key=len, reverse=True):
                    if valid(C_): got = C_; break
            if not got: nd_fail['invalid'] += 1
        s_ = got
        if not s_: nd_fail['zero'] += 1; continue
        sigma[q] = s_
        for S in tg[q]:
            ch_ = chain[S]
            if s_ not in ch_:
                ch_.append(s_); ch_.sort(key=len)
    if dp:
        # single-target slots: per target, choose the chain of new frames by DP over candidates (their first frames
        # and pairwise meets), keeping the frames already placed by multi-target slots. A chain C_1 < ... < C_k gives
        # slot u the largest C_i <= X0(u); its value is sum_i sum_{u: C_i <= X0(u)} (gain_u(C_i) - gain_u(C_{i-1})).
        byS = defaultdict(list)
        for q in singles: byS[next(iter(tg[q]))].append(q)
        for S, qs in byS.items():
            E = chain[S]                                     # fixed elements
            X0s = {q: F0(q) for q in qs}
            cand = set(X0s.values())
            vals = list(cand)
            for i in range(len(vals)):
                for j in range(i + 1, len(vals)):
                    mt = X.meet(vals[i], vals[j])
                    if mt: cand.add(mt)
            for e in E:
                for x in vals:
                    mt = X.meet(e, x)
                    if mt: cand.add(mt)
            cand |= set(E)
            okc = [C_ for C_ in cand if X.nd(C_) and ok_res(list(C_)) and X.le(C_, tperp[S])
                   and all(X.le(C_, e) or X.le(e, C_) for e in E)]
            okc.sort(key=len)
            Eset = set(E); Edims = [len(e) for e in E]
            cover = {C_: [q for q in qs if X.le(C_, X0s[q]) and X.step_ok(C_, X0s[q])] for C_ in okc}
            gq = {q: (lambda d, f0=len(X0s[q]): slot_gain(f0, d)) for q in qs}
            f = {(): 0.0}; prev = {(): None}; nodes = [()] + okc
            for C_ in okc:
                best = None
                for Cp in nodes:
                    if len(Cp) >= len(C_): break
                    if Cp not in f: continue
                    if not (X.le(Cp, C_) and X.step_ok(Cp, C_)): continue
                    if any(len(Cp) < de < len(C_) for de in Edims): continue      # must pass through fixed ones
                    if Cp == () and any(de < len(C_) for de in Edims): continue
                    w = sum(gq[q](len(C_)) - gq[q](len(Cp)) for q in cover[C_])
                    if best is None or f[Cp] + w > best[0]: best = (f[Cp] + w, Cp)
                if best is not None: f[C_] = best[0]; prev[C_] = best[1]
            ends = [C_ for C_ in f if X.step_ok(C_, tperp[S]) and not any(de > len(C_) for de in Edims)
                    and (not E or C_ != () )] + ([()] if not E else [])
            if not ends: continue
            top = max(ends, key=lambda C_: f[C_])
            ch_new = []
            while top: ch_new.append(top); top = prev[top]
            ch_new.reverse()
            chain[S] = ch_new
            for q in qs:
                bestC = None
                for C_ in ch_new:
                    if X.le(C_, X0s[q]) and X.nd(C_) and X.step_ok(C_, X0s[q]): bestC = C_
                if bestC: sigma[q] = bestC
                else: nd_fail['dp-zero'] += 1
    if verbose and defer:
        print('free', len(free), 'deferred', len(sigma), dict(nd_fail),
              'sigma dims', sorted(Counter(len(s_) for s_ in sigma.values()).items()), flush=True)
    # ---------------- the stage-one word ----------------
    ops = []      # (roles, frame, kind)
    S_ = lambda q: ('s', q)
    # garbage readouts at 0 for non-deferred slots (direct readouts; pass 1 at frame 0)
    for q in range(R):
        if q in sigma: continue
        ops.append((tuple([S_(q)] + [('y', t) for t in sorted(targets(q))]), (), 'read'))
    # pass two, phase 1: cone leaves and cone gates (late copies excluded)
    def gate_op(n, excl=()):
        ins, outs = gate_of[n]
        rs = [q for q in dict.fromkeys(ins + outs) if q not in excl]
        return (tuple(S_(q) for q in rs), M[n], 'gate', (n, tuple(q for q in outs[1:] if q in excl)))
    def leaf_copy(n):
        S = T[n - 1]; return ((('x', S), S_(src[S])), line[S], 'xcopy')
    p1 = [n for n, _, _ in gates if n in cone] if defer else []
    vgate = lambda q, n, F: ((('x', T[n - 1]), S_(q)), F, 'vgate')
    for n in p1:
        if not args[n]:
            if n in vleafs:
                for q in vorder[n]:
                    if phase1_touch(q): ops.append(vgate(q, n, vls[q]))
                continue
            ops.append(leaf_copy(n))
        for q in cbefore.get(n, ()): ops.append(((S_(climbcp[q][0]), S_(q)), climbcp[q][1], 'latecopy'))
        ops.append(gate_op(n, latecp))
    # centre reads (copied retained totals, at 0)
    cread_at = len(ops)
    ops.append((tuple(S_(q) for q in sorted(retslot)) + tuple(('y', t) for t in T), (), 'cread'))
    # deferred readouts, sorted by dim sigma
    for q in sorted(sigma, key=lambda q: len(sigma[q])):
        ops.append((tuple([S_(q)] + [('y', t) for t in sorted(tg[q])]), sigma[q], 'read'))
    # late copies of cone nodes (pivot and copy slot at M_n)
    for q in sorted(latecp - set(climbcp)):
        n = writer[q]; ins, outs = gate_of[n]
        ops.append(((S_(outs[0]), S_(q)), M[n], 'latecopy'))
    # V gates of phase-two uses (after their deferred readouts), then the phase-two gates
    if vleaf:
        for n, seq in vorder.items():
            for q in seq:
                if not (defer and n in cone and q in reader and reader[q] in cone): ops.append(vgate(q, n, vls[q]))
    for n, ins, outs in gates:
        if n in cone and defer: continue
        if not args[n]:
            if n in vleafs: continue
            ops.append(leaf_copy(n))
        for q in cbefore.get(n, ()): ops.append(((S_(climbcp[q][0]), S_(q)), climbcp[q][1], 'latecopy'))
        ops.append(gate_op(n, latecp))
    if not defer:          # the centre reads stay after the cone (base word: after all of pass two's mix)
        op = ops.pop(cread_at); ops.append(op)
    # injects at t^perp, uncompute at F, source copies at F
    for q in cpieces: ops.append(((S_(climbcp[q][0]), S_(q)), climbcp[q][1], 'latecopy'))
    for i, q in pout.items():
        S = c.pieces[i][0]; ops.append(((('y', S), S_(q)), tperp[S], 'inject', i))
    for S in T: ops.append(((('y', S),), tperp[S], 'inject'))      # every y_S reaches t^perp
    ops.append((tuple(S_(q) for q in range(R)), X.F, 'high'))
    for n, seq in vorder.items():
        for q in seq: ops.append(vgate(q, n, X.F))
    for S in T:
        if c.tid[S] + 1 not in vleafs: ops.append(((('x', S), S_(src[S])), X.F, 'xcopy'))
    return dict(c=c, kk=kk, X=X, ops=ops, sigma=sigma, M=M, m=m, v=v, R=R, h=h, T=T, retslot=retslot, cone=cone,
                latecp=latecp, tg=tg, targets=targets, vls=vls, vorder=vorder)


def walk(B, check=True):
    """role frame chains, histogram, validity."""
    X = B['X']; h = B['h']; m = B['m']; v = B['v']; T = B['T']; R = B['R']; ops = B['ops']
    fr = defaultdict(list); start = {}
    Hs = Counter(); copies = Counter(); bad = Counter()
    for op in ops:
        roles, F, kind = op[:3]
        if kind == 'cread':
            for r in roles:
                if r[0] == 's':
                    p = fr[r][-1]; copies[len(p)] += 1      # copy made at p, descends to 0, read, erased
                    if check and p and not X.step_ok(p, ()): bad['copy'] += 1
                else:
                    fr[r].append(())
            continue
        for r in roles: fr[r].append(F)
    # auxiliary slots: chain from start = first frame, wrap child F_last -> start + full
    for q in range(R):
        ch = fr[('s', q)]
        assert ch and ch[-1] == X.F
        for a, b in zip(ch, ch[1:]):
            if a != b:
                Hs[abs(len(b) - len(a))] += 1
                if check and not X.step_ok(a, b): bad['slot'] += 1
        w = m - h + len(ch[0]); Hs[w] += 1
        if check and not X.nd(ch[0]): bad['wrap'] += 1
    for k_, n_ in copies.items():
        if k_: Hs[k_] += n_
    # data chains (stage one); y starts at 0, x at the line
    yc = {}; xc = {}
    for S in T:
        ys = [()] + fr[('y', S)]; xs = [X.canon([sum(1 << p for p in S)])] + fr[('x', S)]
        for nm, ch in (('y', ys), ('x', xs)):
            for a, b in zip(ch, ch[1:]):
                if check and a != b and not X.step_ok(a, b): bad[nm] += 1
        assert ys[-1] == X.perp(X.canon([sum(1 << p for p in S)])) and xs[-1] == X.F
        yc[S] = [len(f) for f in ys]; xc[S] = [len(f) for f in xs]
    H = Counter()
    for k_, n_ in Hs.items(): H[k_] += 2 * v * n_           # both stages (complement time reversal)
    # stage-two mirrors: target (X wire) gets the mirror of the y chain, source (Y wire) the mirror of the x chain
    mir = lambda ch: [m - d for d in reversed(ch)]
    p1 = Counter((tuple(xc[S]), tuple(yc[S])) for S in T)
    p2 = Counter((tuple(mir(yc[S])), tuple(mir(xc[S]))) for S in T)
    # each pair's chain = (stage-one part) + (stage-two part); steps inside a part count v times, the junction
    # step depends only on the two boundary dims
    def inner(ch, mult):
        for a, b in zip(ch, ch[1:]):
            if a != b: H[abs(b - a)] += mult
    n2tot = sum(p2.values()); n1tot = sum(p1.values())
    J = [Counter(), Counter()]; K = [Counter(), Counter()]
    for (x1, y1), n1 in p1.items():
        inner(list(x1), n1 * n2tot); inner([0] + list(y1), n1 * n2tot)
        J[0][x1[-1]] += n1; J[1][y1[-1]] += n1
    for (x2, y2), n2 in p2.items():
        inner(list(x2) + [m], n2 * n1tot); inner(list(y2) + [m - 1], n2 * n1tot)
        K[0][x2[0]] += n2; K[1][y2[0]] += n2
    for k in (0, 1):
        for a, na in J[k].items():
            for b, nb in K[k].items():
                if a != b: H[abs(b - a)] += na * nb
    N = v * v; H[1] += N
    Wn = 2 * N + 2 * v * R; L = 2 * v * sum(k_ * n_ for k_, n_ in copies.items())
    s = Wn * m - N + L; tot = sum(k_ * n_ for k_, n_ in H.items())
    return dict(H={k_: n_ for k_, n_ in sorted(H.items()) if n_}, W=Wn, s=s, sum=tot, sum_ok=tot == s, bad=dict(bad), L=L,
                maxrank=max(H))


if __name__ == '__main__':
    from cert import cert
    h = int(sys.argv[1]); fl = set(sys.argv[2:])
    B = build(h, allE='allE' in fl, lift='lift' in fl, defer='defer' in fl, late='late' in fl)
    w = walk(B, check='nocheck' not in fl)
    m = B['m']; H = w['H']
    D = w['W'] * m - w['s']; lam = sum(n * r * math.log(m / r) for r, n in H.items())
    print(dict(h=h, flags=sorted(fl), R=B['R'], W=w['W'], s=w['s'], sum_ok=w['sum_ok'], bad=w['bad'], maxrank=w['maxrank'],
               est=D / lam), flush=True)
    if 'cert' in fl:
        a, gap = cert(H, w['W'], m); print(dict(a_c=str(a), a_c_float=float(a), gap=float(gap)), flush=True)
    if 'dump' in fl:
        pickle.dump(dict(H=H, W=w['W'], s=w['s'], m=m), open('h7_%d_%s.pkl' % (h, '_'.join(sorted(fl - {'dump', 'cert'}))), 'wb'))
