"""Signed dependence reclamation (the mechanism of upstream PR #112, read as data and reimplemented) on our complex
two-stage word (c7.build without lift/defer/vleaf: direct garbage readouts at frame 0, pass two at node frames).

Pass two (the middle mixer L) is replayed symbolically: every auxiliary slot carries an exact signal, a set of
triples with coefficient +1 (the producer never cancels). When a non-terminal slot q is read for the last time
(it "dies" at op d, frame E = frame of op d), its signal b may be CLEARED by exact signed shears from controls:
  dup   b is held by another slot (or b = x_t: read from the data wire X_t, a V event)
  diff  b = n - e for a user n = b + e of b, with n and e held (or e = x_t from X_t)
  sum   b = c + d for b's arguments (c, d), each held (or a leaf read from its data wire)
Each control r is read at E, so its frame chain must allow it: now(r) <= E <= next(r), both steps valid residuals;
a data wire X_t needs now(X_t) <= E (its next frame is F). The cleared slot (signal 0, arbitrary scratch kept, the
shears are invertible on the full dirty space) waits at E and is REUSED by a later fresh slot (a fan-out copy or a
piece/terminal copy) born at op b > d with frame f_b >= E (valid step E -> f_b). Each reuse removes one auxiliary.
The pass-two shears are part of L, so 'high' (L^-1 at F) reverts them too; the garbage readouts G = A L (x = 0)
are recomputed from the actual word. Net effect on the targets: A (L V + K) x where K are the V clears, i.e. every
terminal register holds its exact terminal signal (checked symbolically here, exactly over Q(i) in e2e_rc.py).
Usage: python3 reclaim.py H cfg [cert]   (cfg = comma flags for c7.build: allE, dual, alt, links, pasm, ivec ...)"""
import sys, os, math, bisect
from collections import defaultdict, Counter
LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib')
sys.path.insert(0, LIB)
import c7, e2e7
DEEP = 6
VJMAX = 64
VJ_EAGER = False


def parse_cfg(cfg):
    fl = set() if cfg in ('-', '') else set(cfg.split(','))
    ps = (6 if 'ival+' in fl else 5 if 'ival' in fl else 4 if 'dual+' in fl else 3 if 'dual' in fl else 2 if 'pstar2' in fl
          else 1 if 'pstar' in fl else 0)
    kw = dict(allE='allE' in fl, pstar=ps, alt='alt' in fl, links='links' in fl, pasm='pasm' in fl, ivec='ivec' in fl,
              lift='lift' in fl, vleaf='vleaf' in fl, defer='defer' in fl, clos='clos' in fl)
    if 'dag117' in fl: kw['dag'] = os.path.join(LIB, '..', '..', '..', 'certificates', 'round8', 'pr117_dag.json.gz')
    return fl, kw


def base(h, cfg):
    fl, kw = parse_cfg(cfg)
    import links as _lk; _lk.PIVOT_SMALL_DEAD[0] = 'pivsmall' in fl
    B = c7.build(h, verbose=False, **kw)
    W = e2e7.updates(B)
    return B, W


def reclaim(B, W, use_v=True, rels=('dup', 'diff', 'sum'), pick='maxdim', verbose=True, diag=None, lazy=0, desc=0):
    dcache = {}; premove = defaultdict(list)
    c = B['c']; X = B['X']; args = c.args; R = B['R']; T = B['T']
    kk = B['kk']; pslots = set(kk['pout'].values()); rslots = set(kk['rout'].values())
    term = pslots | rslots
    i0 = next(i for i, w in enumerate(W) if w[2] != 'read')
    ih = next(i for i, w in enumerate(W) if w[2] == 'high')
    F = X.F
    leafnode = {1 << i: i + 1 for i in range(len(T))}
    if not hasattr(c, 'look'): c.look = {c.sup[n]: n for n in c.active if c.args[n]}
    def node_of(sup):
        return leafnode.get(sup) or c.look.get(sup)
    gusers = defaultdict(list)
    for n in c.active:
        if args[n]:
            for a in args[n]: gusers[a].append(n)
    # per-role event frames in pass two (times are floats so clears can be inserted at i + 0.5)
    evt = defaultdict(list); evf = defaultdict(list)
    for i in range(i0, ih):
        roles, lab, kind, ups = W[i]
        for r in roles:
            if kind == 'cread' and r[0] == 's': continue
            evt[r].append(float(i)); evf[r].append(lab)
    last = {}; first_write = {}
    for i in range(i0, ih):
        roles, lab, kind, ups = W[i]
        if kind == 'cread': continue
        for r in roles:
            if r[0] == 's': last[r[1]] = i
    def now_next(r, t):
        ts = evt[r]; k = bisect.bisect_right(ts, t)
        now = evf[r][k - 1] if k else (B['X'].canon([sum(1 << p for p in r[1])]) if r[0] == 'x' else ())
        nxt = evf[r][k] if k < len(ts) else F
        return now, nxt
    def can_read(r, t, E):
        now, nxt = now_next(r, t)
        if not (X.le(now, E) and X.le(E, nxt)): return False
        return X.step_ok(now, E) and X.step_ok(E, nxt)
    def insert(r, t, E):
        ts = evt[r]; k = bisect.bisect_right(ts, t); ts.insert(k, t); evf[r].insert(k, E)
    sig = defaultdict(int)                                    # slot -> support (0 = fresh)
    holders = defaultdict(set)                                # support -> slots
    def setsig(q, s_):
        o = sig[q]
        if o: holders[o].discard(q)
        sig[q] = s_
        if s_: holders[s_].add(q)
    pool = defaultdict(list)                                  # frame -> [(phys, t_clear)]
    ccache = {}
    def compat(E, f):
        k = (E, f)
        if k not in ccache: ccache[k] = X.le(E, f) and X.step_ok(E, f)
        return ccache[k]
    phys = list(range(R)); clears = defaultdict(list); pre = defaultdict(list)          # op index -> list of clear ops
    stat = Counter(); retired = set()
    def find_ctrl(sup, t, E, excl):
        """a control for signal sup at time t, frame E: ('s', slot) or ('x', triple)."""
        for q in holders.get(sup, ()):
            if q in excl or q in rslots or q in retired: continue
            if can_read(('s', q), t, E): return ('s', q)
        if use_v and sup in leafnode:
            S = T[leafnode[sup] - 1]
            if can_read(('x', S), t, E): return ('x', S)
        return None
    def plan_for(q, t, E):
        b = sig[q]; nb = node_of(b)
        if 'dup' in rels:
            ctl = find_ctrl(b, t, E, {q})
            if ctl: return [(ctl, -1)], 'dup'
        if 'diff' in rels and nb:
            for n in gusers.get(nb, ()):
                if n not in c.active: continue
                a0, a1 = args[n]; e = a1 if a0 == nb else a0
                c1 = find_ctrl(c.sup[n], t, E, {q})
                if not c1: continue
                c2 = find_ctrl(c.sup[e], t, E, {q, c1[1]} if c1[0] == 's' else {q})
                if c2: return [(c2, 1), (c1, -1)], 'diff'
        if 'sum' in rels and nb and args[nb]:
            cc, dd = args[nb]
            c1 = find_ctrl(c.sup[cc], t, E, {q})
            if c1:
                c2 = find_ctrl(c.sup[dd], t, E, {q, c1[1]} if c1[0] == 's' else {q})
                if c2: return [(c1, -1), (c2, -1)], 'sum'
        if 'deep' in rels and nb and args[nb]:
            # b as a disjoint sum of held controls and data-wire leaves, down the producer DAG (depth <= DEEP)
            excl = {q}; out = []
            def dec(n, d):
                if d < DEEP:
                    ct = find_ctrl(c.sup[n], t, E, excl)
                    if ct:
                        if ct[0] == 's': excl.add(ct[1])
                        out.append((ct, -1)); return True
                if not args[n] or d == 0: return False
                return dec(args[n][0], d - 1) and dec(args[n][1], d - 1)
            if dec(nb, DEEP): return out, 'deep'
        if 'vjoin' in rels and E is not None:
            # pure data-wire clear at the join of E and every leaf wire's current frame (raised, never lowered)
            leaves = []; bb = b
            while bb:
                lo = bb & -bb; leaves.append(T[lo.bit_length() - 1]); bb ^= lo
            if len(leaves) <= VJMAX:
                J = list(E)
                for S in leaves: J += list(now_next(('x', S), t)[0])
                J = X.canon(J)
                if len(J) < len(F) and (J == E or X.step_ok(now_next(('s', q), t)[0], J)) and \
                        all(can_read(('x', S), t, J) for S in leaves):
                    return [(('x', S), -1) for S in leaves], 'vjoin', J
        return None, None
    def do_clear(q, plan, t, E, where):
        for ctl, cf in plan:
            insert(ctl, t, E)
            if ctl[0] == 's' and ctl[1] in dkey:
                dead[dkey[ctl[1]]].discard(ctl[1]); dead[E].add(ctl[1]); dkey[ctl[1]] = E
        if q in dkey: dead[dkey.pop(q)].discard(q)
        deadset.discard(q)
        insert(('s', q), t, E)
        where.append((q, plan, E))
        setsig(q, 0); retired.add(q)
    dead = defaultdict(set)                                   # frame -> dead dirty slots
    deadset = set(); dkey = {}; used = set()
    for i in range(i0, ih):
        roles, lab, kind, ups = W[i]
        # births at op i: fresh slots written from zero
        if kind != 'cread':
            born = []
            for d_, s_, cf in ups or ():
                if d_[0] == 's' and sig[d_[1]] == 0 and d_[1] not in born and (d_[1] not in first_write):
                    born.append(d_[1])
            for q in born:
                first_write[q] = i
                # a V injection (xcopy / vgate) must land on a fresh register: V is outside L, so L^-1 cannot
                # undo it, and a reused register's earlier L updates would be reverted against the wrong value
                if kind in ('xcopy', 'vgate'): stat['new'] += 1; continue
                best = None
                for E in list(pool):
                    ent = pool[E]
                    if not ent or ent[0][1] >= i or not compat(E, lab): continue
                    if best is None or len(E) > len(best): best = E
                if best is not None:
                    p, tc, qc = pool[best].pop(0); used.add(qc)
                    if not pool[best]: del pool[best]
                    phys[q] = p; stat['reused'] += 1; continue
                if lazy:
                    t = i - 0.25; tried = 0; done = False
                    for Eq in list(dead):
                        if not dead[Eq] or not compat(Eq, lab): continue
                        for q2 in list(dead[Eq]):
                            tried += 1
                            res = plan_for(q2, t, lab); plan, why = res[0], res[1]
                            if len(res) > 2 and res[2] != lab: plan = None
                            if plan:
                                do_clear(q2, plan, t, lab, pre[i])
                                used.add(q2); phys[q] = phys[q2]; stat['lazy_' + why] += 1; stat['reused'] += 1; done = True; break
                            if tried > lazy: break
                        if done or tried > lazy: break
                    if done: continue
                stat['new'] += 1
        # apply op i symbolically
        if kind != 'cread':
            for d_, s_, cf in ups or ():
                if d_[0] != 's': continue
                sv = (1 << c.tid[s_[1]]) if s_[0] == 'x' else sig[s_[1]]
                assert cf == 1 and not (sig[d_[1]] & sv), 'unexpected update'
                setsig(d_[1], sig[d_[1]] | sv)
        # deaths at op i
        for r in roles:
            if r[0] != 's' or kind == 'cread': continue
            q = r[1]
            if last.get(q) != i or q in term or q in retired: continue
            b = sig[q]; t = i + 0.5; E = lab; nb = node_of(b)
            res = plan_for(q, t, E); plan, why = res[0], res[1]
            if len(res) > 2 and res[2] != E and not VJ_EAGER: plan = None; why = None
            elif len(res) > 2: E = res[2]
            if plan: stat[why] += 1
            if not plan:
                plan = None
            if plan is None:
                stat['stuck'] += 1
                if diag is not None:
                    hs = [x for x in holders.get(b, ()) if x != q]
                    nu = sum(1 for n in gusers.get(nb, ()) if n in c.active and holders.get(c.sup[n]))
                    why = []
                    for x in hs:
                        nw, nx = now_next(('s', x), t)
                        why.append(('now<=E' if X.le(nw, E) else 'now!<=E', 'E<=next' if X.le(E, nx) else 'E!<=next(%d vs %d)' % (len(E), len(nx)),
                                    'dead' if x in deadset else 'alive', 'term' if x in term else ''))
                    diag[('sz', min(bin(b).count('1'), 9), tuple(sorted(set(why)))[:2])] += 1
                deadset.add(q)
                if lazy: dead[lab].add(q); dkey[q] = lab
                continue
            do_clear(q, plan, t, E, clears[i])
            pool[E].append((phys[q], i, q))
    # physical chains: resolve phys through reuse chains
    def root(q):
        while phys[q] != q: q = phys[q]
        return q
    rmap = {}; ids = {}
    for q in range(R):
        rq = root(q)
        if rq not in ids: ids[rq] = len(ids)
        rmap[q] = ids[rq]
    Rn = len(ids)
    if verbose: print('reclaim', dict(stat), 'R', R, '->', Rn, flush=True)
    roots = {rmap[q] : q for q in range(R) if phys[q] == q}
    return dict(clears=clears, pre=pre, premove=premove, rmap=rmap, roots=roots, used=used, R=Rn, stat=stat, i0=i0, ih=ih)


def assemble(B, W, rc, h_targets=False):
    """the full stage-one word with clears inserted and slots renamed; explicit ups on every op.
    Returns a B-like dict for c7.walk and e2e_rc (ops = [roles, frame, kind, ups])."""
    rmap = rc['rmap']; Rn = rc['R']; T = B['T']; X = B['X']
    ren = lambda r: ('s', rmap[r[1]]) if r[0] == 's' else r
    rup = lambda ups: [(ren(d), ren(s_), cf) for d, s_, cf in ups] if ups else ups
    ops = []
    roots = set(rc['roots'].values())
    # garbage readouts: each physical slot is read once, by its first logical slot's readout (pass one at frame 0,
    # or a deferred readout in pass two); readouts of later logical slots of the same register are dropped
    def keep_read(roles):
        return roles[0][1] in roots
    for roles, lab, kind, ups in W[:rc['i0']]:
        assert kind == 'read'
        if keep_read(roles):
            ops.append([(ren(roles[0]),) + (tuple(('y', S) for S in T) if h_targets and not lab else tuple(roles[1:])), lab, 'read', None])
    p2 = []; xs_used = set()
    for i in range(rc['i0'], rc['ih']):
        roles, lab, kind, ups = W[i]
        for q, Mt in rc['premove'].get(i, ()):
            ops.append([(('s', rmap[q]),), Mt, 'move', []])
        for q, plan, E in rc['pre'].get(i, ()):
            cu = [(('s', rmap[q]), ren(ctl), cf) for ctl, cf in plan]
            croles = tuple(dict.fromkeys([('s', rmap[q])] + [ren(ctl) for ctl, _ in plan]))
            for ctl, _ in plan:
                if ctl[0] == 'x': xs_used.add(ctl)
            ops.append([croles, E, 'shear', cu]); p2 += cu
        if kind == 'read':
            if keep_read(roles): ops.append([tuple(ren(r) for r in roles), lab, 'read', None])
            continue
        rr = tuple(dict.fromkeys(ren(r) for r in roles))
        u = rup(ups)
        ops.append([rr, lab, kind, u])
        if kind in ('gate', 'latecopy', 'xcopy', 'vgate') and u:
            if kind in ('xcopy', 'vgate'): pass          # V (+): undone by the final -V copies
            else: p2 += u
        for q, plan, E in rc['clears'].get(i, ()):
            if q not in rc['used']: continue           # a clear whose register is never reused is not needed
            cu = [(('s', rmap[q]), ren(ctl), cf) for ctl, cf in plan]
            croles = tuple(dict.fromkeys([('s', rmap[q])] + [ren(ctl) for ctl, _ in plan]))
            for ctl, _ in plan:
                if ctl[0] == 'x': xs_used.add(ctl)
            ops.append([croles, E, 'shear', cu]); p2 += cu
    # injects and the inject-all, then high (L^-1 at F, reverting gates and clears), then -V copies at F
    tail = W[rc['ih'] + 1:]
    hi = W[rc['ih']]
    ops.append([tuple(('s', q) for q in range(Rn)) + tuple(sorted(xs_used)), X.F, 'high',
                [(d, s_, -cf) for d, s_, cf in reversed(p2)]])
    for roles, lab, kind, ups in tail:
        ops.append([tuple(ren(r) for r in roles), lab, kind, rup(ups)])
    out = ops
    Bn = dict(B); Bn['ops'] = out; Bn['R'] = Rn
    Bn['retslot'] = {rmap[q]: nm for q, nm in B['retslot'].items()}; Bn['rmap'] = rmap
    return Bn


def gfix(Bn, P=(1 << 61) - 1, verbose=True):
    """garbage map G = A L on the scratch (x = 0), mod a large prime (support only), from the assembled word; every
    readout then declares exactly its slot's support. A deferred readout whose support is not inside its original
    targets (a reused register carries garbage to more targets) is moved to pass one at frame 0."""
    from fractions import Fraction
    ops = Bn['ops']; inv2 = pow(2, P - 2, P)
    def md(cf):
        cf = Fraction(cf); return cf.numerator % P * pow(cf.denominator % P, P - 2, P) % P
    val = {}
    for roles, lab, kind, ups in ops:
        if kind == 'read' or not ups: continue
        if kind == 'high': break
        for d, s_, cf in ups:
            if s_[0] == 'x': continue
            src = val.get(s_)
            if src is None: src = {s_[1]: 1} if s_[0] == 's' else {}
            if d not in val: val[d] = {d[1]: 1} if d[0] == 's' else {}
            dv = val[d]; c_ = md(cf)
            for k, v in src.items():
                w = (dv.get(k, 0) + c_ * v) % P
                if w: dv[k] = w
                else: dv.pop(k, None)
    G = defaultdict(set)
    for r, dv in val.items():
        if r[0] == 'y':
            for q in dv: G[q].add(r[1])
    moved = 0; first = []; rest = []
    for o in ops:
        if o[2] != 'read': rest.append(o); continue
        q = o[0][0][1]; tg = set(r[1] for r in o[0][1:]) if o[1] else None
        if o[1] and not G[q] <= tg:
            moved += 1; o = [o[0][:1], (), 'read', None]
        if not o[1]:
            first.append([(('s', q),) + tuple(('y', S) for S in sorted(G[q])), (), 'read', None])
        else: rest.append(o)
    Bn['ops'] = first + rest
    if verbose: print('gfix', dict(moved_to_frame0=moved, reads=len(first) + sum(1 for o in rest if o[2] == 'read')), flush=True)
    return moved


def check_signals(Bn, dbg=None):
    """symbolic replay of the whole stage-one word with z = 0 (supports of +-1 sums of triples): every update is a
    +-1 shear between distinct roles that adds a disjoint support or removes a contained one; at the end of L every
    piece / retained slot holds exactly its terminal node; after 'high' every slot is back to signal 0."""
    c = Bn['c']; kk = Bn['kk']; rmap = Bn['rmap']
    val = defaultdict(int); bad = Counter(); chk = {}
    for roles, lab, kind, ups in Bn['ops']:
        if kind == 'high':
            want = {rmap[q]: c.sup[n] for nm, n in c.retained for q in [kk['rout'][nm]]}
            want.update({rmap[q]: c.sup[c.pieces[i][1]] for i, q in kk['pout'].items()})
            chk['terminal_bad'] = sum(val[q] != s_ for q, s_ in want.items()); chk['terminals'] = len(want)
        if kind in ('read', 'cread', 'inject') or not ups: continue
        for d_, s_, cf in ups:
            if d_ == s_: bad['self'] += 1
            if d_[0] == 's' and cf == -1 and dbg is not None and (val[d_[1]] & ((1 << c.tid[s_[1]]) if s_[0] == 'x' else val[s_[1]])) != ((1 << c.tid[s_[1]]) if s_[0] == 'x' else val[s_[1]]): dbg.append((kind, d_, s_, lab))
            if d_[0] != 's': continue
            sv = (1 << c.tid[s_[1]]) if s_[0] == 'x' else val[s_[1]]
            if cf == 1:
                if val[d_[1]] & sv: bad['overlap'] += 1
                val[d_[1]] |= sv
            elif cf == -1:
                if (val[d_[1]] & sv) != sv: bad['notsub'] += 1
                val[d_[1]] &= ~sv
            else: bad['coef'] += 1
        if kind == 'high': chk['nonzero_after_high'] = sum(1 for v in val.values() if v)
    chk.update(bad); return chk


if __name__ == '__main__':
    h = int(sys.argv[1]); cfg = sys.argv[2]; extra = set(sys.argv[3:])
    B, W = base(h, cfg)
    w0 = c7.walk(B, check=False)
    from cert import cert
    H0 = w0['H']; m0 = B['m']
    D0 = w0['W'] * m0 - w0['sum']; lam0 = sum(n * r * math.log(m0 / r) for r, n in H0.items())
    print('base', dict(R=B['R'], W=w0['W'], sum_ok=w0['sum_ok'], est=D0 / lam0), flush=True)
    if 'cert' in extra and 'basecert' in extra:
        a0, _ = cert(H0, w0['W'], m0); print('base_a_c', str(a0), float(a0), flush=True)
    lz = next((int(x[4:]) for x in extra if x.startswith('lazy')), 0)
    for x in extra:
        if x.startswith('deep') and x[4:]: DEEP = int(x[4:])
    rl = ('dup', 'diff', 'sum') + (('deep',) if any(x.startswith('deep') for x in extra) else ()) + (('vjoin',) if 'vjoin' in extra else ())
    ds = next((int(x[4:]) for x in extra if x.startswith('desc')), 0)
    dg = Counter() if 'diag' in extra else None
    rc = reclaim(B, W, use_v='nov' not in extra, lazy=lz, rels=rl, desc=ds, diag=dg)
    if dg is not None:
        for k_, v_ in sorted(dg.items(), key=lambda kv: -kv[1])[:30]: print('diag', v_, k_)
    Bn = assemble(B, W, rc)
    print('signals', check_signals(Bn), flush=True)
    if 'defer' in cfg: gfix(Bn)
    w = c7.walk(Bn, check='nocheck' not in extra)
    m = B['m']; H = w['H']
    D = w['W'] * m - w['sum']; lam = sum(n * r * math.log(m / r) for r, n in H.items())
    out = dict(h=h, cfg=cfg, R=Bn['R'], W=w['W'], s=w['sum'], sum_ok=w['sum_ok'], bad=w['bad'], maxrank=w['maxrank'], est=D / lam)
    print(out, flush=True)
    if 'cert' in extra:
        a, gap = cert(H, w['W'], m); print(dict(a_c=str(a), a_c_float=float(a), gap=float(gap)), flush=True)
