"""Independent exact Q(i) dirty-scratch replay of a complete reclaimed complex word WITH completed-core sharing.
Extends our earlier unshared replay (same dual-point model, same Lift identities, same endpoint) to shared banks; imports
nothing from the construction (reads build_word.py's pickle: ops with explicit updates + producer supports).

Sharing model (PR #128 Sec. 2, our merge variant): per stage, the 2024 cores are partitioned into groups J with binary
Gram I; the cores of a group run CONSECUTIVELY on ONE bank of R dirty Gaussian-rational streams (round j runs the j-th
core of every group, so each core really receives the physical state the previous core of its group left). A core
takes the incoming scratch at its first frames as a pure label (no child). After the group, each stream gets ONE
exterior child on Rfix = (sum_b R_b)^perp, R_b = lifted(last_b) minus lifted(start_b), computed LITERALLY on 576-bit
vectors: the R_b must be pairwise orthogonal (checked on every vector pair) and form a direct sum; Rfix must be
nondegenerate and proper; its phase is wt(P_Rfix eta) mod 4 from an exact Gram inverse (valid for alternating Rfix
too; alternating ones are counted, since they rely on the PR #24 normal form). The phase is NOT taken from the
reference; the reference ratio is compared and mismatches are counted, and restoration is checked on the values.
Mode 'merge': the exterior of a mergeable stream is fused into one existing step of one core (stage one: the step
into F at the stream's last touch in the group's LAST core; stage two: the stream's first step in the group's FIRST
core), applied at the time of that step, i.e. before that op's updates; one child of rank e + |Rfix|.
Mode 'single': every core is its own group (reproduces the per-core wrap; histogram must equal c7.walk's).
PASS: every pair (x, y) -> (F y, F x) (endpoint as complex-gate), every stream of every bank ends at C_full r0,
histogram == the construction's share.shared(..., mode) histogram.
Usage: python3 sreplay.py word.pkl SEED MODE [control]   MODE in merge, own, single"""
import sys, os, time, pickle, random, json, itertools
from collections import Counter, defaultdict
from fractions import Fraction as Q
import numpy as np
import gf2

COS = np.array([1, 0, -1, 0], np.int64); SIN = np.array([0, 1, 0, -1], np.int64)
LIM = 1 << 50
KNOWN = ('gate', 'latecopy', 'xcopy', 'vgate', 'inject', 'shear', 'move', 'high')


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
    """h-level frames with per-column lifted phases for one stage (as in our earlier unshared replay)."""
    def __init__(s, h, stage, eta, T, stats):
        s.h = h; s.stage = stage; s.st = stats
        E = [[(eta >> (i * h + j)) & 1 for j in range(h)] for i in range(h)]
        tv = [sum(1 << p for p in t) for t in T]
        if stage == 1: col = [sum((sum(E[i][j] for j in range(h) if t >> j & 1) & 1) << i for i in range(h)) for t in tv]
        else: col = [sum((sum(E[i][j] for i in range(h) if t >> i & 1) & 1) << j for j in range(h)) for t in tv]
        s.eta = np.array(col, np.uint32)
        s.ids = {}; s.basis = []; s.wtc = {}; s.mv = {}; s.mlist = []; s.mrev = []
    def fid(s, vs):
        k = gf2.key(vs)
        if k not in s.ids: s.ids[k] = len(s.basis); s.basis.append(list(k))
        return s.ids[k]
    def par(s, b): return (np.bitwise_count(np.uint32(b) & s.eta) & 1).astype(np.int64)
    def wt(s, f):
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
        O = gf2.orthonormal(B)
        if O is None: s.st['alt_distinct'] += 1; return None
        ph = np.zeros(len(s.eta), np.int64)
        for o in O: ph += s.par(o) * ((3 * (gf2.pc(o) % 4)) % 4)
        return ph % 4
    def move(s, u, w):
        """move id of the child u -> w (phase per column, rank, alternating?)."""
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
            if (ph != ref).any(): raise Fail('child phase != reference ratio')
        s.st['moves_distinct'] += 1
        s.mlist.append((ph.astype(np.int8), len(hi) - len(lo), alt)); s.mrev.append(k); s.mv[k] = len(s.mlist) - 1
        return s.mv[k]


# ---------------------------------------------------------------- partitions
def ag_partition(h, T):
    """triples of AG(k,2), h = 2^k: each parallel class of planes (4-sets) gives h orthonormal faces."""
    k = h.bit_length() - 1; assert 1 << k == h
    tid = {t: i for i, t in enumerate(T)}; seen = set(); out = []
    for a in range(1, h):
        for b in range(a + 1, h):
            D = frozenset([0, a, b, a ^ b])
            if D in seen: continue
            seen.add(D); cos = set()
            for x in range(h): cos.add(frozenset(x ^ y for y in D))
            grp = []
            for Qd in sorted(cos, key=sorted):
                for f in itertools.combinations(sorted(Qd), 3): grp.append(tid[f])
            out.append(grp)
    return out


def partition(h, T, sizes, spec):
    if spec == 'pr128':
        d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'certificates', 'round8', 'pr128_partition_mrp24.json')))
        lex = list(itertools.combinations(range(24), 3)); tid = {t: i for i, t in enumerate(T)}
        G = [[tid[lex[i]] for i in g] for g in d['groups']]
    elif spec == 'single': G = [[b] for b in range(len(T))]
    else:
        cls = ag_partition(h, T); want = sorted(sizes, reverse=True); G = []
        full = [s for s in want if s == h]; part = [s for s in want if s != h]
        ci = 0
        for _ in full: G.append(cls[ci]); ci += 1
        pools = [list(c) for c in cls[ci:]]                      # first-fit decreasing inside whole classes
        for sz in part:
            k = next(i for i, pl in enumerate(pools) if len(pl) >= sz)
            G.append(pools[k][:sz]); pools[k] = pools[k][sz:]
        assert not any(pools), 'partial sizes do not fill whole classes'
    return G


def check_groups(G, tvec, v):
    flat = sorted(b for g in G for b in g)
    if flat != list(range(v)): raise Fail('partition does not cover the cores exactly')
    bad = 0
    for g in G:
        for i, a in enumerate(g):
            if gf2.pc(tvec[a]) & 1 == 0: bad += 1
            for b in g[i + 1:]:
                if gf2.pc(tvec[a] & tvec[b]) & 1: bad += 1
    return bad


# ---------------------------------------------------------------- literal 576-bit spaces
class Space:
    def __init__(s, vecs):
        s.B = list(gf2.reduce_basis(vecs).values()); s.n = len(s.B); s._gi = None
    def nondeg(s): return gf2.gram_rank(s.B) == s.n
    def alternating(s): return s.n > 0 and not gf2.has_odd(s.B)
    def phase(s, eta):
        """wt(P_U eta) mod 4 (exact, Gram inverse; any nondegenerate U)."""
        if not s.B: return 0
        if s._gi is None: s._gi = gf2.gram_inverse(s.B)
        y = 0
        for j, b in enumerate(s.B):
            if gf2.dot(b, eta): y |= 1 << j
        p = 0
        for i, b in enumerate(s.B):
            if gf2.pc(s._gi[i] & y) & 1: p ^= b
        return gf2.pc(p) % 4


def cross_orth(A, B):
    return all(not gf2.dot(a, b) for a in A for b in B)


# ---------------------------------------------------------------- main
def run(path, seed, mode, ctl=None):
    t0 = time.time()
    d = pickle.load(open(path, 'rb'))
    h = d['h']; T = d['T']; v = len(T); m = h * h; N = v * v; R = d['R']
    tid = {t: i for i, t in enumerate(T)}; Fh = [1 << p for p in range(h)]
    tvec = [sum(1 << p for p in t) for t in T]
    rnd = random.Random(seed); eta = rnd.getrandbits(m)
    stats = Counter(); hist = Counter(); flags = Counter()
    retslot = d['retslot']
    names = sorted(retslot.items(), key=lambda kv: kv[1][1])
    if [nm[0] for _, nm in names] != ['E'] * h: raise Fail('retained basis is not E_0..E_h-1')
    ret_slots = [q for q, _ in names]; rset = set(ret_slots)
    # ---------- instruction list with explicit updates ----------
    den = h - 3; L42 = 2 * (h - 3)
    scat = {S: {('s', ret_slots[i]): Q(1, den) - (Q(1, 2) if i in S else 0) for i in range(h)} for S in T}
    scat42 = np.array([[int(scat[S][('s', ret_slots[i])] * L42) for i in range(h)] for S in T], np.int64)
    fwd = []
    for op in d['ops']:
        roles, F, kind, ups = op[:4]
        if kind == 'cread':
            got = defaultdict(dict)
            for dst, src, cf in ups: got[dst[1]][tuple(src)] = got[dst[1]].get(tuple(src), 0) + Q(cf)
            if any(got[S] != scat[S] for S in T): raise Fail('cread updates differ from the scatter identity')
            fwd.append([list(roles), F, kind, None]); continue
        if kind == 'read': fwd.append([list(roles), F, kind, None]); continue
        if kind not in KNOWN: raise Fail('unknown op ' + kind)
        u = []
        for dst, src, cf in ups or ():
            cf = Q(cf)
            if tuple(dst) == tuple(src): raise Fail('self update')
            if cf.denominator not in (1, 2): raise Fail('coefficient denominator')
            if tuple(dst)[0] == 'x': raise Fail('write into a source data role')
            u.append((tuple(dst), tuple(src), cf))
        rs = set(tuple(r) for r in roles)
        if any(a not in rs or b not in rs for a, b, _ in u): raise Fail('update role not among the op roles')
        fwd.append([list(roles), F, kind, u])
    stats['ops'] = len(fwd)
    # ---------- garbage map by the adjoint (x = 0), scaled by L42 ----------
    ih = next(i for i, w in enumerate(fwd) if w[2] == 'high')
    adj = np.zeros((R, v), np.int64); adjy = np.eye(v, dtype=np.int64) * L42
    for i in range(ih - 1, -1, -1):
        roles, F, kind, ups = fwd[i]
        if kind == 'read': continue
        if kind == 'cread':
            for j, q in enumerate(ret_slots): adj[q] += exdiv(scat42[:, j] @ adjy, L42)
            continue
        for dst, src, c in reversed(ups):
            if src[0] == 'x': continue
            a_d = adj[dst[1]] if dst[0] == 's' else adjy[tid[dst[1]]]
            inc = a_d * c.numerator
            if c.denominator != 1: inc = exdiv(inc, c.denominator)
            if src[0] == 's': adj[src[1]] += inc
            else: adjy[tid[src[1]]] += inc
    if int(np.abs(adj).max()) >= 2 ** 15: raise Fail('G entries too large for int16')
    G42 = adj.T.astype(np.int16); del adj
    stats['G_nnz'] = int((G42 != 0).sum())
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
        elif kind != 'cread' and ups:
            for dst, _, _ in ups:
                if dst[0] == 's': written.add(dst[1])
    if bad: raise Fail('readout checks %s' % dict(bad))
    if ctl == 'skipread':
        k = next(i for i, w in enumerate(fwd) if w[2] == 'read' and G42[:, w[3][0]].any()); fwd[k][2] = 'noop'
    # forward chains (for the merge rule: last touch is a step into F, and the first frame lies inside X_last)
    ch = defaultdict(list); lastop = {}
    for i, (roles, F, kind, ups) in enumerate(fwd):
        if kind == 'noop': continue
        for r in roles:
            if r[0] == 's' and kind != 'cread': ch[r[1]].append(F); lastop[r[1]] = i
    Fkey = gf2.key(Fh)
    mergeable = {}
    for q in range(R):
        c_ = ch[q]
        if gf2.key(c_[-1]) != Fkey: raise Fail('slot does not end at F')
        if len(c_) >= 2 and gf2.key(c_[-2]) != Fkey and gf2.contains(list(c_[-2]), list(c_[0])):
            mergeable[q] = h - gf2.dim(list(c_[-2]))
    stats['mergeable'] = len(mergeable)
    # ---------- partition ----------
    spec = 'single' if mode == 'single' else ('pr128' if h == 24 else 'ag')
    G = partition(h, T, d['sizes'], spec)
    if ctl == 'badgroup':                   # swap two triples between groups so one group gets an odd-overlap pair
        done = False
        for gi, ga in enumerate(G):
            for gj in range(gi + 1, len(G)):
                for x in range(len(G[gj])):
                    b = G[gj][x]
                    for y in range(len(ga)):
                        rest = [a for k, a in enumerate(ga) if k != y]
                        if any(gf2.pc(tvec[a] & tvec[b]) & 1 for a in rest):
                            ga[y], G[gj][x] = b, ga[y]; done = True; break
                    if done: break
                if done: break
            if done: break
    flags['group_gram_bad'] = check_groups(G, tvec, v)
    if mode != 'single' and sorted(len(g) for g in G) != sorted(d['sizes']): raise Fail('group sizes differ from the construction')
    stats['groups'] = len(G)
    print('instructions ready', round(time.time() - t0), 's', dict(stats), dict(flags), flush=True)
    # ---------- inputs ----------
    rng = np.random.default_rng(seed)
    S0 = (h - 3) ** 2 * 4
    gi = lambda *sh: (rng.integers(-1000, 1001, sh, dtype=np.int64) * S0, rng.integers(-1000, 1001, sh, dtype=np.int64) * S0)
    xin = gi(v, v); yin = gi(v, v)
    Fe = gf2.pc(eta) % 4
    cpl = {}
    def comp(F):
        if F not in cpl: cpl[F] = tuple(gf2.perp_within(Fh, list(F)))
        return cpl[F]
    Erow = [(eta >> (i * h)) & ((1 << h) - 1) for i in range(h)]
    Ecol = [sum(((eta >> (i * h + j)) & 1) << i for i in range(h)) for j in range(h)]
    Ot = [gf2.orthonormal(gf2.perp_within(Fh, [tv])) for tv in tvec]
    def phiD(stage):
        rows = Erow if stage == 1 else Ecol
        return np.array([sum(gf2.dot(r, o) * (gf2.pc(o) % 4) for r in rows for o in O) % 4 for O in Ot], np.int64)
    ctx = {}

    # ================= compile one stage: frames are identical in every core, so the move/update program is too
    def compile_stage(stage, ins, lf):
        lab = (lambda F: F) if stage == 1 else comp
        PHI = phiD(stage); ctx['PHI%d' % stage] = PHI
        fF = lf.fid(Fh)
        if ((PHI + lf.wt(fF)) % 4 != Fe).any(): raise Fail('orthogonal split of the full space')
        first = {}
        for roles, F, kind, ups in ins:
            if kind == 'noop': continue
            for r in roles:
                if r[0] == 's' and r[1] not in first and not (kind == 'cread' and r[1] in rset): first[r[1]] = F
        if len(first) != R: raise Fail('a slot is never touched')
        afr = np.array([lf.fid(lab(first[q])) for q in range(R)], np.int64); start = afr.copy()
        if stage == 1:
            dfr = {'x': np.array([lf.fid([tvec[t]]) for t in range(v)], np.int64), 'y': np.full(v, lf.fid([]), np.int64)}
            dpd = {'x': np.zeros(v, bool), 'y': np.zeros(v, bool)}
        else:
            dummy = lf.fid([])
            dfr = {'x': np.full(v, dummy, np.int64), 'y': np.full(v, dummy, np.int64)}
            dpd = {'x': np.ones(v, bool), 'y': np.ones(v, bool)}
        ci = next(i for i, w in enumerate(ins) if w[2] == 'cread')
        nxt = {}
        for w in ins[ci + 1:]:
            for r in w[0]:
                if r[0] == 's' and r[1] in rset and r[1] not in nxt: nxt[r[1]] = w[1]
        prog = []
        # merge sites: stage one = the move at the slot's last touch (a step into F); stage two = the slot's first move
        site = {}
        if stage == 1:
            for q, e in mergeable.items(): site[q] = ('last', lastop[q])
        firstmove_seen = set()
        def amove(q, f, oi=None):
            u = int(afr[q])
            if u == f: return
            mid = lf.move(u, f); tag = None
            if q in mergeable:
                if stage == 1 and oi == site[q][1]: tag = q
                if stage == 2 and q not in firstmove_seen:
                    tag = q
                    if lf.mlist[mid][1] != mergeable[q]: raise Fail('stage-two first step rank != forward last step rank')
            firstmove_seen.add(q)
            prog.append(('A', q, mid, tag)); afr[q] = f
        def entr(kind, rows, f):
            want = (lambda a2: [tvec[a2]]) if kind == 'y' else (lambda a2: [])
            for a2 in rows:
                if gf2.key(lf.basis[f]) != gf2.key(want(a2)): raise Fail('stage-two entrance frame')
            ph = ctx['PSI'][:, rows].T
            ref = (PHI[None, :] + lf.wt(f)[None, :] - ctx['OLDWT'][kind][:, rows].T) % 4
            if (ph % 4 != ref).any(): raise Fail('entrance phase != reference')
            prog.append(('E', kind, np.asarray(rows, np.int64)))
            dfr[kind][rows] = f; dpd[kind][rows] = False
        def dmove(kind, idx, f):
            idx = np.asarray(idx, np.int64)
            if dpd[kind][idx].any():
                e = idx[dpd[kind][idx]]; entr(kind, e, f); idx = idx[~dpd[kind][idx]]
            if not len(idx): return
            cur = dfr[kind][idx]
            for u in np.unique(cur):
                if u == f: continue
                sel = idx[cur == u]; mid = lf.move(int(u), f)
                prog.append(('D', kind, sel, mid)); dfr[kind][sel] = f
        for oi, (roles, F, kind, ups) in enumerate(ins):
            if kind == 'noop': continue
            f = lf.fid(lab(F))
            if kind == 'cread':
                cps = []
                for j, q in enumerate(ret_slots):
                    if stage == 2: amove(q, lf.fid(lab(nxt[q])))
                    u = int(afr[q]); cps.append((j, q, -1 if u == f else lf.move(u, f)))
                ys = [tid[r[1]] for r in roles if r[0] == 'y']
                dmove('y', ys, f)
                prog.append(('C', cps, np.array(ys, np.int64), 1 if stage == 1 else -1)); continue
            if kind == 'read':
                q, tg = ups
                amove(q, f, oi); dmove('y', tg, f)
                prog.append(('R', q, tg, -1 if stage == 1 else 1))
                nx = ins[oi + 1] if oi + 1 < len(ins) else None
                if nx is None or nx[2] != 'read' or lf.fid(lab(nx[1])) != f: prog.append(('FL',))
                continue
            n0 = len(prog)
            for r in roles:
                if r[0] == 's': amove(r[1], f, oi)
                else: dmove(r[0], [tid[r[1]]], f)
            tagged = [x[3] for x in prog[n0:] if x[0] == 'A' and x[3] is not None]
            for dst, src, c in ups:
                ds = ('s', dst[1]) if dst[0] == 's' else (dst[0], tid[dst[1]])
                ss = ('s', src[1]) if src[0] == 's' else (src[0], tid[src[1]])
                prog.append(('U', ds, ss, c.numerator, c.denominator))
            for q in tagged: prog.append(('X', q))
        if stage == 1:
            if any(gf2.key(lf.basis[int(f)]) != Fkey for f in dfr['x']): raise Fail('stage-one X end frame')
            if any(gf2.key(lf.basis[int(f)]) != gf2.key(gf2.perp_within(Fh, [tvec[t]])) for t, f in enumerate(dfr['y'])):
                raise Fail('stage-one Y end frame')
        else:
            if dpd['x'].any() or dpd['y'].any(): raise Fail('data role untouched in stage two')
        ctx['dfr%d' % stage] = dfr
        return prog, afr.copy(), start, lab

    # ================= execute one round (one core per active group) on a bank
    def exec_stage(prog, lf, cols, bank, data, mtab, mctx):
        """bank: (re, im) R x nc; data: kind -> (re, im) v x nc; mctx: q -> (extra phase (nc,), [(rank_add, count)], n_plain)."""
        nc = len(cols); are, aim = bank
        PSI_T = ctx['PSI'][cols] if 'PSI' in ctx and lf.stage == 2 else None
        pend = []
        def arr(spec):
            if spec[0] == 's': return are, aim, spec[1]
            D = data[spec[0]]; return D[0], D[1], spec[1]
        for ins_ in prog:
            t = ins_[0]
            if t == 'A':
                _, q, mid, tag = ins_; ph = mtab[mid]; rk = lf.mlist[mid][1]
                if tag is not None and tag in mctx:
                    extra, adds, nplain = mctx[tag]
                    if ctl != 'mergeafter': ph = (ph.astype(np.int64) + extra) % 4
                    for ra, cnt in adds: hist[rk + ra] += cnt
                    if nplain: hist[rk] += nplain
                else: hist[rk] += nc
                are[q], aim[q] = rot(are[q], aim[q], ph)
                if lf.mlist[mid][2]: stats['alt_children'] += nc
            elif t == 'X':
                if ctl == 'mergeafter' and ins_[1] in mctx:
                    q = ins_[1]; are[q], aim[q] = rot(are[q], aim[q], mctx[q][0])
            elif t == 'U':
                _, ds, ss, num, dn = ins_
                dre, dim_, di = arr(ds); sre, sim, si = arr(ss)
                for dp, sp in ((dre, sre), (dim_, sim)):
                    inc = sp[si] * num
                    if dn != 1: inc = exdiv(inc, dn)
                    dp[di] += inc
                if int(max(np.abs(dre[di]).max(), np.abs(dim_[di]).max())) >= LIM: raise Fail('magnitude bound exceeded')
            elif t == 'D':
                _, kind, sel, mid = ins_; ph = mtab[mid]; D = data[kind]
                D[0][sel], D[1][sel] = rot(D[0][sel], D[1][sel], ph[None, :])
                hist[lf.mlist[mid][1]] += nc * len(sel)
                if lf.mlist[mid][2]: stats['alt_children'] += nc * len(sel)
            elif t == 'E':
                _, kind, rows = ins_; D = data[kind]
                ph = PSI_T[:, rows].T
                D[0][rows], D[1][rows] = rot(D[0][rows], D[1][rows], ph)
                hist[(h - 1) ** 2] += nc * len(rows)
            elif t == 'R':
                pend.append(ins_[1:])
            elif t == 'FL':
                qs = np.array([q for q, _, _ in pend]); sg = pend[0][2]
                rows = np.unique(np.concatenate([tg for _, tg, _ in pend]))
                Gs = G42[np.ix_(rows, qs)].astype(np.float64)
                Yre, Yim = data['y'][0], data['y'][1]
                for part, src in ((Yre, are), (Yim, aim)):
                    Z = src[qs]; chk(Z)
                    if float(np.abs(Gs).sum(axis=1).max()) * float(np.abs(Z).max() + 1) >= 2.0 ** 52: raise Fail('float matmul bound')
                    part[rows] += sg * exdiv(np.rint(Gs @ Z.astype(np.float64)).astype(np.int64), L42)
                pend.clear()
            elif t == 'C':
                _, cps, ys, sg = ins_
                cre = np.zeros((h, nc), np.int64); cim = np.zeros((h, nc), np.int64)
                for j, q, mid in cps:
                    if mid < 0: cre[j], cim[j] = are[q], aim[q]
                    else:
                        cre[j], cim[j] = rot(are[q], aim[q], mtab[mid]); hist[lf.mlist[mid][1]] += nc
                        if lf.mlist[mid][2]: stats['alt_children'] += nc
                Yre, Yim = data['y'][0], data['y'][1]
                for part, cp in ((Yre, cre), (Yim, cim)):
                    chk(cp); part[ys] += sg * exdiv(scat42[ys] @ cp, L42)
            else: raise Fail('bad program op')
        if pend: raise Fail('unflushed reads')

    # ================= literal exterior spaces per (group, slot class)
    def exteriors(stage, lf, start, last, lab):
        """per group: {class: (Rfix Space, |Rfix|)}; checks pairwise orthogonality and directness literally."""
        lift = (lambda t, r: gf2.kron(r, tvec[t], h)) if stage == 1 else (lambda t, r: gf2.kron(tvec[t], r, h))
        classes = defaultdict(list)
        for q in range(R): classes[(int(start[q]), int(last[q]))].append(q)
        res_h = {}
        for (fs, fl_) in classes:
            Ls, Ll = lf.basis[fs], lf.basis[fl_]
            if not gf2.contains(Ll, Ls): raise Fail('core start frame not inside its last frame')
            r = gf2.perp_within(Ll, Ls)
            if len(r) != len(Ll) - len(Ls): raise Fail('core residual dimension')
            res_h[(fs, fl_)] = r
        ext = []; units = [1 << p for p in range(m)]
        cache = {}
        for gi_, g in enumerate(G):
            ent = {}
            for cl, r in res_h.items():
                key = (tuple(sorted(g)) if ctl != 'wrapfix' else tuple(g), cl)
                if key in cache: ent[cl] = cache[key]; continue
                parts = [[lift(t, x) for x in r] for t in g]
                if ctl == 'wrapfix' and len(parts) > 1: parts = parts[1:]          # control: forget one core's residual
                orth_bad = 0
                for i in range(len(parts)):
                    for j in range(i + 1, len(parts)):
                        if not cross_orth(parts[i], parts[j]): orth_bad += 1
                allv = [x for p in parts for x in p]
                Gs = gf2.reduce_basis(allv)
                direct = len(Gs) == len(allv)
                fix = Space(gf2.perp_within(units, list(Gs.values())))
                if fix.n != m - len(Gs): raise Fail('exterior dimension')
                nd = fix.nondeg()
                if not nd and ctl is None: raise Fail('exterior degenerate')
                ent[cl] = dict(fix=fix, n=fix.n, orth_bad=orth_bad, direct=direct, nondeg=nd, alt=fix.alternating())
                cache[key] = ent[cl]
            ext.append(ent)
        return ext, classes, res_h

    # ================= one stage over all groups
    def stage_all(stage, ins, lf, data_full):
        prog, last, start, lab = compile_stage(stage, ins, lf)
        mtab_full = np.stack([x[0] for x in lf.mlist]) if lf.mlist else np.zeros((0, v), np.int8)
        print('stage %d compiled' % stage, len(prog), 'instructions', len(lf.mlist), 'moves', round(time.time() - t0), 's', flush=True)
        ext, classes, res_h = exteriors(stage, lf, start, last, lab)
        cls_of = {q: cl for cl, qs in classes.items() for q in qs}
        for g_ext in ext:
            for e in g_ext.values():
                flags['ext_orth_bad'] += e['orth_bad']; flags['ext_not_direct'] += (not e['direct'])
                flags['ext_degenerate'] += (not e['nondeg']); flags['ext_alternating_nonzero'] += (e['alt'] and e['n'] > 0)
        # literal fix-up phases per group/class at this eta, compared with the reference ratio
        fixph = []
        for gi_, g in enumerate(G):
            ent = {}
            for cl, e in ext[gi_].items():
                ph = e['fix'].phase(eta) if e['nondeg'] else 0
                ref = (Fe - sum(int(lf.wt(cl[1])[b]) - int(lf.wt(cl[0])[b]) for b in g)) % 4
                if e['nondeg'] and ph != ref: flags['fix_phase_ne_reference'] += 1
                ent[cl] = ph
            fixph.append(ent)
        # merge legality on the literal spaces: the step residual of the merge core is orthogonal to Rfix
        merged = set()
        if mode == 'merge' and ctl != 'nomerge':
            lift = (lambda t, r: gf2.kron(r, tvec[t], h)) if stage == 1 else (lambda t, r: gf2.kron(tvec[t], r, h))
            sitecore = (lambda g: g[-1]) if stage == 1 else (lambda g: g[0])
            if ctl == 'm1first' and stage == 1: sitecore = lambda g: g[0]
            if ctl == 'm2last' and stage == 2: sitecore = lambda g: g[-1]
            stepres = {}
            for ins_ in prog:
                if ins_[0] == 'A' and ins_[3] is not None:
                    q = ins_[1]; u, w = None, None
                    stepres[q] = ins_[2]
            mres = {}
            for gi_, g in enumerate(G):
                for q, mid in stepres.items():
                    e = ext[gi_][cls_of[q]]
                    if e['n'] == 0: continue
                    key = (gi_, mid)
                    if key not in mres:
                        uv = lf.mrev[mid]
                        U, W = lf.basis[uv[0]], lf.basis[uv[1]]
                        lo, hi = (U, W) if gf2.contains(W, U) else (W, U)
                        rh = gf2.perp_within(hi, lo); rl = [lift(sitecore(g), x) for x in rh]
                        ok = cross_orth(rl, e['fix'].B) and len(gf2.reduce_basis(rl + e['fix'].B)) == len(rl) + e['n']
                        if len(rl) + e['n'] > m - 1: ok = False
                        mres[key] = ok
                        if not ok: flags['merge_not_orthogonal_or_improper'] += 1
                    merged.add((gi_, q))
            flags['merged_children_stage%d' % stage] = len(merged)
        # banks
        ng = len(G); bre, bim = gi(R, ng); r0 = (bre.copy(), bim.copy())
        L = max(len(g) for g in G)
        for j in range(L):
            gidx = [k for k in range(ng) if len(G[k]) > j]
            if ctl == 'm1first' and stage == 1: site_round = lambda k: j == 0
            elif ctl == 'm2last' and stage == 2: site_round = lambda k: j == len(G[k]) - 1
            else: site_round = (lambda k: j == len(G[k]) - 1) if stage == 1 else (lambda k: j == 0)
            cols = np.array([G[k][j] for k in gidx], np.int64)
            mtab = mtab_full[:, cols]
            bank = (bre[:, gidx], bim[:, gidx])
            if ctl == 'swapbank' and j == 1: bank[0][[0, 1]] = bank[0][[1, 0]]; bank[1][[0, 1]] = bank[1][[1, 0]]
            mctx = {}
            sg_ = [ci_ for ci_, k in enumerate(gidx) if site_round(k)] if merged else []
            if sg_:
                for q in mergeable:
                    extra = np.zeros(len(gidx), np.int64); adds = Counter(); hit = 0
                    for ci_ in sg_:
                        k = gidx[ci_]
                        if (k, q) in merged:
                            extra[ci_] = fixph[k][cls_of[q]]; adds[ext[k][cls_of[q]]['n']] += 1; hit += 1
                    if hit: mctx[q] = (extra, sorted(adds.items()), len(gidx) - hit)
            sub = {kd: (data_full[kd][0][:, cols].copy(), data_full[kd][1][:, cols].copy()) for kd in ('x', 'y')}
            exec_stage(prog, lf, cols, bank, sub, mtab, mctx)
            for kd in ('x', 'y'): data_full[kd][0][:, cols] = sub[kd][0]; data_full[kd][1][:, cols] = sub[kd][1]
            bre[:, gidx] = bank[0]; bim[:, gidx] = bank[1]
        # exteriors that were not merged
        e_re, e_im = rot(r0[0], r0[1], Fe)
        for k in range(ng):
            for cl, qs in classes.items():
                e = ext[k][cl]
                if e['n'] == 0 or ctl == 'nofix': continue
                qs_ = [q for q in qs if (k, q) not in merged]
                if not qs_: continue
                bre[qs_, k], bim[qs_, k] = rot(bre[qs_, k], bim[qs_, k], fixph[k][cl])
                hist[e['n']] += len(qs_)
                if e['alt']: stats['alt_exteriors'] += len(qs_)
        notrest = int(((bre != e_re) | (bim != e_im)).sum())
        stats['aux_not_restored_s%d' % stage] = notrest
        print('stage %d done' % stage, round(time.time() - t0), 's', dict(stats), flush=True)

    LF1 = Lift(h, 1, eta, T, stats); LF2 = Lift(h, 2, eta, T, stats)
    data1 = {'x': (xin[0].copy(), xin[1].copy()), 'y': (yin[0].copy(), yin[1].copy())}
    stage_all(1, fwd, LF1, data1)
    # entrance phases (as in our earlier unshared replay)
    Em = np.array([[(eta >> (i * h + j)) & 1 for j in range(h)] for i in range(h)], np.float64)
    Om = [np.array([[(o >> p) & 1 for p in range(h)] for o in O], np.float64) for O in Ot]
    Wv = [np.array([gf2.pc(o) % 4 for o in O], np.int64) for O in Ot]
    Qall = np.concatenate([Em @ Om[a2].T for a2 in range(v)], axis=1)
    W2 = np.stack(Wv)
    PSI = np.zeros((v, v), np.int64)
    for a1 in range(v):
        M = (np.rint(Om[a1] @ Qall).astype(np.int64) % 2).reshape(h - 1, v, h - 1)
        PSI[a1] = (Wv[a1] @ (M * W2[None, :, :]).sum(axis=2)) % 4
    wtF1 = LF1.wt(LF1.fid(Fh))
    ctx['PSI'] = PSI
    ctx['OLDWT'] = {'y': np.tile(wtF1, (v, 1)), 'x': np.stack([LF1.wt(LF1.fid(gf2.perp_within(Fh, [tvec[a1]]))) for a1 in range(v)])}
    data2 = {'y': (data1['x'][0].T.copy(), data1['x'][1].T.copy()), 'x': (data1['y'][0].T.copy(), data1['y'][1].T.copy())}
    del data1
    inv = [[w[0], w[1], w[2], ([(a, b, -c) for a, b, c in reversed(w[3])] if w[2] not in ('read', 'cread', 'noop') else w[3])]
           for w in reversed(fwd)]
    stage_all(2, inv, LF2, data2)
    Xre, Xim = data2['y']; Yre, Yim = data2['x']
    dfr = ctx['dfr2']; fF2 = LF2.fid(Fh)
    for a2 in range(v):
        for (re_, im_, fr_), tgt in (((Xre, Xim, dfr['y']), fF2), ((Yre, Yim, dfr['x']), LF2.fid(gf2.perp_within(Fh, [tvec[a2]])))):
            u = int(fr_[a2])
            if u != tgt:
                mid = LF2.move(u, tgt); ph, rk, alt = LF2.mlist[mid]
                re_[a2], im_[a2] = rot(re_[a2], im_[a2], ph); hist[rk] += v; fr_[a2] = tgt
    hist[1] += N
    A = (Xre.T, Xim.T); Bv = (Yre.T, Yim.T)
    Tm = np.array([[(tv >> p) & 1 for p in range(h)] for tv in tvec], np.float64)
    cw = (np.rint(Tm @ Em @ Tm.T).astype(np.int64)) % 2
    Pe = (9 * cw) % 4
    res = {}
    eq = lambda a, b: int(((a[0] != b[0]) | (a[1] != b[1])).sum())
    add = lambda a, b: (a[0] + b[0], a[1] + b[1])
    res['A'] = eq(A, rot(*yin, Fe + 2))
    res['B'] = eq(Bv, add(rot(*xin, Fe - 2 * Pe), rot(*yin, Fe - Pe)))
    Bc = add(Bv, rot(*A, -cw))
    pw = (Fe + 9 - 2 * cw) % 4
    for _ in range(500):
        a1, a2 = rnd.randrange(v), rnd.randrange(v)
        if gf2.pc(eta ^ gf2.kron(tvec[a1], tvec[a2], h)) % 4 != pw[a1, a2]: raise Fail('popc(eta+w) formula')
    res['Yout'] = eq(rot(*Bc, 9), rot(*xin, pw))
    res['Xout'] = eq(rot(*A, 2), rot(*yin, Fe))
    H = {k: n for k, n in sorted(hist.items()) if n}
    ng = len(G); Wn = 2 * N + 2 * ng * R
    s = sum(k * n for k, n in H.items())
    want = d['walkH'] if mode == 'single' else d['shared'][mode].get('H')
    ok = not any(res.values()) and stats['aux_not_restored_s1'] == 0 and stats['aux_not_restored_s2'] == 0
    clean = ok and flags['ext_orth_bad'] == 0 and flags['ext_not_direct'] == 0 and flags['ext_degenerate'] == 0 and \
        flags['fix_phase_ne_reference'] == 0 and flags['merge_not_orthogonal_or_improper'] == 0 and flags['group_gram_bad'] == 0
    out = dict(h=h, seed=seed, mode=mode, control=ctl, PASS=clean and H == want, wrong=res,
               aux_not_restored=(stats['aux_not_restored_s1'], stats['aux_not_restored_s2']), hist_equals_construction=(H == want),
               W=Wn, s=s, D=Wn * m - s, maxrank=max(H), flags=dict(flags), moves_distinct=stats['moves_distinct'],
               alt_distinct=stats['alt_distinct'], alt_children=stats['alt_children'], alt_exteriors=stats['alt_exteriors'],
               mergeable=stats['mergeable'], groups=ng, G_nnz=stats['G_nnz'], secs=round(time.time() - t0))
    if want is not None and H != want:
        dk = sorted(set(H) | set(want)); out['hist_diff'] = {k: (H.get(k, 0), want.get(k, 0)) for k in dk if H.get(k, 0) != want.get(k, 0)}
    tag = os.path.basename(path).replace('.pkl', '')
    pickle.dump(dict(H=H, W=Wn, R=R, groups=ng, m=m), open('hist_%s_%d_%s%s.pkl' % (tag, seed, mode, '_' + ctl if ctl else ''), 'wb'))
    return out


if __name__ == '__main__':
    path = sys.argv[1]; seed = int(sys.argv[2]); mode = sys.argv[3]; ctl = sys.argv[4] if len(sys.argv) > 4 else None
    try: print(run(path, seed, mode, ctl), flush=True)
    except Fail as e: print(dict(path=path, seed=seed, mode=mode, control=ctl, PASS=False, FAIL=str(e)), flush=True)
