"""Items 1-2 (algebra/F2 part), own code. Paired-cube motif of #144 at h = 2p: ports S = one coordinate from each of three
distinct pairs (coordinate 2i+b = pair i, selector b), q_S = chi_S. Exact checks over Q and F2:
  I = K + H + B with B_{T,S} = (|T cap S|-1)/2, K = I - B_block, H = B_block - B;
  K = (P - A)/2 inside a cube, K^2 = I, K maps each parity class to the other, its 4x4 blocks are Hadamard/2;
  H rebuilt from the signed channels F, A, E, G and the 35-call root-group schedule equals B_block - B;
  B equals the star scatter 1/2 sum_{i in T} S_i - 1/6 sum_i S_i;
  every nonzero H and K entry joins orthogonal addresses; every channel aggregate's support is orthogonal to its whole
  target group; group caps nest with dims h-4, h-3, h-2, h-1 (target chain [h-4,1,1,1]);
  source itinerary <q_S> < U_par < q_T^perp < F_A with Lagrangian distances 2, h-4, 1 (formula of the note and brute force
  on F2^{2h} for small h), U_par Gram rank;
  star spans U_c: dim h-2, equal to {u: u_cbar = 0, outside part even}, radical dim; l = sum dim U_c = h(h-2).
Usage: python3 cube_alg.py p [p ...]"""
import sys, itertools
from fractions import Fraction as Q
from collections import Counter
import numpy as np

def rank2(vs):
    b = {}
    for x in vs:
        for p_, y in b.items():
            if x >> p_ & 1: x ^= y
        if x:
            p_ = x.bit_length() - 1
            for k in list(b):
                if b[k] >> p_ & 1: b[k] ^= x
            b[p_] = x
    return b
def basis(vs): return list(rank2(vs).values())
def dim(vs): return len(rank2(vs))
def dot(a, b): return bin(a & b).count('1') & 1
def perp(vs, h):
    B = basis(vs)
    return basis([x for x in range(1 << h) if all(dot(x, b) == 0 for b in B)]) if h <= 14 else perp_big(B, h)
def perp_big(B, h):
    # solve B x = 0 over F2 (rows = basis vectors)
    rows = [(b, 0) for b in B]; piv = {}
    for b in B:
        x = b
        for c, y in piv.items():
            if x >> c & 1: x ^= y
        if x:
            c = x.bit_length() - 1
            for k in list(piv):
                if piv[k] >> c & 1: piv[k] ^= x
            piv[c] = x
    free = [j for j in range(h) if j not in piv]; out = []
    for f in free:
        x = 1 << f
        for c, y in piv.items():
            if y >> f & 1: x |= 1 << c
        out.append(x)
    return out
def inside(A, B):
    bb = rank2(B)
    for x in A:
        for p_, y in bb.items():
            if x >> p_ & 1: x ^= y
        if x: return False
    return True
def meet_dim(U, V): return dim(U) + dim(V) - dim(list(U) + list(V))
def dL(U, V):                       # note's formula d(L_U, L_V) = dim U + dim V - 2 dim(U cap V)
    return dim(U) + dim(V) - 2 * meet_dim(U, V)
def L_of(U, h):                     # L_U = {(u, u+v)} as 2h-bit ints (x in low h bits, z in high h bits)
    Up = perp(U, h); gens = [u | (u << h) for u in basis(U)] + [v << h for v in Up]
    return gens
def symp(a, b, h):
    m = (1 << h) - 1
    return (dot(a & m, b >> h) + dot(a >> h, b & m)) & 1
def dL_brute(U, V, h):
    LU, LV = L_of(U, h), L_of(V, h)
    assert dim(LU) == h and all(symp(a, b, h) == 0 for a in LU for b in LU), 'L_U not Lagrangian'
    return h - (dim(LU) + dim(LV) - dim(LU + LV))

def run(p):
    h = 2 * p; ok = True
    def rep(name, c):
        nonlocal ok
        ok &= bool(c); print('  %-100s %s' % (name, 'OK' if c else 'FAIL'), flush=True)
    cubes = list(itertools.combinations(range(p), 3))
    ports = [(I, bits) for I in cubes for bits in itertools.product(range(2), repeat=3)]
    idx = {pt: k for k, pt in enumerate(ports)}; v = len(ports)
    coords = lambda pt: frozenset(2 * i + b for i, b in zip(*pt))
    q = [sum(1 << c for c in coords(pt)) for pt in ports]
    print('p=%d h=%d v=%d cubes=%d' % (p, h, v, len(cubes)))
    Bm = np.array([[Q(len(coords(T) & coords(S)) - 1, 2) for S in ports] for T in ports], dtype=object)
    blk = np.zeros((v, v), dtype=object); blk[:] = Q(0)
    for T in range(v):
        for S in range(v):
            if ports[T][0] == ports[S][0]: blk[T, S] = Bm[T, S]
    I_ = np.zeros((v, v), dtype=object); I_[:] = Q(0)
    for t in range(v): I_[t, t] = Q(1)
    K = I_ - blk; Hm = blk - Bm
    rep('K + H + B == I (definitions)', np.all(K + Hm + Bm == I_))
    # K inside a cube
    Kc = K[:8, :8]; Pm = np.zeros((8, 8), dtype=object); Am = np.zeros((8, 8), dtype=object); Pm[:] = Q(0); Am[:] = Q(0)
    bl = list(itertools.product(range(2), repeat=3))
    for a in range(8):
        for b in range(8):
            dd = sum(x != y for x, y in zip(bl[a], bl[b]))
            if dd == 3: Pm[a, b] = Q(1)
            if dd == 1: Am[a, b] = Q(1)
    rep('K|cube == (P - A)/2', np.all(Kc == (Pm - Am) / 2))
    rep('K^2 == I (exact, K is cube-block-diagonal: checked per block)', all(np.all(K[8*c:8*c+8, 8*c:8*c+8].dot(K[8*c:8*c+8, 8*c:8*c+8]) == I_[:8, :8]) for c in range(len(cubes))) and all(K[T, S] == 0 for T in range(v) for S in range(v) if ports[T][0] != ports[S][0]))
    par = [sum(b) % 2 for _, b in ports]
    rep('K maps each parity class to the other (no same-parity entry)', all(K[T, S] == 0 for T in range(v) for S in range(v) if par[T] == par[S]))
    blocks_ok = True
    for c in range(len(cubes)):
        for pa in (0, 1):
            rows = [8 * c + k for k in range(8) if par[8 * c + k] != pa]; cols = [8 * c + k for k in range(8) if par[8 * c + k] == pa]
            M = np.array([[K[r, s] for s in cols] for r in rows], dtype=object)
            blocks_ok &= all(abs(x) == Q(1, 2) for x in M.flat) and np.all((2 * M).dot((2 * M).T) == 4 * np.eye(4, dtype=object))
    rep('every opposite-parity 4x4 block of K is (signed Hadamard)/2', blocks_ok)
    rep('every nonzero H or K entry joins orthogonal addresses', all(dot(q[T], q[S]) == 0 for T in range(v) for S in range(v)
                                                                     if (Hm[T, S] != 0 or K[T, S] != 0)))
    rep('distinct non-orthogonal addresses have B = 0', all(Bm[T, S] == 0 for T in range(v) for S in range(v) if T != S and dot(q[T], q[S])))
    # ---- channels and the 35-call schedule ----
    def src(I, bits): return idx[(I, tuple(bits))]
    F = {I: Counter({src(I, b): 1 for b in bl}) for I in cubes}
    A = {(I, I[k], a): Counter({src(I, b): 1 for b in bl if b[k] == a}) for I in cubes for k in range(3) for a in (0, 1)}
    E = {}
    for I in cubes:
        for k, l in itertools.combinations(range(3), 2):
            for a, b in itertools.product(range(2), repeat=2):
                E[I, I[k], I[l], a, b] = Counter({src(I, bb): 1 for bb in bl if bb[k] == a and bb[l] == b})
    def G(I, i, j, mode):
        if mode == 0: return E[I, i, j, 0, 0] - E[I, i, j, 1, 1] if False else Counter({**E[I, i, j, 0, 0], **{s: -1 for s in E[I, i, j, 1, 1]}})
        return Counter({**E[I, i, j, 0, 1], **{s: -1 for s in E[I, i, j, 1, 0]}})
    Hc = np.zeros((v, v), dtype=object); Hc[:] = Q(0)
    calls = Counter(); grp_ok = True; chain_ok = True; tchain = Counter()
    for J in cubes:
        groups = []          # (targets, aggregate counter, coefficient)
        agg = Counter()
        for I in cubes:
            if not set(I) & set(J): agg.update(F[I])
        groups.append(('F', [src(J, b) for b in bl], agg, Q(1, 2)))
        def face(k, key_bits):
            out = []
            for fix in itertools.product(range(2), repeat=len(key_bits)):
                tg = [src(J, b) for b in bl if all(b[kk] == f for kk, f in zip(key_bits, fix))]
                eta_k = fix[key_bits.index(k)]
                ag = Counter()
                for I in cubes:
                    if set(I) & set(J) == {J[k]}: ag.update(A[I, J[k], 1 - eta_k])
                out.append(('face%d' % k, tg, ag, Q(1, 2)))
            return out
        def edge(k, l, key_bits):
            out = []
            for fix in itertools.product(range(2), repeat=len(key_bits)):
                tg = [src(J, b) for b in bl if all(b[kk] == f for kk, f in zip(key_bits, fix))]
                ek, el = fix[key_bits.index(k)], fix[key_bits.index(l)]
                ag = Counter()
                for I in cubes:
                    if set(I) & set(J) == {J[k], J[l]}:
                        for s, c in G(I, J[k], J[l], ek ^ el).items(): ag[s] += c
                out.append(('edge%d%d' % (k, l), tg, ag, Q(2 * ek - 1, 2)))
            return out
        groups += face(0, (0,)) + face(1, (0, 1)) + edge(0, 1, (0, 1)) + face(2, (0, 1, 2)) + edge(0, 2, (0, 1, 2)) + edge(1, 2, (0, 1, 2))
        frames = {t: [] for t in [src(J, b) for b in bl]}
        for name, tg, ag, coef in groups:
            calls[name] += 1
            cap = perp([q[t] for t in tg], h)
            sup = [q[s] for s, c in ag.items() if c]
            grp_ok &= inside(sup, cap)                       # aggregate support orthogonal to the whole group
            for t in tg:
                for s, c in ag.items(): Hc[t, s] += coef * c
                frames[t].append(cap)
        order = {'F': 0, 'face0': 1, 'face1': 2, 'edge01': 2, 'face2': 3, 'edge02': 3, 'edge12': 3}
        for t, fr in frames.items():
            # chain in schedule order, ending at q_T^perp
            seq = [Uc for Uc in fr]
            ds = [dim(Uc) for Uc in seq]
            chain_ok &= all(inside(a, b) for a, b in zip(seq, seq[1:])) and ds[0] == h - 4 and ds[-1] == h - 1
            dd = sorted(set([0] + ds + [h - 1])); tchain[tuple(b - a for a, b in zip(dd, dd[1:]))] += 1
    rep('channels + root groups reproduce H = B_block - B exactly', np.all(Hc == Hm))
    rep('calls per target cube %s (total 35)' % dict((k, c // len(cubes)) for k, c in calls.items()), sum(calls.values()) == 35 * len(cubes))
    rep('every aggregate support lies in its group cap (orthogonal to the whole group)', grp_ok)
    rep('target chains nested, start at h-4, end at q_T^perp; step multiset %s' % dict(tchain), chain_ok and set(tchain) == {(h - 4, 1, 1, 1)})
    # ---- stars ----
    Sv = [Counter({s: 1 for s in range(v) if c in coords(ports[s])}) for c in range(h)]
    Bs = np.zeros((v, v), dtype=object); Bs[:] = Q(0)
    for t in range(v):
        for c in range(h):
            co = (Q(1, 2) if c in coords(ports[t]) else 0) - Q(1, 6)
            for s in Sv[c]: Bs[t, s] += co
    rep('star scatter 1/2 sum_{i in T} S_i - 1/6 sum_i S_i == B', np.all(Bs == Bm))
    dims = []; eqset = True; rad = Counter()
    for c in range(h):
        U = basis([q[s] for s in Sv[c]]); dims.append(len(U))
        cb = c ^ 1; outside = [j for j in range(h) if j not in (c, cb)]
        if h <= 14:
            claim = [u for u in range(1 << h) if not (u >> cb & 1) and sum(u >> j & 1 for j in outside) % 2 == 0]
            eqset &= len(claim) == 1 << len(U) and inside(claim, U)
        Ub = basis(U); G_ = [[dot(a, b) for b in Ub] for a in Ub]
        rad[len(Ub) - dim([sum(r[j] << j for j in range(len(Ub))) for r in G_])] += 1
    rep('star spans: dims %s (want h-2 = %d), l = %d (want h(h-2) = %d)' % (set(dims), h - 2, sum(dims), h * (h - 2)),
        set(dims) == {h - 2} and sum(dims) == h * (h - 2))
    if h <= 14: rep('star span == {u : u_cbar = 0, outside coordinates even}', eqset)
    print('    star span radical dims (degenerate frames, need the arbitrary-subspace interface): %s' % dict(rad))
    # ---- source itinerary for K ----
    it_ok = True; dl = Counter(); gram = Counter()
    for c in range(len(cubes)):
        for pa in (0, 1):
            cls = [8 * c + k for k in range(8) if par[8 * c + k] == pa]; opp = [8 * c + k for k in range(8) if par[8 * c + k] != pa]
            U = basis([q[s] for s in cls])
            it_ok &= len(U) == 3 and all(inside(U, perp([q[t]], h)) for t in opp)
            G_ = [[dot(a, b) for b in U] for a in U]; gram[dim([sum(r[j] << j for j in range(3)) for r in G_])] += 1
            for s in cls:
                for t in opp:
                    qp = perp([q[t]], h)
                    dl[(dL([q[s]], U), dL(U, qp), dL(qp, [(1 << h) - 1] and [1 << j for j in range(h)]))] += 1
    rep('U_par (4 same-parity addresses) has dim 3 and lies in q_T^perp for every opposite-parity T', it_ok)
    rep('source itinerary distances <q>->U->q_T^perp->F: %s (want (2, h-4, 1))' % dict(dl), set(dl) == {(2, h - 4, 1)})
    print('    Gram rank of U_par over F2: %s (rank 1: U_par is DEGENERATE, radical dim 2)' % dict(gram))
    # brute-force check of the Lagrangian distance formula on F2^{2h} (small h only)
    if h <= 10:
        import random
        rnd = random.Random(p); bad = 0; tried = 0
        samples = [[q[0]], basis([q[s] for s in range(4) if par[s] == 0] or [q[0]])]
        for _ in range(60):
            k1, k2 = rnd.randrange(h + 1), rnd.randrange(h + 1)
            U = basis([rnd.randrange(1 << h) for _ in range(k1)]); V = basis([rnd.randrange(1 << h) for _ in range(k2)])
            samples.append(U); samples.append(V)
        samples += [basis([q[s] for s in Sv[0]]), perp([q[0]], h), [1 << j for j in range(h)], []]
        for U in samples:
            for V in samples[:12]:
                tried += 1; bad += dL(U, V) != dL_brute(U, V, h)
        rep('Lagrangian distance formula == brute force on F2^%d (%d pairs, incl. degenerate)' % (2 * h, tried), bad == 0)
    print('p=%d: %s' % (p, 'ALL OK' if ok else 'SOME FAIL'), flush=True)
    return ok

if __name__ == '__main__':
    res = [run(int(a)) for a in sys.argv[1:]]
    print('ALL OK' if all(res) else 'SOME FAIL')
