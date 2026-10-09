"""Binary linear algebra and binary symplectic Clifford labels for nullity-dependent sharing (own code).

Vectors of F2^m are ints. A point of the symplectic space E (+) E is the int a | (b << m) for the pair (a, b), with
pairing <(a,b),(a',b')> = a.b' + b.a'. A symplectic matrix is the list of its 2m column images.
Conventions follow #130's general-Clifford-frame note: F = [[I,I],[0,I]] (the Gaussian-dyadic C^(x)m), L_0 = {(0,z)},
frame label L(T) = T^-1 L_0, and the Fourier (Bruhat) rank of an exact Clifford X is the rank of its B block,
X(0,z) = (Bz, Dz).  Canonical frame of a nondegenerate U: S_U = [[I,P_U],[0,I]], P_U the orthogonal projector
(its label is L_U = {(u, u+v): u in U, v in U-perp}); it is the symplectic image of C_U = H diag(i^wt(p_U x)) H."""
import random

popc = lambda x: bin(x).count('1')
dot = lambda a, b: popc(a & b) & 1

def echelon(vs):
    piv = {}
    for x in vs:
        while x:
            p = x.bit_length() - 1
            if p in piv: x ^= piv[p]
            else: piv[p] = x; break
    return piv
def rank(vs): return len(echelon(vs))
def basis(vs): return list(echelon(vs).values())
def reduce(piv, x):
    while x:
        p = x.bit_length() - 1
        if p in piv: x ^= piv[p]
        else: return x
    return 0
def contains(V, x): return reduce(echelon(V), x) == 0
def perp_in(V, W, m):
    """{x in span V : x.w = 0 for all w in W}."""
    V = basis(V); W = basis(W)
    if not W: return V
    # solve: combinations c of V with (sum c_i V_i).w_j = 0
    rows = []
    for i, v in enumerate(V):
        rows.append((sum(dot(v, w) << j for j, w in enumerate(W)), v))
    # gaussian elimination on the syndrome part
    out = []; piv = {}
    for s, v in rows:
        while s:
            p = s.bit_length() - 1
            if p in piv: s ^= piv[p][0]; v ^= piv[p][1]
            else: piv[p] = (s, v); break
        if not s: out.append(v)
    return basis(out)
def gram(B): return [[dot(a, b) for b in B] for a in B]
def nondeg(B):
    B = basis(B); n = len(B)
    return rank([sum(dot(a, b) << j for j, b in enumerate(B)) for a in B]) == n
def radical(B): B = basis(B); return perp_in(B, B, 0)
def alternating(B): return not any(popc(x) & 1 for x in basis(B))

def orthonormal_basis(V):
    """orthonormal basis of a nondegenerate non-alternating space (asserted)."""
    W = basis(V); out = []
    while W:
        e = next((x for x in W if popc(x) & 1), None)
        if e is None:
            # W is alternating: borrow the last orthonormal vector o and a hyperbolic pair (x, y) in W:
            # (o, x, y) -> (o+x+y, o+x, o+y) is orthonormal and spans the same space.
            assert out, 'space is alternating'
            o = out.pop(); x = W[0]; y = next(z for z in W if dot(x, z))
            new = [o ^ x ^ y, o ^ x, o ^ y]
            assert all(dot(a, b) == (a == b) for a in new for b in new)
            out += new; W = perp_in(W, [x, y], 0); continue
        out.append(e); W = perp_in(W, [e], 0)
    assert all(dot(a, b) == (i == j) for i, a in enumerate(out) for j, b in enumerate(out))
    return out

# ---------------- matrices (column images) ----------------
def mapply(M, x):
    y = 0; k = 0
    while x:
        if x & 1: y ^= M[k]
        x >>= 1; k += 1
    return y
def mmul(A, B): return [mapply(A, c) for c in B]
def from_basis(src, dst, n):
    """the linear map sending src[i] -> dst[i] (src a basis of F2^n)."""
    piv = {}
    for k, b in enumerate(src):
        x, c = b, 1 << k
        for p in sorted(piv, reverse=True):
            if x >> p & 1: x ^= piv[p][0]; c ^= piv[p][1]
        assert x, 'dependent'
        p = x.bit_length() - 1
        for pp in list(piv):
            y, cy = piv[pp]
            if y >> p & 1: piv[pp] = (y ^ x, cy ^ c)
        piv[p] = (x, c)
    M = []
    for j in range(n):
        c = piv[j][1]; y = 0; k = 0
        while c:
            if c & 1: y ^= dst[k]
            c >>= 1; k += 1
        M.append(y)
    return M
def minv(M, n): return from_basis(M, [1 << j for j in range(n)], n)
def transpose(M, n): return [sum(((M[j] >> i) & 1) << j for j in range(n)) for i in range(n)]
def is_orth(M, n): return all(dot(M[i], M[j]) == (i == j) for i in range(n) for j in range(n))

# ---------------- symplectic blocks on E (+) E ----------------
def block(Am, Bm, Cm, Dm, m):
    """[[A,B],[C,D]]: (a,b) -> (A a + B b, C a + D b); blocks given as column lists of length m."""
    cols = []
    for j in range(m): cols.append(Am[j] | (Cm[j] << m))
    for j in range(m): cols.append(Bm[j] | (Dm[j] << m))
    return cols
def zero(m): return [0] * m
def ident(m): return [1 << j for j in range(m)]
def proj(U, m):
    """orthogonal projector onto a nondegenerate U (columns)."""
    U = basis(U); n = len(U)
    G = [sum(dot(U[i], U[j]) << j for j in range(n)) for i in range(n)]
    Gi = minv([sum(((G[r] >> c) & 1) << r for r in range(n)) for c in range(n)], n)    # columns of G^-1
    cols = []
    for j in range(m):
        y = sum(dot(U[i], 1 << j) << i for i in range(n))          # coordinates <U_i, e_j>
        c = mapply(Gi, y)
        cols.append(mapply(U, c) if c else 0)
    return cols
def shear(P, m): return block(ident(m), P, zero(m), ident(m), m)
def F_on(coords_mask, m):
    """C on the given coordinates: [[I, I_S],[0, I]]."""
    return shear([(1 << j) if coords_mask >> j & 1 else 0 for j in range(m)], m)
def canon(U, m): return shear(proj(U, m), m)
def KG(G, m):
    """#130's exact lift K_G = [[G,0],[G+G^-T, G^-T]]."""
    Git = transpose(minv(G, m), m)
    return block(G, zero(m), [a ^ b for a, b in zip(G, Git)], Git, m)
def lift_orth(k, m): return block(k, zero(m), zero(m), k, m)      # orthogonal k acts as diag(k, k)
def T_general(G, s, m):
    """#130's general frame T_U = K_G T_{E_s} K_G^-1 for U = G E_s."""
    K = KG(G, m); return mmul(mmul(K, F_on((1 << s) - 1, m)), minv(K, 2 * m))
def brank(X, m):
    """Fourier rank: rank of the B block (images of (0, e_j))."""
    return rank([X[m + j] & ((1 << m) - 1) for j in range(m)])
def label(T, m):
    Ti = minv(T, 2 * m); return basis([mapply(Ti, 1 << (m + j)) for j in range(m)])
def is_symp(X, m):
    n = 2 * m
    w = lambda x, y: dot(x & ((1 << m) - 1), y >> m) ^ dot(x >> m, y & ((1 << m) - 1))
    return all(w(X[i], X[j]) == w(1 << i, 1 << j) for i in range(n) for j in range(i + 1, n))

def support(X, m):
    """smallest Z (span) with X - I supported on Z (+) Z and trivial on Z-perp (+) Z-perp: Z contains both components of
    im(X - I) and the orthogonal complement of {x : (x,0),(0,x) in ker(X - I)}."""
    n = 2 * m; mask = (1 << m) - 1
    im = [X[j] ^ (1 << j) for j in range(n)]
    comps = basis([v & mask for v in im] + [v >> m for v in im])
    # K = {x : (X - I)(x,0) = 0 and (X - I)(0,x) = 0}
    rows = [(im[j] | (im[m + j] << n), 1 << j) for j in range(m)]
    piv = {}; ker = []
    for s, c in rows:
        while s:
            p = s.bit_length() - 1
            if p in piv: s ^= piv[p][0]; c ^= piv[p][1]
            else: piv[p] = (s, c); break
        if not s: ker.append(c)
    Kp = perp_in(ident(m), ker, m)
    return basis(comps + Kp)
