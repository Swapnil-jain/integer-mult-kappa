"""Own GF(2) linear algebra for the gate (written independently of independent/round8-coreshare/lib/frames.py).
Vectors are Python ints (bit i = coordinate i). A subspace is any list of spanning vectors; key() gives a canonical
form (fully reduced row echelon keyed by the highest set bit), so two spans are equal iff their keys are equal."""


def dot(a, b): return bin(a & b).count('1') & 1
def pc(a): return bin(a).count('1')


def reduce_basis(vs):
    """fully reduced echelon basis {top bit: vector}."""
    piv = {}
    for x in vs:
        for t in sorted(piv, reverse=True):
            if x >> t & 1: x ^= piv[t]
        if x:
            t = x.bit_length() - 1
            for s in piv:
                if piv[s] >> t & 1: piv[s] ^= x
            piv[t] = x
    return piv


def key(vs): return tuple(sorted(reduce_basis(vs).values()))
def dim(vs): return len(reduce_basis(vs))


def contains(big, small):
    piv = reduce_basis(big)
    for x in small:
        for t in sorted(piv, reverse=True):
            if x >> t & 1: x ^= piv[t]
        if x: return False
    return True


def kernel_combos(rows, n):
    """rows: list of n-bit ints (a matrix with len(rows) rows, n columns).  Return a basis of {c in F2^len(rows) :
    sum_i c_i rows_i = 0}, as ints over row indices."""
    aug = [(r, 1 << i) for i, r in enumerate(rows)]
    out = []; piv = {}
    for r, c in aug:
        for t in sorted(piv, reverse=True):
            if r >> t & 1: r ^= piv[t][0]; c ^= piv[t][1]
        if r: piv[r.bit_length() - 1] = (r, c)
        else: out.append(c)
    return out


def perp_within(V, U):
    """basis of span(V) intersect U^perp."""
    V = list(reduce_basis(V).values()); U = list(U)
    rows = [sum(dot(v, u) << j for j, u in enumerate(U)) for v in V]
    res = []
    for c in kernel_combos(rows, len(U)):
        x = 0
        for i, v in enumerate(V):
            if c >> i & 1: x ^= v
        res.append(x)
    return res


def gram_rank(B):
    return len(reduce_basis([sum(dot(a, b) << j for j, b in enumerate(B)) for a in B]))


def nondegenerate(B):
    B = list(reduce_basis(B).values()); return gram_rank(B) == len(B)


def has_odd(B): return any(pc(x) & 1 for x in B)       # a span has an odd vector iff some basis vector is odd


def orthonormal(B):
    """an orthonormal basis of a nondegenerate nonalternating span (asserted), else None if alternating."""
    W = list(reduce_basis(B).values())
    assert gram_rank(W) == len(W), 'degenerate'
    if W and not has_odd(W): return None
    out = []
    while W:
        e = next((x for x in W if pc(x) & 1), None)
        if e is not None:
            out.append(e); W = perp_within(W, [e]); continue
        # W alternating: take a hyperbolic pair u, w (u.w = 1) and merge with the last odd vector f:
        # {f+u, f+w, f+u+w} is orthonormal and spans <f, u, w>.
        f = out.pop(); u = W[0]; w = next(x for x in W if dot(u, x))
        out += [f ^ u, f ^ w, f ^ u ^ w]; W = perp_within(W, [u, w])
    n = len(out)
    assert all(dot(out[i], out[j]) == (i == j) for i in range(n) for j in range(n)), 'not orthonormal'
    assert dim(out) == n == dim(B), 'basis does not span'
    assert contains(B, out)
    return out


def gram_inverse(B):
    """inverse of the Gram matrix of basis B over F2, rows as ints (bit j = column j)."""
    n = len(B)
    G = [sum(dot(B[i], B[j]) << j for j in range(n)) for i in range(n)]
    A = [(G[i], 1 << i) for i in range(n)]
    for col in range(n):
        k = next(r for r in range(col, n) if A[r][0] >> col & 1)
        A[col], A[k] = A[k], A[col]
        for r in range(n):
            if r != col and A[r][0] >> col & 1: A[r] = (A[r][0] ^ A[col][0], A[r][1] ^ A[col][1])
    return [A[i][1] for i in range(n)]


def wt_point(B, eta):
    """wt(p_U eta) mod 4 for nondegenerate span B, eta a vector (exact, scalar version)."""
    B = list(reduce_basis(B).values())
    if not B: return 0
    Gi = gram_inverse(B); y = [dot(b, eta) for b in B]
    p = 0
    for i, b in enumerate(B):
        c = 0
        for j in range(len(B)):
            if Gi[i] >> j & 1: c ^= y[j]
        if c: p ^= b
    return pc(p) % 4


def kron(a, b, h):
    """a (x) b with bit index i*h + j for a_i b_j."""
    x = 0; i = 0
    while a:
        if a & 1: x |= b << (i * h)
        a >>= 1; i += 1
    return x
