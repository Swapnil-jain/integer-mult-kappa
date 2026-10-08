"""Exact linear algebra over Q on integer row vectors of length h (stdlib only).
A subspace is stored as its canonical basis: the RREF over Q with every row scaled to a primitive integer vector
with positive pivot (unique for the subspace), as a tuple of tuples."""
from math import gcd


def prim(r):
    g = 0
    for x in r:
        if x: g = gcd(g, x)
    if g > 1: r = [x // g for x in r]
    return r


def rref(rows, h):
    M = [list(r) for r in rows if any(r)]; piv = []; r = 0; n = len(M)
    for c in range(h):
        if r == n: break
        i = next((i for i in range(r, n) if M[i][c]), None)
        if i is None: continue
        M[r], M[i] = M[i], M[r]; R_ = M[r]
        if R_[c] < 0: R_ = M[r] = [-x for x in R_]
        R_ = M[r] = prim(R_); p = R_[c]
        for k in range(n):
            if k != r and M[k][c]:
                q = M[k][c]; g = gcd(p, q); a, b = p // g, q // g
                M[k] = prim([a * x - b * y for x, y in zip(M[k], R_)])
        piv.append(c); r += 1
    out = []
    for row in M[:r]:
        c = next(j for j, x in enumerate(row) if x)
        if row[c] < 0: row = [-x for x in row]
        out.append(tuple(prim(row)))
    return tuple(out)


def pivots(B): return [next(j for j, x in enumerate(r) if x) for r in B]


def ann(B, h):
    """canonical basis of the annihilator {z : B z = 0} of a canonical basis B."""
    if not B: return tuple(tuple(int(i == j) for j in range(h)) for i in range(h))
    pv = pivots(B); free = [j for j in range(h) if j not in pv]; out = []
    L = 1
    for r, c in zip(B, pv): L = L * r[c] // gcd(L, r[c])
    for f in free:
        z = [0] * h; z[f] = L
        for r, c in zip(B, pv): z[c] = -r[f] * (L // r[c])
        out.append(z)
    return rref(out, h)


def inside_ann(A, Z):
    """rows of A killed by every row of Z (Z = annihilator of the bigger space)."""
    return all(sum(a * z for a, z in zip(ra, rz)) == 0 for ra in A for rz in Z)


def meet_ann(ZA, ZB, h):
    return ann(rref(list(ZA) + list(ZB), h), h)


def red(B, p):
    """reduction mod p of a canonical basis, normalised to the mod-p RREF (None if a pivot vanishes mod p)."""
    out = []
    for r in B:
        c = next(j for j, x in enumerate(r) if x)
        if r[c] % p == 0: return None
        iv = pow(r[c], p - 2, p); out.append(tuple(x * iv % p for x in r))
    return tuple(out)


def rref_mod(rows, p, h):
    M = [[x % p for x in r] for r in rows]; r = 0; n = len(M)
    for c in range(h):
        if r == n: break
        i = next((i for i in range(r, n) if M[i][c]), None)
        if i is None: continue
        M[r], M[i] = M[i], M[r]; inv = pow(M[r][c], p - 2, p); M[r] = [x * inv % p for x in M[r]]
        for k in range(n):
            if k != r and M[k][c]:
                g = M[k][c]; M[k] = [(x - g * y) % p for x, y in zip(M[k], M[r])]
        r += 1
    return tuple(tuple(x) for x in M[:r])


def rank_mod(rows, p, k):
    M = [[x % p for x in r] for r in rows]; r = 0; n = len(M)
    for c in range(k):
        if r == n: break
        i = next((i for i in range(r, n) if M[i][c]), None)
        if i is None: continue
        M[r], M[i] = M[i], M[r]; inv = pow(M[r][c], p - 2, p); M[r] = [x * inv % p for x in M[r]]
        for kk in range(r + 1, n):
            if M[kk][c]:
                g = M[kk][c]; M[kk] = [(x - g * y) % p for x, y in zip(M[kk], M[r])]
        r += 1
    return r


def nondeg(B, p=2**31 - 1):
    """det(B (9I - J) B^T) != 0 mod p  =>  != 0 over Z, i.e. G = I - J/9 nondegenerate on span B."""
    if not B or len(B) == len(B[0]): return True
    sm = [sum(b) for b in B]
    Gr = [[9 * sum(x * y for x, y in zip(a, b)) - sa * sb for b, sb in zip(B, sm)] for a, sa in zip(B, sm)]
    return rank_mod(Gr, p, len(B)) == len(B)
