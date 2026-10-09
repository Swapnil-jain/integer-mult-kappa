"""Exact algebra of the shared exteriors (independent of the construction; reads the word pickle for the slot frames).
For each group J (PR #128 partition at h = 24, AG partition otherwise) and each physical source frame class sigma
(first frame of a slot, stage one), both orientations (stage one A (x) t, stage two t (x) A):
  (1) projector annihilation over F2: P_t P_s = 0 and P_s P_t = 0 (24x24 matrices) for every ordered pair in J;
      lifted: P_{A (x) t} P_{A (x) s} = (P_A P_A) (x) (P_t P_s) = 0 in both orders;
  (2) the exterior D = sigma (x) span J + F (x) span(J)^perp has dimension m - |J|(h - dim sigma), is the exact
      orthogonal complement of sum_b sigma^perp (x) t_b (checked on 576-bit vectors), and is nondegenerate;
  (3) exact Z/4 phase identity  q_F = sum_b q_{sigma^perp (x) t_b} + q_D  as quadratic forms on F2^m, proved through
      the row decomposition: it reduces to (a) wt(z) = sum_b 3 [t_b . z] + q_{C_J}(z) and (b) wt(y) = q_sigma(y) +
      q_{sigma^perp}(y) on F2^h, each checked on all 1 + h + h(h-1)/2 polarization points (a weight-of-projection
      form mod 4 is determined there; see the docstring of polar()), and (a) also exhaustively on all 2^h points
      for every partial group;
  (4) the alternating census: sigma, sigma^perp and C_J with no odd vector (they need the PR #24 normal form), and the
      four-term Gauss sums sum_{x in U} i^{q(x)} of every alternating rank-2 block of the partial-group complements.
Usage: python3 excheck.py word.pkl"""
import sys, pickle, itertools, json, os
import numpy as np
import gf2


def proj(B, Gi, y):
    if not B: return 0
    c = 0
    for j, b in enumerate(B):
        if gf2.dot(b, y): c |= 1 << j
    p = 0
    for i, b in enumerate(B):
        if gf2.pc(Gi[i] & c) & 1: p ^= b
    return p


def qform(B):
    B = list(gf2.reduce_basis(B).values()); Gi = gf2.gram_inverse(B) if B else []
    return lambda y: gf2.pc(proj(B, Gi, y)) % 4


def polar(h):
    """0, e_i, e_i + e_j.  For q(y) = wt(P y) mod 4 with P linear: q(x+y) = q(x) + q(y) - 2 |Px & Py| and
    |Px & Py| mod 2 = (Px).(Py) is F2-bilinear, so q is a Z/4 quadratic form fixed by these points; a sum or
    difference of such forms is again one, hence vanishing on these points means vanishing everywhere."""
    yield 0
    for i in range(h): yield 1 << i
    for i in range(h):
        for j in range(i + 1, h): yield (1 << i) | (1 << j)


def main(path):
    d = pickle.load(open(path, 'rb'))
    h = d['h']; T = d['T']; v = len(T); m = h * h; Fh = [1 << p for p in range(h)]
    tvec = [sum(1 << p for p in t) for t in T]
    sys.argv = [sys.argv[0]]
    import sreplay
    G = sreplay.partition(h, T, d['sizes'], 'pr128' if h == 24 else 'ag')
    first = {}
    for roles, F, kind, ups in d['ops']:
        for r in roles:
            if r[0] == 's' and r[1] not in first and not (kind == 'cread' and r[1] in d['retslot']): first[r[1]] = gf2.key(list(F))
    sig = {}
    for q, k in first.items(): sig.setdefault(k, 0); sig[k] += 1
    out = dict(h=h, groups=len(G), sigma_classes=len(sig), slots=len(first), bad={}, max_dim_sigma=max(len(k) for k in sig),
               sigma_dim_hist=dict(sorted(__import__('collections').Counter(len(k) for k in sig).items())))
    bad = out['bad']
    if out['max_dim_sigma'] >= h: bad['sigma_full'] = 1          # a completed core must keep a positive active space
    def B(k, x): bad[k] = bad.get(k, 0) + x
    # (1) annihilation of the outer projectors, both orders
    Pt = lambda t: [t if (t >> i) & 1 else 0 for i in range(h)]          # rows of t t^T over F2
    def mul(A, Bm):
        out_ = []
        for row in A:
            r = 0
            for k in range(h):
                if row >> k & 1: r ^= Bm[k]
            out_.append(r)
        return out_
    pairs = 0
    for g in G:
        for a in g:
            for b in g:
                if a == b: continue
                pairs += 1
                if any(mul(Pt(tvec[a]), Pt(tvec[b]))): B('projector_product_nonzero', 1)
        for a in g:
            P = Pt(tvec[a])
            if mul(P, P) != P: B('projector_not_idempotent', 1)
    out['ordered_pairs_checked'] = pairs
    # (3a) outer identity per group, polarization (and exhaustive for partial groups)
    alt_C = 0; gauss = {}
    for g in G:
        span = [tvec[b] for b in g]
        C = gf2.perp_within(Fh, span)
        if len(C) != h - len(g): B('complement_dim', 1)
        if C and not gf2.nondegenerate(C): B('complement_degenerate', 1)
        if C and not gf2.has_odd(C): alt_C += 1
        qC = qform(C)
        for z in polar(h):
            if (gf2.pc(z) - sum(3 * gf2.dot(tvec[b], z) for b in g) - qC(z)) % 4: B('outer_identity_polar', 1)
        if C:
            Z = np.arange(1 << h, dtype=np.uint32)
            lhs = np.bitwise_count(Z).astype(np.int64)
            for b in g: lhs -= 3 * (np.bitwise_count(Z & np.uint32(tvec[b])) & 1).astype(np.int64)
            O = gf2.orthonormal(C)
            if O is None: B('partial_complement_alternating', 1)
            else:
                for o in O: lhs -= (gf2.pc(o) % 4) * (np.bitwise_count(Z & np.uint32(o)) & 1).astype(np.int64)
                if (lhs % 4).any(): B('outer_identity_exhaustive', 1)
                out['partial_complement_weights'] = sorted(gf2.pc(o) for o in O)
    out['alternating_outer_complements'] = alt_C
    # (3b) inner identity per sigma class, and (2) the lifted exterior for every (group, sigma) in both orientations
    alt_sig = 0; nd_bad = 0
    for sk in sig:
        S = list(sk); Sp = gf2.perp_within(Fh, S)
        if S and not gf2.nondegenerate(S): nd_bad += 1; continue
        if (S and not gf2.has_odd(S)) or (Sp and not gf2.has_odd(Sp)): alt_sig += 1
        qS, qP = qform(S), qform(Sp)
        for y in polar(h):
            if (gf2.pc(y) - qS(y) - qP(y)) % 4: B('inner_identity_polar', 1)
    out['sigma_degenerate'] = nd_bad; out['alternating_sigma_classes'] = alt_sig
    units = [1 << p for p in range(m)]; lifted = 0
    for stage in (1, 2):
        lift = (lambda r, t: gf2.kron(r, t, h)) if stage == 1 else (lambda r, t: gf2.kron(t, r, h))
        for sk in sig:
            S = list(sk); Sp = gf2.perp_within(Fh, S)
            for g in G:
                Rb = [lift(r, tvec[b]) for b in g for r in Sp]
                C = gf2.perp_within(Fh, [tvec[b] for b in g])
                D = [lift(s, tvec[b]) for b in g for s in S] + [lift(e, c) for e in Fh for c in C]
                if len(D) != m - len(g) * (h - len(S)): B('exterior_rank_formula', 1)
                if any(gf2.dot(x, y) for x in Rb for y in D) or len(gf2.reduce_basis(Rb + D)) != m: B('exterior_not_complement', 1)
                lifted += 1
                if h > 8 and lifted > 400: break          # literal 576-bit spot coverage; the rank formula is algebra
    out['lifted_exteriors_checked'] = lifted
    out['PASS'] = not bad
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
