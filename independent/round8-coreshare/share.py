"""Completed-core scratch sharing (mechanism of upstream PR #128, reimplemented from its description) on our complex
two-stage words.

Per stage, the cores of an orthogonal triple group J (Gram I) run consecutively on ONE bank of R streams. Each core is
individually transparent on arbitrary dirty scratch (checked by our exact replays), and its net physical action on
stream q is C_{F (x) t_b} C_{sigma_q (x) t_b}^{-1} (chain from its first frame sigma_q to F). All C_U are diagonal in
one Walsh basis, so they commute, and for orthogonal nondegenerate U, V: C_U C_V = C_{U+V} with the exact phase
(pairwise support overlaps are even). The group's product is C_full C_Sigma^{-1} C_{G_J^perp}^{-1}, with
Sigma = sigma_q (x) span(J) (dim |J| sigma) and G_J^perp the complement (dim (24-|J|)*24). The missing transform
C_{Sigma + G_J^perp} is applied ONCE per (stage, group, stream) at a core boundary, where the stream holds only
scratch (no computed value), so any extra commuting transform there is absorbed by the next core's dirty contract.
It is merged into one core's last step X_last -> F when that step is the stream's last op in the core (no op at F
after it): one child of rank e + X (orthogonal sum, nondegenerate). Otherwise it is its own child of rank X.
The per-core wrap child (F -> sigma + full, rank m - h + sigma) disappears. Stage two is the complement time
reversal: the same multiset (the merge sits on the core's first step instead).
Data chains and copied-centre copies are unchanged; W = 2N + 2 * (#groups) * R.
Gotcha: the merge is only legal at the 'high' step if no later op of the SAME core touches the stream at F (V-leaf
gates and source xcopies at F): an extra transform there would move a live computed value out of its frame."""
import math
from collections import Counter, defaultdict

GROUPS = [24] * 83 + [8] * 4


def slot_chains(B):
    X = B['X']; R = B['R']; fr = defaultdict(list)
    for op in B['ops']:
        roles, F, kind = op[:3]
        if kind == 'cread':
            for r in roles:
                if r[0] != 's': fr[r].append(())
            continue
        for r in roles: fr[r].append(F)
    return [fr[('s', q)] for q in range(R)]


def lamf(r, m): return r * math.log(m / r) if r else 0.0


def shared(B, w, groups=GROUPS, mode='merge'):
    """w = c7.walk(B) output. Returns new (H, W, s) for the shared allocation."""
    X = B['X']; h = B['h']; m = B['m']; v = B['v']; R = B['R']; F = X.F
    assert sum(groups) == v
    chs = slot_chains(B)
    Hs = Counter(); per = []
    for ch in chs:
        steps = [abs(len(b) - len(a)) for a, b in zip(ch, ch[1:]) if a != b]
        sig = len(ch[0]); wrap = m - h + sig
        for r in steps: Hs[r] += 1
        Hs[wrap] += 1
        merge = len(ch) >= 2 and ch[-1] == F and ch[-2] != F
        if merge and not X.le(ch[0], ch[-2]): merge = False      # need sigma inside X_last for an orthogonal sum
        e = len(F) - len(ch[-2]) if merge else 0
        per.append((steps, sig, merge, e))
    H = Counter(w['H'])
    for r, n in Hs.items():
        H[r] -= 2 * v * n
        assert H[r] >= 0
    stat = Counter()
    for k in groups:
        comp = (h - k) * h
        for steps, sig, merge, e in per:
            for r in steps: H[r] += 2 * k
            Xr = k * sig + comp
            if Xr == 0: stat['nofix'] += 1; continue
            if merge and mode == 'merge':
                H[e] -= 2; H[e + Xr] += 2; stat['merged'] += 1
            else:
                H[Xr] += 2; stat['own'] += 1
    H = {r: n for r, n in sorted(H.items()) if n}
    assert all(n > 0 for n in H.values()) and max(H) <= m - 1
    N = v * v
    Wn = 2 * N + 2 * len(groups) * R
    s = sum(r * n for r, n in H.items())
    return dict(H=H, W=Wn, s=s, D=Wn * m - s, stat=dict(stat), maxrank=max(H),
                defer_slots=sum(1 for p in per if p[1]), merge_slots=sum(1 for p in per if p[2]))


def est(H, W, m, s):
    lam = sum(n * r * math.log(m / r) for r, n in H.items()); return (W * m - s) / lam
