"""Physical (per-operation) frame descent on the paired word, with the pairs and gauges fixed (own code).

#168 moves operation frames by "physical descent in connected bundles of operations that share a frame, chosen by
the charged moment". Our node-level descent (descent.py) runs before gauges and pairs exist, on #144's node ledger.
Here the variables are the frames of the individual operations and the objective is exactly what recount() charges:
every role's spliced chain  start -> op frames -> root frame -> [recipient gauge -> its op frames -> ...] -> full
pays phi(step) for every positive dimension step, phi(t) = t * expm1(a ln(m/t)) (the moment surplus at a).

Legal frames: an operation i on roles (a, b) at node x needs span(x) <= U_i and, on both chains, prev <= U_i <= next.
Two move types per sweep:
  * single op:  U_i anywhere in [lo, hi] (lo = span + both prevs, hi = both nexts), as lo + a prefix of a fixed
    complement (all dims tried);
  * bundle:     a connected set of ops joined by chain links whose two frames are EQUAL moves as one frame. Inside a
    bundle every step is 0 and a single op cannot move (its lo and hi both equal the shared frame), so this is the only
    way to shrink or grow a run of equal frames. lo/hi/cost use only the links leaving the bundle.
Gotcha: the donor -> recipient splice is a chain link like any other, so the donor's last frame stays inside the
recipient gauge automatically; recount() re-checks every nesting from scratch, so nothing here is trusted."""
import math, sys, os
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cube_word import basis, perp, contained
from cbirth import spans_of
SUB = os.environ.get('PD_SUB') == '1'
SEG = int(os.environ.get('PD_SEG', '0'))

def build_seqs(W, F, pairs):
    R = W['pr']['R']; h = W['h']
    donor = {a: b for a, b in pairs}; recip = {b: a for a, b in pairs}
    full = tuple(1 << j for j in range(h - 1, -1, -1))
    seqs = []
    for sr in range(R):
        if sr in recip: continue
        seq = []; cur = sr
        while cur is not None:
            seq.append(('F', F['start'][cur]))
            seq += [('O', i) for i in F['rops'][cur]]
            if cur in F['rootframe']: seq.append(('F', F['rootframe'][cur]))
            cur = donor.get(cur)
        seq.append(('F', full))
        seqs.append(seq)
    return seqs

def run(W, F, pairs, a=6.1e-4, sweeps=8, verbose=True, log=print):
    h, m = W['h'], W['m']
    phi = lambda t: t * math.expm1(a * math.log(m / t)) if t > 0 else 0.0
    frames = list(F['frames']); ops = F['ops']; sp = spans_of(W['wit']['g'])
    seqs = build_seqs(W, F, pairs)
    # links: for each op, list of (prev item, next item); an item is ('F', frame) or ('O', j)
    prevs = defaultdict(list); nexts = defaultdict(list)
    for seq in seqs:
        for k, it in enumerate(seq):
            if it[0] == 'O':
                i = it[1]; prevs[i].append(seq[k - 1]); nexts[i].append(seq[k + 1])
    nops = len(ops)
    assert all(len(prevs[i]) == 2 for i in range(nops) if prevs[i]), 'every op sits on two chains'
    fr = lambda it: it[1] if it[0] == 'F' else frames[it[1]]
    def cost_all():
        c = 0.0
        for seq in seqs:
            d = [len(fr(it)) for it in seq]
            c += sum(phi(y - x) for x, y in zip(d, d[1:]) if y > x)
        return c
    def inter(fs):
        return perp(basis(tuple(z for f in fs for z in perp(f, h))), h)
    def best_frame(lo, hi, P, N):
        """P, N: dims of external prev / next frames; returns (delta-free cost, frame) minimising the link cost"""
        comp = []; cur = list(lo); r0 = len(lo)
        for z in hi:
            if len(basis(cur + [z])) > len(cur): comp.append(z); cur = list(basis(cur + [z]))
        best = None
        for k in range(len(comp) + 1):
            r = r0 + k
            c = sum(phi(r - p) for p in P) + sum(phi(n - r) for n in N)
            if best is None or c < best[0] - 1e-12: best = (c, k)
        return best[0], (basis(tuple(lo) + tuple(comp[:best[1]])) if best[1] else tuple(lo))
    tot = cost_all()
    if verbose: log('pdescent start: ops %d seqs %d cost %.4f' % (nops, len(seqs), tot))
    for sw in range(sweeps):
        moved = 0; gain = 0.0
        # bundles: union ops along links with equal frames
        par = list(range(nops))
        def find(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for i in range(nops):
            for it in nexts[i]:
                if it[0] == 'O' and frames[it[1]] == frames[i]:
                    a_, b_ = find(i), find(it[1])
                    if a_ != b_: par[a_] = b_
        comps = defaultdict(list)
        for i in range(nops):
            if prevs[i]: comps[find(i)].append(i)
        groups = list(comps.values()) + [[i] for c in comps.values() if len(c) > 1 for i in c]
        if SUB:
            # sub-bundles: inside a component, the ops reachable from i through equal-frame prev links (its
            # in-bundle ancestors, which can shrink with i) and through next links (descendants, which can grow)
            for c in comps.values():
                if len(c) < 3: continue
                cs = set(c)
                for i in c:
                    for links in (prevs, nexts):
                        seen = {i}; st = [i]
                        while st:
                            z = st.pop()
                            for it in links[z]:
                                if it[0] == 'O' and it[1] in cs and it[1] not in seen: seen.add(it[1]); st.append(it[1])
                        if 1 < len(seen) < len(c): groups.append(sorted(seen))
        for B in groups:
            Bs = set(B)
            if any(frames[i] != frames[B[0]] for i in B): continue   # a member moved earlier this sweep
            U = frames[B[0]]
            Pf = [fr(it) for i in B for it in prevs[i] if not (it[0] == 'O' and it[1] in Bs)]
            Nf = [fr(it) for i in B for it in nexts[i] if not (it[0] == 'O' and it[1] in Bs)]
            lo = basis(tuple(z for i in B for z in sp[ops[i][2]]) + tuple(z for f in Pf for z in f))
            hi = inter(Nf) if Nf else tuple(1 << j for j in range(h - 1, -1, -1))
            if not contained(lo, hi): continue
            P = [len(f) for f in Pf]; N = [len(f) for f in Nf]; r = len(U)
            c0 = sum(phi(r - p) for p in P) + sum(phi(n - r) for n in N)
            c1, U1 = best_frame(lo, hi, P, N)
            if c1 < c0 - 1e-9 and U1 != U:
                for i in B: frames[i] = U1
                moved += len(B); gain += c1 - c0
        if SEG > 1:
            # segment merges: a run of 2..SEG consecutive ops on one chain, with DIFFERENT frames, set to one common
            # frame U in [span + external prevs, intersection of external nexts]. This forms new bundles (internal
            # steps become 0), which single-op and equal-bundle moves cannot do. Current cost counts every incident
            # link once: external links from both sides, internal links once (from the destination's prevs).
            for seq in seqs:
                opsq = [k for k, it in enumerate(seq) if it[0] == 'O']
                for a_ in range(len(opsq)):
                    for L in range(2, SEG + 1):
                        if a_ + L > len(opsq) or opsq[a_ + L - 1] - opsq[a_] != L - 1: break
                        B = [seq[k][1] for k in opsq[a_:a_ + L]]; Bs = set(B)
                        if len(Bs) < L or len({frames[i] for i in B}) == 1: continue
                        Pf = [fr(it) for i in B for it in prevs[i] if not (it[0] == 'O' and it[1] in Bs)]
                        Nf = [fr(it) for i in B for it in nexts[i] if not (it[0] == 'O' and it[1] in Bs)]
                        lo = basis(tuple(z for i in B for z in sp[ops[i][2]]) + tuple(z for f in Pf for z in f))
                        hi = inter(Nf) if Nf else tuple(1 << j for j in range(h - 1, -1, -1))
                        if not contained(lo, hi): continue
                        c0 = sum(phi(len(frames[i]) - len(fr(it))) for i in B for it in prevs[i])
                        c0 += sum(phi(len(fr(it)) - len(frames[i])) for i in B for it in nexts[i] if not (it[0] == 'O' and it[1] in Bs))
                        c1, U1 = best_frame(lo, hi, [len(f) for f in Pf], [len(f) for f in Nf])
                        if c1 < c0 - 1e-9:
                            for i in B: frames[i] = U1
                            moved += len(B); gain += c1 - c0
        tot2 = cost_all()
        if verbose: log('pdescent sweep %d: moved %d ops, predicted %.4f, cost %.4f -> %.4f' % (sw, moved, gain, tot, tot2))
        tot = tot2
        if not moved: break
    return frames

def with_frames(F, frames):
    F2 = dict(F); F2['frames'] = frames; return F2
