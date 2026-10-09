"""Word loader + ledger for frozen paired-cube complex words (check_word2 JSON). Own code, read-only on the word.
The ledger reproduces check_word2.ledger exactly (verified against the claim), and exposes the per-role chains so
accounting variants (collective reuse, priced non-monotone steps) can be recounted from the same objects."""
import gzip, json, math
from fractions import Fraction as Q
from collections import Counter, defaultdict

def basis(rows):
    b = {}
    for x in rows:
        for q in sorted(b, reverse=True):
            if x >> q & 1: x ^= b[q]
        if x:
            q = x.bit_length() - 1
            for k in list(b):
                if b[k] >> q & 1: b[k] ^= x
            b[q] = x
    return tuple(b[q] for q in sorted(b, reverse=True))

def inside(A, B):
    piv = {y.bit_length() - 1: y for y in basis(B)}
    for x in A:
        for q in sorted(piv, reverse=True):
            if x >> q & 1: x ^= piv[q]
        if x: return False
    return True

def meet(A, B, h):
    # A cap B = perp(perp A + perp B)
    return perp(perp(A, h) + perp(B, h), h)

def perp(rows, h):
    rows = basis(rows); piv = {r.bit_length() - 1: r for r in rows}; out = []
    for j in range(h):
        if j in piv: continue
        x = 1 << j
        for q, r in piv.items():
            if r >> j & 1: x |= 1 << q
        out.append(x)
    return basis(out)

class Word:
    def __init__(self, path):
        J = json.load(gzip.open(path, 'rt')); self.J = J
        self.h, self.v, self.m, self.R = J['h'], J['v'], J['m'], J['R']; h = self.h
        self.args, self.signs, self.inputs = J['args'], J['signs'], J['inputs']
        self.roots = J['roots']; self.ops = [tuple(o) for o in J['ops']]; self.sched = J['sched']; self.cut = J['cut']
        self.sources = {int(x): s for x, s in J['sources'].items()}; self.rootroles = J['rootroles']
        self.gauges = J['gauges']; self.tau = {int(b): t for b, t in J['tau'].items()}
        self.pairs = [tuple(q) for q in J['pairs']]; self.frames = [tuple(f) for f in J['frames']]
        self.sinks = {z['role'] for z in J.get('sinks', [])}
        self.FULL = tuple(1 << j for j in range(h - 1, -1, -1))
        n = len(self.args); self.n = n; sp = [()] * n
        for x in range(1, n): sp[x] = (self.inputs[x - 1],) if self.args[x] is None else basis(sp[self.args[x][0]] + sp[self.args[x][1]])
        self.spans = sp
        self.pos = {i: k for k, i in enumerate(self.sched)}
        self.rops = [[] for _ in range(self.R)]
        for i in self.sched: a, b, _ = self.ops[i]; self.rops[a].append(i); self.rops[b].append(i)
        # value (node id) held by each role after each op; None = fresh dirty start
        hold = {sr: xx for xx, sr in self.sources.items()}; self.after = {}
        for i in self.sched:
            d, c, xx = self.ops[i]; hold[d] = xx; self.after[i] = xx
        self.final = hold
        self.rootframe, self.rootkind = {}, {}
        for r, sr in zip(self.roots, self.rootroles):
            self.rootkind[sr] = r['kind']
            if r['kind'] == 'center': self.rootframe[sr] = sp[r['node']]
            else: self.rootframe[sr] = perp(basis([self.inputs[t] for t in r['targets']]), h)
        self.start = [()] * self.R
        for xx, sr in self.sources.items(): self.start[sr] = (self.inputs[xx - 1],)
        self.gsig = {}
        for z in self.gauges: self.gsig[z['role']] = perp(z['A'], h); self.start[z['role']] = self.gsig[z['role']]
        self.donor = {a: b for a, b in self.pairs}; self.recip = {b: a for a, b in self.pairs}

    def chain(self, sr):
        return [self.start[sr]] + [self.frames[i] for i in self.rops[sr]] + ([self.rootframe[sr]] if sr in self.rootframe else [])

    def slots(self):
        """physical slots: list of spliced role sequences (sinks dropped), as in check_word2.ledger"""
        out = []
        for sr in range(self.R):
            if sr in self.recip or sr in self.sinks: continue
            seq = [sr]; cur = self.donor.get(sr)
            while cur is not None: seq.append(cur); cur = self.donor.get(cur)
            out.append(seq)
        return out

    def target_hist(self):
        """target-chain children exactly as check_word2.ledger with the frozen sinks (pre/post shears, pivot writes)"""
        h, v = self.h, self.v; S = {z['role']: z for z in self.J.get('sinks', [])}
        ev = defaultdict(list)
        for k, z in enumerate(self.gauges):
            A = basis(z['A'])
            for t in z['targets']: ev[t].append(((self.tau[z['role']], 0, len(self.gauges) - 1 - k), A))
        rootidx = {sr: j for j, sr in enumerate(self.rootroles)}
        for s, z in S.items():
            T = z['targets']; L = self.pos[self.rops[s][-1]]
            ann = basis([self.inputs[t] for t in self.roots[rootidx[s]]['targets']])
            if len(T) > 1:
                for t in T: ev[t].append(((self.cut, -1, 0), self.FULL))
            for i in self.rops[s]: ev[z['pivot']].append(((self.pos[i], 1, 0), perp(self.frames[i], h)))
            for t in T: ev[t].append(((L, 2, 0), ann))
        for sr, j in rootidx.items():
            if self.rootkind[sr] == 'side' and sr not in S:
                ann = basis([self.inputs[t] for t in self.roots[j]['targets']])
                for t in self.roots[j]['targets']: ev[t].append(((float('inf'), j, 0), ann))
        tgt = Counter(); bad = 0
        for t in range(v):
            cur = self.FULL
            for key, A in sorted(ev.get(t, []), key=lambda e: e[0]):
                if not inside(A, cur): bad += 1
                tgt[len(cur) - len(A)] += 1; cur = A
            tgt[len(cur) - 1] += 1
        return tgt, bad

    def local_hist(self, seqfn=None):
        """aux/source role children: per slot, the nested chain of its tenants, increments counted.
        seqfn(slot) may return a replacement list of frames (for accounting variants); steps are then charged
        as dim(B) - dim(A cap B) (free descent, paid climb), which equals d1-d0 on nested steps."""
        h = self.h; local = Counter(); bad = 0
        for slot in self.slots():
            seq = []
            for sr in slot: seq += self.chain(sr)
            seq.append(self.FULL)
            if seqfn: seq = seqfn(slot, seq)
            sr = slot[0]
            if len(self.start[sr]) == 1 and sr not in self.gsig: local[1] += 1
            for A_, B_ in zip(seq, seq[1:]):
                if inside(A_, B_):
                    if len(B_) > len(A_): local[len(B_) - len(A_)] += 1
                else:
                    bad += 1; k = len(B_) - len(meet(A_, B_, h))
                    if k: local[k] += 1
            for r in slot:
                if self.rootkind.get(r) == 'center': local[len(self.rootframe[r])] += 1
        return local, bad

    def ledger(self, local=None, tgt=None, dW=0, extraC=None):
        h, v = self.h, self.v
        if local is None: local, _ = self.local_hist()
        if tgt is None: tgt, _ = self.target_hist()
        C = Counter()
        for hist in (local, {1: v, 2: v, h - 4: v}, tgt):
            for r_, k in hist.items():
                if r_: C[r_] += 3 * k
        C[2] += 2 * v
        for z in self.gauges:
            if z['role'] not in self.recip: C[3 * (h - len(basis(z['A'])))] += 1
        for r_, k in (extraC or {}).items(): C[r_] += k
        Wv = 2 * v + self.R - len(self.pairs) - len(self.sinks) + dW
        return +C, Wv

def froot(C, W, m):
    """float root a of (1/W) sum n (w/m)^(1-a) = 1"""
    def F(a): return math.fsum(n * (w / m) ** (1 - a) for w, n in C.items() if n) / W
    lo, hi = 0.0, 0.05
    for _ in range(200):
        mid = (lo + hi) / 2; lo, hi = (mid, hi) if F(mid) < 1 else (lo, mid)
    return lo

def g(r, a, m): return r * math.expm1(a * math.log(m / r)) if r > 0 else 0.0
