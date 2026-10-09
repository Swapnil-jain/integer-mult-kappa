"""#144-style bit (m=3h=69, one role per slot through all three stages, W=2v+R per vertex), with gauge omission
(an omitted deferred role reads its old value at frame 0: first step f->r becomes 0->r, exterior 3f disappears,
its readout levels leave the target chains), optionally with nullity-dependent sharing (a kept slot of nullity
u=h-f is shared by t=floor(m/u)/3 vertex orbits: t*u total visits... see vis()). Greedy omission selection."""
import pickle, json, os, sys, math, bisect
from collections import Counter, defaultdict
from fractions import Fraction as Q
from coverlib import certify, with_fallback, stopped, root
HERE = os.path.dirname(os.path.abspath(__file__))
import gzip
def _slots():
    with gzip.open(os.path.join(HERE, '..', '..', 'certificates', 'round9', 'w2_slots.json.gz'), 'rt') as fh:
        d = json.load(fh)
    d['rest'] = dict(enumerate(d['rest'])); d['tg'] = {s: t for s, t in d['tg']}; d['sel'] = set(d['sel'])
    return d
D = {'w2': _slots()}   # frozen witness-2 slot data (certificates/round9/w2_slots.json.gz)

def vis(u, m, share):
    """visits per physical role: 3 (cross-stage, #144) or floor(m/u) (nullity sharing; >= 3 since u <= h = m/3)"""
    return m // u if share else 3

def build(w, omitted, share=False):
    h, v, R, f, rest, tg, sel = w['h'], w['v'], w['R'], w['f'], w['rest'], w['tg'], w['sel']; m = 3 * h
    aux = Counter(); ext = Counter(); Wv = Q(2 * v); lev = defaultdict(set)
    for s in range(R):
        keep = s in sel and s not in omitted
        ds = ([f[s]] if keep else [0]) + rest[s]
        for a, b in zip(ds, ds[1:]):
            if b > a: aux[b - a] += 1
        u = h - f[s] if keep else h; t = vis(u, m, share); Wv += Q(3, t)
        if m - t * u: ext[m - t * u] += Q(3, t)
        if keep:
            for x in tg[s]: lev[x].add(f[s])
    tgt = Counter()
    for x in range(v):
        ds = sorted(set([0] + list(lev[x]) + [h - 1]))
        for a, b in zip(ds, ds[1:]): tgt[b - a] += 1
    src = Counter(r for p in w['xd'] for r in p if r)
    H = Counter()
    for part in (aux, tgt, src, Counter({h - 1: h})):
        for r, n in part.items(): H[r] += 3 * n
    H[2] += 2 * v
    for r, n in ext.items(): H[r] += n
    K = math.lcm(Wv.denominator, *[Q(n).denominator for n in H.values()])
    Hk = {r: int(n * K) for r, n in H.items() if n}; Wk = int(Wv * K); s_ = sum(r * n for r, n in Hk.items())
    return dict(m=m, W=Wk, s=s_, D=Wk * m - s_, hist=Hk, k=K, Wv=Wv, aux=aux, tgt=tgt, src=src, ext=ext)

def greedy(w, share=False, a0=4.6e-4, passes=6, start=None):
    """omission set minimising F/W at saving a0 (exact deltas incl. target-level merges); first-improvement passes."""
    h, v, R, f, rest, tg, sel = w['h'], w['v'], w['R'], w['f'], w['rest'], w['tg'], w['sel']; m = 3 * h; tau = 1 - a0
    c = lambda r: (r / m) ** tau if r > 0 else 0.0
    om = set(start or ())
    base = build(w, om, share); F = sum(n * c(r) for r, n in base['hist'].items()) / base['k']; Wv = float(base['Wv'])
    cnt = defaultdict(Counter)
    for s in sel:
        if s not in om:
            for x in tg[s]: cnt[x][f[s]] += 1
    lvls = {x: sorted(set([0, h - 1]) | set(cnt[x])) for x in range(v)}
    def slotF(s, keep):
        ds = ([f[s]] if keep else [0]) + rest[s]; u = h - f[s] if keep else h; t = vis(u, m, share)
        val = 3 * sum(c(b - a) for a, b in zip(ds, ds[1:]) if b > a) + 3 / t * c(m - t * u); return val, 3 / t
    def tdelta(s, remove):
        """F change on targets when slot s's readouts leave (remove) or join the chains"""
        d = 0.0; fs = f[s]
        for x in tg[s]:
            L = lvls[x]
            if remove:
                if cnt[x][fs] != 1 or fs in (0, h - 1): continue
                i = bisect.bisect_left(L, fs); a, b = L[i - 1], L[i + 1]
                d += 3 * (c(b - a) - c(fs - a) - c(b - fs))
            else:
                if cnt[x][fs] > 0 or fs in L: continue
                i = bisect.bisect_left(L, fs); a, b = L[i - 1], L[i]
                d += 3 * (c(fs - a) + c(b - fs) - c(b - a))
        return d
    def apply(s, remove):
        for x in tg[s]:
            if remove:
                cnt[x][f[s]] -= 1
                if cnt[x][f[s]] == 0 and f[s] not in (0, h - 1): lvls[x].remove(f[s])
            else:
                if cnt[x][f[s]] == 0 and f[s] not in (0, h - 1): bisect.insort(lvls[x], f[s])
                cnt[x][f[s]] += 1
    for p in range(passes):
        changed = 0
        for s in sorted(sel, key=lambda s: f[s]):
            keep = s not in om
            (fk, wk), (fo, wo) = slotF(s, True), slotF(s, False)
            if keep: dF = fo - fk + tdelta(s, True); dW = wo - wk
            else: dF = fk - fo + tdelta(s, False); dW = wk - wo
            if (F + dF) / (Wv + dW) < F / Wv - 1e-15:
                F += dF; Wv += dW; changed += 1
                if keep: om.add(s); apply(s, True)
                else: om.discard(s); apply(s, False)
        print(f'  pass {p}: changes {changed}, omitted {len(om)}, F/W {F / Wv:.12f}', flush=True)
        if not changed: break
    return om

AC = Q(4856569, 10**10)
def kappa144(ab):
    beta, eta = Q(1, 10**6), Q(1, 10**8)
    a = min(Q(ab), (1 - beta) * AC - Q(1, 10**10)); q0 = a * (1 - 2 * eta); c0 = q0 * (1 + eta)
    return (1 - eta) / (1 + c0 + q0) * q0

def report(tag, c):
    cc = with_fallback(c); a, r = certify(cc); ab = stopped(a); k = kappa144(ab)
    print(f'{tag:58s} W/v={float(c["Wv"]):.4f} D/vertex={float(Q(c["D"], c["k"]))} a*={a} ({float(a):.9e}) a_bit={float(ab):.10e} kappa={float(k):.10e} floor1e-10={math.floor(k * 10**10)}', flush=True)
    return dict(tag=tag, astar=str(a), a_bit=str(ab), a_bit_f=float(ab), kappa=str(k), kappa_f=float(k), maxchild=max(c['hist']), W_per_vertex=str(c['Wv']))

if __name__ == '__main__':
    res = []
    sel144 = json.load(open(os.path.join(HERE, 'in', 'pr144_selection.json')))
    w1, w2 = D['w1'], D['w2']
    c = build(w1, set(sel144['omitted_readout_order']))
    P = json.load(open(os.path.join(HERE, 'in', 'pr144_bit_input.json')))
    theirs = {int(r): n for r, n in eval(P['child_histogram']).items()} if isinstance(P['child_histogram'], str) else {int(r): n for r, n in P['child_histogram'].items()}
    print('w1 + #144 selection: child histogram equals #144:', c['hist'] == theirs and c['k'] == 1, 'W', c['W'], 'rank', c['s'], 'D', c['D'])
    res.append(report('w1 = PR97, #144 selection (reproduction)', c))
    res.append(report('w1, no omission (= #137)', build(w1, set())))
    for name, w in (('w1', w1), ('w2', w2)):
        for share in (False, True):
            print(f'greedy {name} share={share}', flush=True)
            om = greedy(w, share)
            r = report(f'{name}, greedy omission ({len(om)}), share={share}', build(w, om, share)); r['omitted'] = len(om); res.append(r)
            pickle.dump(sorted(om), open(os.path.join(HERE, f'omit_{name}_{int(share)}.pkl'), 'wb'))
    res.append(report('w2, no omission, share=True', build(w2, set(), True)))
    json.dump(res, open(os.path.join(HERE, 'omit.json'), 'w'), indent=1)
