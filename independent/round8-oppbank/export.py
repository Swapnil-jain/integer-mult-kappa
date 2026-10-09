"""Gate export: re-run check_design.py's exact head (Q1-Q3) on witness 2 (frozen in certificates/round7/), then
collect EVERY charged h-space residual step (all ranks) and every aux sigma, verify O1 exactly over Z
(nesting A <= B with a verified annihilator; dims), and write the frames + steps for the C checker.
Also rebuilds the one-run histogram independently from check_design's exact objects (not c7.build)."""
import sys, os, time, pickle, struct
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, '..', '..'))
OUT = os.environ.get('GATE_OUT', '/tmp/round8-bitgate'); os.makedirs(OUT, exist_ok=True)
JFS = os.path.join(REPO, 'independent', 'joint-frame-stack')
src = open(os.path.join(JFS, 'check_design.py')).read().split('\n')
cut = next(i for i, l in enumerate(src) if l.startswith("report('Q5. per-role chain rank"))
code = '\n'.join(src[:cut])
sys.path[:0] = [JFS, os.path.join(JFS, '..', 'deferred-readout'), os.path.join(JFS, '..', 'two-stage-bit'), os.path.join(REPO, 'scripts')]
NS = dict(__file__=os.path.join(JFS, 'check_design.py'), __name__='cd')
exec(compile(code, 'check_design_head', 'exec'), NS)
assert all(NS['OK'].values()), NS['OK']
log = NS['log']
h = NS['h']; v = NS['v']; trip = NS['trip']; R = NS['R']; f = NS['f']; sel = NS['sel']; SIG = NS['SIG']
walk = NS['walk']; key_frame = NS['key_frame']; lev = NS['lev']; tperp = NS['tperp']; xsteps = NS['xsteps']
FULL = NS['FULL']; ZERO = (); ann = NS['ann']; rref = NS['rref']; yd = NS['yd']; xdata = NS['xdata']
m = h * h

# ---- every charged h-space step, independently collected (all ranks), with the classes it occurs in
hs = {}
def add(A, B, kind):
    if len(B) > len(A): hs.setdefault((A, B), set()).add(kind)
rk = Counter(); corner_rank = Counter(); sig_used = set()
for s in range(R):
    fr = [key_frame(k) for k in walk[s]] + [FULL]
    if s not in sel: fr = [ZERO] + fr
    for A, B in zip(fr, fr[1:]):
        assert len(B) >= len(A)
        if len(B) > len(A): rk[len(B) - len(A)] += 1; add(A, B, 'chain')
    sg = SIG[s] if s in sel else ZERO
    assert len(sg) == (f[s] if s in sel else 0) and h - len(sg) == h - f[s], s
    corner_rank[h - f[s]] += 1; sig_used.add(sg)
ych = 0
for t in range(v):
    q = [ZERO] + list(lev.get(t, [])) + [tperp[t]]
    for A, B in zip(q, q[1:]): add(A, B, 'ydata'); ych += 1
for A, B in xsteps: add(A, B, 'xdata')
cent = [ann(rref([[1 - 3 * (q == c) for q in range(h)]], h), h) for c in range(h)]
for B in cent: add(ZERO, B, 'centre')
kinds = Counter(k for ks in hs.values() for k in ks)
log('distinct h-steps %d, by class (a step may occur in several) %s; aux sigmas %d; slots %d' % (len(hs), dict(kinds), len(sig_used), R))
log('rank distribution of h-steps %s' % sorted(Counter(len(B) - len(A) for A, B in hs).items()))
log('aux corner rank r_u = h - dim sigma distribution %s' % sorted(corner_rank.items()))
assert min(corner_rank) >= 1, 'r_u = 0 would give a width-m child'

# ---- O1 exactly over Z: A <= B (annihilator verified: Z B^T = 0 and rank Z = h - dim B mod a prime)
P = (1 << 61) - 1
def rank_mod(M):
    M = [[x % P for x in r] for r in M]; rk_ = 0; n = len(M)
    for c in range(h):
        i = next((i for i in range(rk_, n) if M[i][c]), None)
        if i is None: continue
        M[rk_], M[i] = M[i], M[rk_]; iv = pow(M[rk_][c], P - 2, P)
        for k in range(rk_ + 1, n):
            if M[k][c]:
                fct = M[k][c] * iv % P; M[k] = [(x - fct * y) % P for x, y in zip(M[k], M[rk_])]
        rk_ += 1
    return rk_
frames = {}
def fid(B):
    if B not in frames: frames[B] = len(frames)
    return frames[B]
for (A, B) in hs: fid(A); fid(B)
for sg in sig_used: fid(sg)
fid(FULL)
Zc = {}; bad_ann = 0; bad_rank = 0
for B in frames:
    if len(B) in (0, h): continue
    assert rank_mod(B) == len(B)            # basis rows independent (mod P => over Q)
    Z = ann(B, h)
    if any(sum(a * b for a, b in zip(z, r)) for z in Z for r in B): bad_ann += 1
    if len(Z) != h - len(B) or rank_mod(Z) != h - len(B): bad_rank += 1
    Zc[B] = Z
bad_nest = 0
for (A, B) in hs:
    if len(B) == h or not A: continue
    Z = Zc[B]
    if any(sum(a * z for a, z in zip(ra, rz)) for ra in A for rz in Z): bad_nest += 1
log('O1 exact: frames %d, annihilator failures %d, rank failures %d, nesting failures %d' % (len(frames), bad_ann, bad_rank, bad_nest))
assert bad_ann == bad_rank == bad_nest == 0
mx = max(abs(x) for B in frames for r in B for x in r)
log('max |entry| in frame bases %d' % mx)

# ---- one-run histogram rebuilt from these exact objects
H = Counter()
for _ in range(2):
    for r, c in corner_rank.items(): H[m - r] += v * c
    for r, c in rk.items(): H[r] += v * c
    H[h - 1] += v * h
    for st in yd:
        for x in st: H[x] += v
    for xs in xdata:
        for x in xs: H[x] += v
N = v * v
H[m - 2 * h + 1] += 2 * N; H[1] += N
W = 2 * N + 2 * v * R; L = 2 * v * h * (h - 1); s = W * m - N + L
assert sum(w * c for w, c in H.items()) == s
pickle.dump(dict(hist=dict(H), W=W, s=s, m=m, N=N, L=L, h=h, v=v, R=R, rk=dict(rk), corner_rank=dict(corner_rank)),
            open(os.path.join(OUT, 'hist_gate.pkl'), 'wb'))
log('one-run histogram: W=%d s=%d D=%d widths=%d max=%d' % (W, s, W * m - s, len(H), max(H)))

# ---- binary export for the C checker (entries are small integers; written mod P)
ids = sorted(frames.items(), key=lambda kv: kv[1])
with open(os.path.join(OUT, 'gate_data.bin'), 'wb') as out:
    out.write(struct.pack('<5q', h, v, len(frames), len(hs), len(sig_used)))
    for t in trip: out.write(struct.pack('<3q', *t))
    for B, i in ids:
        out.write(struct.pack('<q', len(B)))
        for r in B: out.write(struct.pack('<%dq' % h, *[x % P for x in r]))
    KC = {'chain': 1, 'ydata': 2, 'xdata': 4, 'centre': 8}
    for (A, B), ks in sorted(hs.items(), key=lambda kv: (frames[kv[0][0]], frames[kv[0][1]])):
        out.write(struct.pack('<4q', frames[A], frames[B], len(B) - len(A), sum(KC[k] for k in ks)))
    for sg in sorted(sig_used, key=lambda b: frames[b]): out.write(struct.pack('<q', frames[sg]))
log('wrote gate_data.bin')
