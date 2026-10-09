"""Round-seven witness 2, part A (stdlib, exact over Q): every frame of the lifted + late-copy program.

From the definitions and the frozen program (certificates/round7/jfstack_23.json.gz):
  span(g)  = span{t_S : S in sup(g)} for every region representative g of the compiled program;
  K_g      = span of the root covectors reachable from g along the roles ((9 t_T - 3) for a role ending in an
             output to Y_T, 1 - 3 e_c for a retained role), U_g = ann K_g;
  M_g      = span(g) + (U_g cap rowspace F_j(g)), j(g) from the program, F = the frozen flag; leaves <t_S>;
  roots    t_T^perp, U_c;  late-copy frames C = exact intersections of their recorded component frames.
Checks: j monotone along every role; K_g kills span(g) over Z; regions sharing a frame id have the same exact frame;
every recomputed frame equals the frozen exact basis (jfstack_frames_23.json.gz, used by parts B and C) and its
recorded dimension; every consecutive pair of every role chain nested exactly (equal dims => equal); G = I - J/9
nondegenerate on every frame used (exact integer Gram det != 0 mod a prime). Negative controls at the end.
Usage: python3 check_frames.py"""
import time, sys
from collections import defaultdict, Counter
from xq import rref, ann, inside_ann, meet_ann, red, rref_mod, nondeg
T0 = time.time()
def log(*a): print('[%5.0fs]' % (time.time() - T0), *a, flush=True)
OK = {}
def report(name, ok):
    OK[name] = bool(ok); log('%-100s %s' % (name, 'PASS' if ok else 'FAIL'))
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jfdata, random
KD = jfdata.load(); EXF, _ = jfdata.load_frames()
D = KD['D']; X = KD['X']; meta = KD['meta']
h = 23; G = D['G']; C = D['C']; trip = [tuple(t) for t in G['trip']]; v = len(trip); args = G['args']
hold = C['hold']; j = X['j']; Mid = X['Mid']
tv = lambda T: [int(q in T) for q in range(h)]
Fl = KD['Fl']
report('0. flag F invertible over Q (full rank, exact)', len(rref(Fl, h)) == h)
annF = [ann(rref(Fl[:x], h), h) for x in range(h + 1)]

# ---------------------------------------------------------------- reachability and K_g
reps = sorted({n for hs in hold for n in hs} | set(meta['opreg']), key=lambda n: (G['dim'][n], n))
succ = defaultdict(set); rootcov = defaultdict(list); badj = 0
for s, hs in enumerate(hold):
    for a, b in zip(hs, hs[1:]): succ[a].add(b); badj += j.get(a, 0) > j.get(b, 0)
    if s in C['out']: rootcov[hs[-1]].append(tuple(9 * x - 3 for x in tv(C['out'][s][1])))
    if s in C['ret']: rootcov[hs[-1]].append(tuple(1 - 3 * (q == C['ret'][s]) for q in range(h)))
report('3. j monotone along every jf role (%d roles, %d successor pairs)' % (len(hold), sum(len(x) for x in succ.values())), badj == 0)
order = sorted(reps, key=lambda n: (G['dim'][n], n))
pos = {n: i for i, n in enumerate(order)}
assert all(pos[a] < pos[b] for a in succ for b in succ[a]), 'successor not later in region order'
Kid = {}; Ktab = []; Kix = {}; ucache = {}
def kid(B):
    if B not in Kix: Kix[B] = len(Ktab); Ktab.append(B)
    return Kix[B]
for n in reversed(order):
    kids = tuple(sorted({Kid[b] for b in succ.get(n, ())})); rc = tuple(sorted(set(rootcov.get(n, ()))))
    key = (kids, rc)
    if key not in ucache:
        if not rc and len(kids) == 1: ucache[key] = kids[0]
        else:
            rows = list(rc)
            for k in kids: rows.extend(Ktab[k])
            ucache[key] = kid(rref(rows, h))
    Kid[n] = ucache[key]
log('K_g: %d distinct covector spaces over %d regions' % (len(Ktab), len(order)))
spanB = {}; badkill = 0; baddim = 0
for n in order:
    B = rref([tv(trip[i]) for i in range(v) if G['sup'][n] >> i & 1], h); spanB[n] = B
    baddim += len(B) != G['dim'][n]
    K = Ktab[Kid[n]]
    badkill += any(sum(a * b for a, b in zip(kr, sr)) for kr in K for sr in B)
report('3. exact dim span(g) == recorded dim for all %d regions' % len(order), baddim == 0)
report('3. every reachable root covector kills span(g) over Z (exact dot products)', badkill == 0)
log('spans done')

# ---------------------------------------------------------------- exact region frames
Mex = {}; Uann = {}
for n in order:
    if args[n] is None or not j.get(n, 0): Mex[n] = spanB[n]; continue
    kk = (Kid[n], j[n])
    if kk not in Uann: Uann[kk] = ann(rref(list(annF[j[n]]) + list(Ktab[Kid[n]]), h), h)   # U_g cap rowspace F_j
    Mex[n] = rref(list(spanB[n]) + list(Uann[kk]), h)
del Uann
log('exact region frames: %d (lifted %d)' % (len(Mex), sum(1 for n in order if args[n] is not None and j.get(n, 0))))
EX = {}; bad_share = 0; src = {}
for n in order:
    i = Mid[n]
    if i in EX: bad_share += EX[i] != Mex[n]
    else: EX[i] = Mex[n]; src[i] = ('reg', n)
report('3. regions sharing a frame id have the same exact frame', bad_share == 0)
del Mex
for s, (c, T) in C['out'].items():
    i = X['rootid']['out'][X['phys'][s]]; B = ann(rref([[9 * x - 3 for x in tv(T)]], h), h)
    if i in EX: bad_share += EX[i] != B
    else: EX[i] = B; src[i] = ('out', T)
for s, c in C['ret'].items():
    i = X['rootid']['ret'][X['phys'][s]]; B = ann(rref([[1 - 3 * (q == c) for q in range(h)]], h), h)
    if i in EX: bad_share += EX[i] != B
    else: EX[i] = B; src[i] = ('ret', c)
for u, d in enumerate(meta['usedest']):        # use frames D(u) (late-copy components)
    i = X['dfr'][u]
    if d[0] == 'reg': B = EX[Mid[d[1]]]
    elif d[0] == 'out': B = ann(rref([[9 * x - 3 for x in tv(d[1][1])]], h), h)
    else: B = ann(rref([[1 - 3 * (q == d[1]) for q in range(h)]], h), h)
    if i in EX: bad_share += EX[i] != B
    else: EX[i] = B; src[i] = ('use', u)
report('3. root and use frames sharing ids agree exactly', bad_share == 0)
ANN = {}
def Z(i):
    if i not in ANN: ANN[i] = ann(EX[i], h)
    return ANN[i]
ncap = 0
for i, ids in sorted(X['capdef'].items()):
    assert all(k in EX for k in ids), 'cap component without exact frame'
    Zs = []
    for k in ids: Zs.extend(Z(k))
    B = ann(rref(Zs, h), h)
    if i in EX: bad_share += EX[i] != B
    else: EX[i] = B; src[i] = ('cap', ids); ncap += 1
log('late-copy intersections exact: %d' % ncap)
used = sorted({i for ch in X['chains'] for i in ch} | set(X['opfr']))
missing = [i for i in used if i not in EX]
report('2. every frame id used by a chain or a gate has an exact definition (%d used, %d missing)' % (len(used), len(missing)), not missing)
report('2. recomputed exact frame == frozen exact basis, and recorded dim (every used id)',
       sorted(EXF) == used and all(EX[i] == EXF[i] and len(EX[i]) == KD['dims'][i] for i in used))
nd = [i for i in used if not nondeg(EX[i])]
nd2 = [i for i in nd if not nondeg(EX[i], 2**61 - 1)]
report('2. G nondegenerate on every used frame (%d frames; Gram det != 0 mod a prime => over Z)' % len(used), not nd2)
badn = badeq = 0
for ch in X['chains']:
    for a, b in zip(ch, ch[1:]):
        if len(EX[a]) > len(EX[b]) or not inside_ann(EX[a], Z(b)): badn += 1
        elif len(EX[a]) == len(EX[b]) and EX[a] != EX[b]: badeq += 1
report('3. every consecutive pair of every final role chain nested exactly (%d roles)' % len(X['chains']), badn == 0 and badeq == 0)
bits = Counter(max((abs(x).bit_length() for r in EX[i] for x in r), default=0) // 16 * 16 for i in used)
log('   entry sizes (bits, bucketed): %s; distinct exact subspaces among used ids: %d' % (sorted(bits.items()), len(set(EX[i] for i in used))))

# ---------------------------------------------------------------- negative controls
line = rref([[1] * 9 + [0] * 14], h)
report('N. NEG the G test rejects the G-degenerate line e0+..+e8', not nondeg(line) and not nondeg(line, 2**61 - 1))
rng = random.Random(11); tried = det = 0
for s in rng.sample(range(len(X['chains'])), 3000):
    ch = X['chains'][s]; ks = [k for k in range(len(ch) - 1) if len(EX[ch[k]]) < len(EX[ch[k + 1]])]
    if not ks: continue
    k = rng.choice(ks); tried += 1; det += not inside_ann(EX[ch[k + 1]], ann(EX[ch[k]], h))
report('N. NEG swapped consecutive chain frames detected (%d / %d)' % (det, tried), det == tried and tried)
tried = brk = 0; roots = {**X['rootid']['out'], **X['rootid']['ret']}
for p, ch in rng.sample(list(enumerate(X['chains'])), 4000):
    if p not in roots: continue
    for i in ch[:-1]:
        s_ = src.get(i)
        if s_ and s_[0] == 'reg' and j.get(s_[1], 0) >= 3:
            B = rref(list(EX[i]) + Fl[:j[s_[1]]], h); tried += 1
            brk += not inside_ann(B, ann(EX[roots[p]], h)); break
    if tried >= 300: break
report('N. NEG lift ignoring U_g (span + F_j) leaves the role root frame (%d / %d)' % (brk, tried), brk > 0)
log('ALL %s' % ('PASS' if all(OK.values()) else 'FAIL: ' + str([k for k, x in OK.items() if not x])))
sys.exit(0 if all(OK.values()) else 1)
