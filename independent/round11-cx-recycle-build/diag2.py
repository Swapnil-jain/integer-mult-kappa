import sys
from resel_br import *
import resel_br
W = load_word(11); F0 = gfirst_frames(W); h,m=22,66; a0=5.6e-4
# instrument: replicate the candidate loop quickly, record stats for rejected
G = W['wit']['g']; pr, s = W['pr'], W['s']; v = pr['v']; R = pr['R']; roots=G['roots']; inputs=G['inputs']; ann=W['wit']['ann']
ops=s['ops']; first=s['first']; touched=s['touched']
co=[0]*R
for r, sr in zip(roots, s['rootroles']): co[sr] |= sum(1 << t for t in r['targets'])
for a, b, x in reversed(ops): co[b] |= co[a]
cands = sorted((x for x in range(R) if x not in touched and first[x] is not None), key=lambda x: (len(ann[first[x]]), bin(co[x]).count("1"), x))
gs={z['role'] for z in s['sel']}
st=Counter(); ntg=Counter(); rootc=Counter()
rr=set(s['rootroles'])
for x in cands:
    if x in gs: continue
    r=h-len(ann[first[x]]); st[(r, bin(co[x]).count('1'))]+=1; rootc[x in rr]+=1
print('rejected (r, #co targets):', st.most_common(15)); print('root roles among rejected', rootc)
print('gauged (r, targets)', Counter((z['r'], len(z['targets'])) for z in s['sel']).most_common(5))
limit = [None] * v
for r in roots:
    if r['kind'] == 'center': continue
    A = basis(inputs[t] for t in r['targets'])
    for t in r['targets']:
        if limit[t] is None: limit[t] = A
for t in range(v):
    if limit[t] is None: limit[t] = (inputs[t],)
selA={z['role']:z for z in s['sel']}
rows=Counter(); full=0
for x in cands:
    A = tuple(ann[first[x]]); tg = []; bits = co[x]
    while bits:
        low = bits & -bits; t = low.bit_length() - 1; bits ^= low; tg.append(t)
        A = basis(A + tuple(limit[t]))
        if len(A) == h: break
    if len(A)==h: full+=1; continue
    r = h - len(ann[first[x]]); d = h - len(A)
    base = 3 * (ex(r - d, a0, m) - ex(r, a0, m))
    tgt = 3 * sum(ex(d, a0, m) + ex(h - len(limit[t]) - d, a0, m) - ex(h - len(limit[t]), a0, m) for t in tg)
    bestpair = 3 * (ex(0, a0, m) - ex(h - d, a0, m))   # ideal donor e = d
    if x in selA:
        for t in tg: limit[t] = A
    else: rows[(r,d,len(tg), round((base+tgt+bestpair)/a0,1))]+=1
print('full A (skipped):', full)
print('rejected (r,d,#tg, ideal-paired delta/a0):', rows.most_common(20))
Es, byE = compat_table(F0, R)
pos=F0['pos']; rops=F0['rops']
limit2 = [None] * v
for r in roots:
    if r['kind'] == 'center': continue
    A = basis(inputs[t] for t in r['targets'])
    for t in r['targets']:
        if limit2[t] is None: limit2[t] = A
for t in range(v):
    if limit2[t] is None: limit2[t] = (inputs[t],)
limit=limit2
res=Counter(); shown=0
for x in cands:
    A = tuple(ann[first[x]]); tg = []; bits = co[x]
    while bits:
        low = bits & -bits; t = low.bit_length() - 1; bits ^= low; tg.append(t)
        A = basis(A + tuple(limit[t]))
        if len(A) == h: break
    r = h - len(ann[first[x]]); d = h - len(A)
    if x in selA:
        for t in tg: limit[t] = A
        continue
    if (r,d) not in ((17,16),(3,2)): continue
    dl=pos[rops[x][0]]
    es=Counter(); est=Counter()
    for E in Es:
        if all(bin(a&f).count('1')%2==0 for a in A for f in E):
            es[len(E)]+=len(byE[E]); est[len(E)]+=sum(1 for t,_ in byE[E] if t<dl)
    res[(r,d,max(es) if es else None, max(est) if est else None)]+=1
    if shown<3: print(x, r,d,'deadline',dl,'cut',F0['cut'],'compat by e',dict(es),'timely',dict(est)); shown+=1
print(res)
