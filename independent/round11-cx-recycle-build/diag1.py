import sys; sys.argv=['x','11']
from cbirth import *
from collections import Counter, defaultdict
W = load_word(11); F = gfirst_frames(W); tau = late_tau(W, F)
h,m=W['h'],W['m']; R=W['pr']['R']; rops,frames,pos=F['rops'],F['frames'],F['pos']; gsig=F['gsig']
elig=[sr for sr in range(R) if rops[sr] and sr not in F['rootframe'] and sr not in gsig]
print('eligible donors', len(elig), 'roots', len(F['rootframe']), 'gauged', len(gsig), 'cut', F['cut'])
byE=defaultdict(list)
for sr in elig: byE[frames[rops[sr][-1]]].append(pos[rops[sr][-1]])
print('distinct donor last frames', len(byE), 'dims', Counter(len(E) for E in byE).most_common(12))
nc=Counter(); ncT=Counter(); tauc=Counter()
for z in W['s']['sel']:
    A=z['A']; B=z['role']
    comp=[E for E in byE if all(bin(a&f).count('1')%2==0 for a in A for f in E)]
    n1=sum(len(byE[E]) for E in comp); n2=sum(sum(1 for t in byE[E] if t<tau[B]) for E in comp)
    nc[min(n1,5)]+=1; ncT[min(n2,5)]+=1; tauc['phase1' if tau[B]<=F['cut'] else 'rest']+=1
print('compat donors (cap5):', sorted(nc.items()), ' with timing:', sorted(ncT.items()), tauc)
# death times of donors vs deadlines
import statistics
print('donor deaths in phase1:', sum(1 for sr in elig if pos[rops[sr][-1]]<F['cut']))
