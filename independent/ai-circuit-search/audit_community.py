"""Independent exact scalar, subspace and matched-transition audit of winner.

Labels are interpreted as actual integer bases in Q^23. All basis vectors
satisfy sum(x)=3*x_c for a common centre, proving positivity for I-J/9.
No floating point or modular-zero inference occurs in this audit.
Inherited physical basis/stream/tape/analytic transfer is still conditional.
"""
import json
import gzip
import argparse
import struct
import subprocess
import sys
from collections import Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import search_community as sc
from partial_swap import graph, positive
from partial_swap.binary import read_array

ROOT=Path(__file__).resolve().parent
OUT=ROOT
sys.path.insert(0,str(OUT))
from verify import require


def audit(g,dag,matchfile,claimed):
    h=g.h;n=len(g.args);trip=list(combinations(range(h),3))
    require(g.inputs==trip,'input enumeration')
    # Reconstruct every scalar sum independently, with dense global supports.
    supports=[0]*n;degree=[0]*n
    for nd in sorted(g.active):
        if g.args[nd]:
            a,b=g.args[nd]
            require(a<nd and b<nd,'DAG order')
            require(not supports[a]&supports[b],'overlapping scalar sums')
            supports[nd]=supports[a]|supports[b]
            degree[a]+=1;degree[b]+=1
        else:
            require(1<=nd<=len(trip),'source identifier')
            supports[nd]=1<<(nd-1)
    stars=[sum(1<<j for j,t in enumerate(trip) if c in t) for c in range(h)]
    roots=list(g.outputs.values());targets=list(g.outputs)
    for (c,t),nd in g.outputs.items():
        expected=stars[c]
        for v in t:
            if v!=c:expected &= ~stars[v]
        require(supports[nd]==expected,'wrong exact output sum')
        degree[nd]+=1
    label=Path(str(dag)+'.positive')
    label=label if label.exists() else label.with_name(label.name+'.gz')
    with (gzip.open(label,'rb') if label.suffix=='.gz' else open(label,'rb')) as stream:
        require(struct.unpack('<2I',stream.read(8))==(h,n),'label header')
        ranks=read_array(stream,'I',n);forced=read_array(stream,'Q',n)
        sy=read_array(stream,'b',n*h)
        require(stream.read()==b'','label trailing data')
    labels={nd:(forced[nd],tuple(sy[nd*h:(nd+1)*h])) for nd in g.active}

    @lru_cache(None)
    def basis(label):
        F,s=label;fixed=[i for i in range(h) if F>>i&1];f=len(fixed)
        require(f in (1,2,3),'label common centre count')
        require(all((s[i]==1)==bool(F>>i&1) for i in range(h)),'fixed symbols')
        require(all(z==0 or z==1 or abs(z)>=2 for z in s),'symbol range')
        if f==3:
            require(all(z in (0,1) for z in s),'source symbols')
            return (tuple(int(i in fixed) for i in range(h)),)
        components=sorted({abs(z) for z in s if abs(z)>=2})
        vectors=[]
        for comp in components:
            signs=[(1 if z>0 else -1) if abs(z)==comp else 0 for z in s]
            mass=sum(signs)
            vectors.append(tuple(mass if i in fixed else (3-f)*signs[i] for i in range(h)))
        require(bool(vectors),'empty active label')
        # Disjoint nonzero component coordinates prove linear independence.
        for v in vectors:
            require(all(sum(v)==3*v[c] for c in fixed),'common-centre positivity identity')
        return tuple(vectors)

    def contains(label,x):
        F,s=label;fixed=[i for i in range(h) if F>>i&1]
        if not fixed:return False
        if any(sum(x)!=3*x[c] for c in fixed):return False
        values={}
        for i,z in enumerate(s):
            if z==0:
                if x[i]:return False
            elif abs(z)>=2:
                value=x[i]*(1 if z>0 else -1);key=abs(z)
                if key in values and values[key]!=value:return False
                values[key]=value
        return True

    @lru_cache(None)
    def included(a,b):return all(contains(b,x) for x in basis(a))

    oldrank={nd:(g.union[nd].bit_count()-g.core[nd].bit_count()) if g.args[nd] else 1 for nd in g.active}
    nodekey=lambda nd:(ranks[nd],oldrank[nd],nd)
    inclusions=0
    for nd in sorted(g.active):
        require(len(basis(labels[nd]))==ranks[nd],'label dimension')
        if g.args[nd]:
            for operand in g.args[nd]:
                require(included(labels[operand],labels[nd]),'gate subspace inclusion')
                require(nodekey(operand)<nodekey(nd),'gate partial order')
                inclusions+=1
        else:
            v=tuple(int(i in trip[nd-1]) for i in range(h))
            require(contains(labels[nd],v),'input vector outside label')
    for (c,t),nd in g.outputs.items():
        if len(t)==1:
            require(ranks[nd]==h-1 and forced[nd]>>c&1,'retained centre space')
        else:
            for v in basis(labels[nd]):
                require(3*sum(v[i] for i in t)==sum(v),'output not orthogonal to target')
    hist=Counter();adds=0;loss=0
    for nd in sorted(g.active):
        r=ranks[nd];require(degree[nd]>0,'unused active node')
        if g.args[nd]:
            adds+=1;hist[r]+=degree[nd]-1;hist[h-r]+=1
            for operand in g.args[nd]:hist[r-ranks[operand]]+=1
        else:hist[1]+=degree[nd]
    for (_,t),nd in g.outputs.items():
        r=ranks[nd]
        if len(t)==1:hist[r]+=1;hist[h]+=1;loss+=r
        else:hist[h-1-r]+=1;hist[1]+=1
    with (gzip.open(matchfile,'rb') if str(matchfile).endswith('.gz') else open(matchfile,'rb')) as stream:
        nn,nm=struct.unpack('<2I',stream.read(8));require(nn==n,'matching header')
        edges=[struct.unpack('<2I',stream.read(8)) for _ in range(nm)]
        require(stream.read()==b'','matching trailing data')
    donors=set();uses=set()
    for donor,encoded in edges:
        require(donor not in donors and encoded not in uses,'not a matching')
        donors.add(donor);uses.add(encoded)
        if encoded>>31:
            j=encoded&0x7fffffff;target=roots[j];value=target;order=n+j
        else:
            target=encoded//2;value=g.args[target][encoded&1];order=target
        require(value in g.args[donor],'unrelated carrier')
        require(included(labels[donor],labels[target]),'matched subspace inclusion')
        require(nodekey(donor)<(ranks[target],oldrank[target],order),'cyclic matched order')
        ru,rv,rt=ranks[donor],ranks[value],ranks[target]
        require(rv<=ru<=rt,'matched dimensions decrease')
        hist[h-ru]-=1;hist[rv]-=1;hist[rt-rv]-=1;hist[rt-ru]+=1
    roles=adds+len(roots)-nm
    require(all(0<=r<=h and count>=0 for r,count in hist.items()),'negative transition count')
    require(sum(r*count for r,count in hist.items())==h*roles+2*loss,'rank mass')
    require([hist[r] for r in range(h+1)]==claimed['histogram'],'C++ histogram mismatch')
    for k,v in dict(c=adds,q=len(roots),R=roles,matched=nm,loss=loss).items():require(claimed[k]==v,k+' count mismatch')
    # An actual source vector outside a target label must be rejected.
    bad=next(nd for nd in g.active if not g.args[nd]);corrupt=[0]*h;corrupt[0]=1
    require(not contains(labels[bad],corrupt),'invalid-vector negative control')
    return dict(status='PASS',active_nodes=len(g.active),unique_subspaces=len(set(labels.values())),
        gate_inclusions=inclusions,matching_inclusions=nm,output_interfaces=len(roots),roles=roles,
        rank_sum=sum(r*count for r,count in hist.items()),all_scalar_outputs_exact=True,
        integer_bases_independent=True,positive_form_identity_exact=True,all_inclusions_exact=True,
        augmented_dependency_order_acyclic=True,histogram_recomputed=True,invalid_vector_rejected=True,
        scope='Finite circuit/subspace/matching/transition audit; physical simultaneous-basis, tape, analytic and exact-recovery transfer remain inherited hypotheses')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-saved',action='store_true',help='Audit supplied labels/matching with Python only; do not rerun C++ producers')
    options=parser.parse_args()
    best=json.loads((ROOT/'community-search-results.json').read_text())['best']
    graph.PairedExclusionCircuit=sc.mutated(best['split'],best['cutoff'],best['vector'])
    g=graph.graph(23);g.verify()
    work=ROOT/'witnesses' if options.verify_saved else sc.WORK
    work.mkdir(exist_ok=True)
    dag=work/'winner.bin'
    matchfile=work/('winner.matching.gz' if options.verify_saved else 'winner.matching')
    if options.verify_saved:
        matched=best['producer']
    else:
        graph.export(g,dag)
        original=json.loads(subprocess.check_output([str(sc.WORK/'match_exported_dag.exe'),str(dag),str(dag)+'.links'],stderr=subprocess.DEVNULL,text=True))
        positive.run(str(dag))
        matched=json.loads(subprocess.check_output([str(sc.WORK/'match_positive_export.exe'),str(dag),str(dag)+'.positive',str(matchfile)],stderr=subprocess.DEVNULL,text=True))
        require(matched==best['producer'],'winner reproduction')
        require(matched['matched']==original['matched'],'matching enlargement')
    log=audit(g,dag,matchfile,matched)
    baseline=json.loads((ROOT/'croc-certificate.json').read_text())
    generated,_=sc.assemble(matched,baseline)
    certificate=json.loads((OUT/'certificate.json').read_text())
    expected=certificate['networks']['community_candidate_bit']
    expected['hist']={int(k):v for k,v in expected['hist'].items()}
    require(generated==expected,'full candidate profile mismatch')
    log['full_child_histogram_reconstructed']=True
    print(json.dumps(log,indent=2))


if __name__=='__main__':main()
