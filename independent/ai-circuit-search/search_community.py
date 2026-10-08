"""Apply new sum trees to the community generic first axis; preserve fixed25.

Exploratory, conditional profiles. Uses attributed community graph, matching,
and positive-label producers, including two compiled C++17 matching helpers.
Does not re-prove simultaneous basis, physical transfer, or analytic lemmas.
"""
import argparse
import gc
import json
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'/'croc'/'scripts'))
from partial_swap import graph, paired, positive

BASE=paired.PairedExclusionCircuit
WORK=ROOT/'build'


def mutated(split,cutoff,vector):
    class Modified(BASE):
        base_threshold=cutoff
        def total(self,vals):
            if not vals:return 0
            if len(vals)==1:return vals[0]
            n=len(vals)
            cut={'balanced':n//2,'power':1<<((n-1).bit_length()-1),
                 'left':1,'right':n-1,'third':max(1,n//3)}[split]
            return self.add(self.total(vals[:cut]),self.total(vals[cut:]))
        def vector(self,vals,two=True):
            if two or vector=='prefix':return super().vector(vals,two)
            # Tree all-exclusions: propagate the sum of disjoint siblings.
            if not vals:return 0,[],{}
            tree=[]
            def build(lo,hi):
                if hi-lo==1:
                    tree.append((vals[lo],None,None,lo));return len(tree)-1
                n=hi-lo
                mid=lo+(n//2 if vector=='balanced' else 1<<((n-1).bit_length()-1))
                a,b=build(lo,mid),build(mid,hi)
                tree.append((self.add(tree[a][0],tree[b][0]),a,b,None));return len(tree)-1
            top=build(0,len(vals));one=[0]*len(vals)
            def visit(nd,outside):
                _,a,b,index=tree[nd]
                if a is None:one[index]=outside;return
                visit(a,self.add(outside,tree[b][0]));visit(b,self.add(outside,tree[a][0]))
            visit(top,0)
            return tree[top][0],one,{}
    return Modified


def assemble(row,baseline):
    old=baseline['bit']['counts'];h=23;copies=2300
    hist=Counter({int(k):v for k,v in old['child_multiplicities'].items()})
    for w,n in old['parts']['internal_23'].items():hist[int(w)]-=n
    for w,n in old['parts']['exterior_23'].items():hist[int(w)]-=n
    local=list(row['histogram']);local[h]-=h;local[1]+=h
    assert all(n>=0 for n in local)
    assert sum(r*n for r,n in enumerate(local))==h*row['R']+row['loss']
    for r,n in enumerate(local):
        n*=copies
        if 2*r>h:hist[2*r-h]+=n;hist[1]+=(h-r)*n
        else:hist[1]+=r*n
    bank=copies*row['R'];hist[h]+=bank;hist[old['m']-2*h]+=bank
    W=old['W']-old['B1']+bank
    rank=W*old['m']-old['N']+old['L']
    assert all(n>=0 for n in hist.values())
    assert sum(w*n for w,n in hist.items())==rank
    net=dict(m=old['m'],W=W,hist=dict(sorted((w,n) for w,n in hist.items() if n)),rank_sum=rank)
    items=[(w*n/(net['m']*W),math.log(net['m']/w)) for w,n in net['hist'].items()]
    deficit=(net['m']*W-rank)/(net['m']*W)
    lo,hi=0.,.001
    for _ in range(70):
        mid=(lo+hi)/2
        if sum(weight*math.expm1(mid*l) for weight,l in items)<deficit:lo=mid
        else:hi=mid
    return net,(lo+hi)/2


def main():
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=45);args=p.parse_args()
    WORK.mkdir(exist_ok=True)
    baseline=json.loads((ROOT/'croc-certificate.json').read_text())
    expected=next(r for r in json.loads((ROOT/'vendor'/'croc'/'certificates'/'copied-centers-bit-axes.json').read_text()) if r['h']==23)
    configs=[(s,c,v) for v in ['prefix','balanced','power'] for c in [2,3,4] for s in ['balanced','power','left','right','third']]
    data=dict(search=[],best=None,scope='Experimental community profiles; only baseline cross-checked, changed transfer remains conditional')
    for i,(split,cutoff,vector) in enumerate(configs[:args.limit]):
        graph.PairedExclusionCircuit=mutated(split,cutoff,vector)
        try:
            g=graph.graph(23);scalar=g.verify();dag=WORK/'current.bin';graph.export(g,dag)
            del g;gc.collect()
            original=json.loads(subprocess.check_output([str(WORK/'match_exported_dag.exe'),str(dag),str(dag)+'.links'],stderr=subprocess.DEVNULL,text=True))
            labels=positive.run(str(dag))
            matched=json.loads(subprocess.check_output([str(WORK/'match_positive_dag.exe'),str(dag),str(dag)+'.positive'],stderr=subprocess.DEVNULL,text=True))
            assert original['matched']==matched['matched']
            assert matched['loss']==23*22
            net,saving=assemble(matched,baseline)
            row=dict(split=split,cutoff=cutoff,vector=vector,saving=saving,producer=matched,scalar=scalar,labels=labels)
            if i==0:
                for key in ('h','v','c','q','R','loss','histogram'):assert matched[key]==expected[key],key
                assert matched['matched']==expected['matching']
                b=baseline['bit']['counts']
                assert net==dict(m=b['m'],W=b['W'],hist={int(w):n for w,n in b['child_multiplicities'].items()},rank_sum=b['total_rank'])
                data['baseline_reproduced']=True
            if data['best'] is None or saving>data['best']['saving']:
                data['best']={**row,'network':net};print('BEST',i,split,cutoff,vector,saving,flush=True)
        except (AssertionError,subprocess.CalledProcessError) as e:
            row=dict(split=split,cutoff=cutoff,vector=vector,rejected=type(e).__name__,reason=str(e))
        data['search'].append(row)
        (ROOT/'community-search-results.json').write_text(json.dumps(data,indent=2))
        print('progress',i+1,len(configs[:args.limit]),'saving',row.get('saving','rejected'),flush=True)


if __name__=='__main__':main()
