"""Broader positive-circuit search: grouping arity, exclusion circuits and tree layout.

All sums have disjoint supports. These changes preserve the exact scalar outputs;
their cost and frame implications are still conditional on the upstream lemmas.
"""
import itertools
import json
import math
import random
import sys
from pathlib import Path
import search
import sidegen, rtgm

BASE=sidegen.Excl


def make_class(group_size=2,cutoff=4,split='power',loo='prefix',join='original'):
    class Experimental(BASE):
        def total(self,vals):
            vals=[v for v in vals if v]
            if not vals:return 0
            if len(vals)==1:return vals[0]
            if split=='left':
                out=0
                for v in vals:out=self.add(out,v)
                return out
            if split=='right':
                out=0
                for v in reversed(vals):out=self.add(v,out)
                return out
            if split=='small':vals=sorted(vals,key=lambda n:(self.support[n].bit_count(),n))
            if split=='large':vals=sorted(vals,key=lambda n:(-self.support[n].bit_count(),n))
            cut=1<<((len(vals)-1).bit_length()-1) if split=='power' else len(vals)//2
            return self.add(self.total(vals[:cut]),self.total(vals[cut:]))

        def loo(self,vals):
            if loo=='prefix':return super().loo(vals)
            if not vals:return 0,[]
            if len(vals)==1:return vals[0],[0]
            if loo=='direct':return self.total(vals),[self.total(vals[:i]+vals[i+1:]) for i in range(len(vals))]
            cut=(1<<((len(vals)-1).bit_length()-1)) if loo=='power' else len(vals)//2
            lt,le=self.loo(vals[:cut]);rt,re=self.loo(vals[cut:])
            return self.add(lt,rt),[self.add(x,rt) for x in le]+[self.add(lt,x) for x in re]

        def block(self,pts,edges,wts):
            def e(a,b):return edges[tuple(sorted((a,b)))]
            def inside(keep):return [wts[a] for a in keep]+[e(a,b) for a,b in itertools.combinations(keep,2)]
            if len(pts)<=cutoff:
                def omit(om):
                    # Preserve the upstream input order for a direct comparison.
                    return self.total([v for p,v in edges.items() if not set(p)&set(om)]+[v for p,v in wts.items() if p not in om])
                return omit(()),{a:omit((a,)) for a in pts},{(a,b):omit((a,b)) for a,b in itertools.combinations(pts,2)}
            groups=[pts[i:i+group_size] for i in range(0,len(pts),group_size)]
            ng=len(groups)
            coarse={(i,j):self.total([e(a,b) for a in groups[i] for b in groups[j]]) for i,j in itertools.combinations(range(ng),2)}
            wt={i:self.total(inside(g)) for i,g in enumerate(groups)}
            total,outside,far=self.block(list(range(ng)),coarse,wt)
            strips={};sums={}
            for i,g in enumerate(groups):
                others=[j for j in range(ng) if j!=i]
                for a in g:
                    keep=[u for u in g if u!=a]
                    carry=self.total(inside(keep))
                    vals=[self.total([e(u,w) for u in keep for w in groups[j]]) for j in others]
                    st,one=self.loo([carry]+vals)
                    strips[a]=dict(zip(others,one[1:]));sums[a]=st
            single={a:self.add(outside[i],sums[a]) for i,g in enumerate(groups) for a in g}
            out={}
            for i,g in enumerate(groups):
                for a,b in itertools.combinations(g,2):
                    keep=[u for u in g if u not in (a,b)]
                    if not keep:out[a,b]=outside[i]
                    else:
                        external=[e(u,w) for u in keep for j,gg in enumerate(groups) if j!=i for w in gg]
                        out[a,b]=self.add(outside[i],self.total(inside(keep)+external))
            for i,j in itertools.combinations(range(ng),2):
                for a in groups[i]:
                    left=self.add(far[i,j],strips[a][j]) if join=='original' else None
                    for b in groups[j]:
                        cross=self.total([e(u,w) for u in groups[i] if u!=a for w in groups[j] if w!=b])
                        A,B,C,D=far[i,j],strips[a][j],strips[b][i],cross
                        if join=='original':value=self.add(left,self.add(C,D))
                        elif join=='reverse':value=self.add(self.add(A,C),self.add(B,D))
                        elif join=='cross':value=self.add(self.add(A,D),self.add(B,C))
                        elif join=='chain':value=self.add(self.add(self.add(A,B),C),D)
                        else:value=self.total([A,B,C,D])
                        out[a,b]=value
            return total,single,out
    return Experimental


def make_order(name,group_size):
    if hasattr(search.orders,name):return getattr(search.orders,name)
    if name=='xor':return lambda h,c:sorted((j for j in range(h) if j!=c),key=lambda j:c^j)
    if name=='xorrev':return lambda h,c:sorted((j for j in range(h) if j!=c),key=lambda j:c^j,reverse=True)
    if name=='group':
        def grouped(h,c):
            own=c//group_size
            return [j for j in range(h) if j//group_size!=own]+[j for j in range(h) if j//group_size==own and j!=c]
        return grouped
    raise ValueError(name)


def exact_outputs(loc):
    for (a,b),nd in loc.outputs.items():
        want=sum(1<<i for i,p in enumerate(loc.inputs) if a not in p and b not in p)
        if loc.support[nd]!=want:raise ValueError('wrong pair-exclusion support')


def score(cfg):
    cls=make_class(**{k:cfg[k] for k in ['group_size','cutoff','split','loo','join']})
    # Check scalar outputs before constructing any asymptotic certificate.
    loc=cls(cfg['h']-1);exact_outputs(loc)
    sidegen.Excl=cls
    c=search.build(cfg['h'],make_order(cfg['order'],cfg['group_size']))
    net=search.simple(c);m,W,hist=search.verify.normalize_network(net)
    weights=[(w*n/(m*W),math.log(m/w)) for w,n in hist.items()]
    low,high=0.,.001
    for _ in range(60):
        a=(low+high)/2
        # expm1 isolates the tiny improvement from the near-one baseline.
        val=sum(weight*math.expm1(a*log) for weight,log in weights)
        if val < (m*W-c['s'])/(m*W):low=a
        else:high=a
    return dict(config=cfg,saving=(low+high)/2,roles=c['R'],network=net)


def main():
    phase=sys.argv[1] if len(sys.argv)>1 else 'broad'
    result=dict(phase=phase,search=[],best=None)
    configs=[]
    if phase=='broad':
        for group,cutoff,loo,join in itertools.product([2,3,4],[4,6],['prefix','balanced','power'],['original','cross','chain']):
            configs.append(dict(group_size=group,cutoff=cutoff,split='power',loo=loo,join=join,h=23,order='group' if group!=2 else 'gm'))
        for split,order,join in itertools.product(['balanced','power','small','large'],['gm','xor','xorrev'],['original','reverse','cross']):
            configs.append(dict(group_size=2,cutoff=4,split=split,loo='prefix',join=join,h=23,order=order))
    elif phase=='smallbase':
        for cutoff,split,loo,join in itertools.product([2,3],['balanced','power','left','right'],['prefix','balanced','power'],['original','cross']):
            configs.append(dict(group_size=2,cutoff=cutoff,split=split,loo=loo,join=join,h=23,order='gm'))
    else:
        prior=json.loads((search.ROOT/'experimental-broad.json').read_text())
        top=sorted([r for r in prior['search'] if 'saving' in r],key=lambda r:r['saving'],reverse=True)
        seen=set();chosen=[]
        for row in top:
            cfg=row['config'];key=tuple(cfg[k] for k in ['group_size','cutoff','split','loo','join','order'])
            if key not in seen:chosen.append(cfg);seen.add(key)
            if len(chosen)==8:break
        for cfg,h in itertools.product(chosen,range(18,34)):
            configs.append({**cfg,'h':h})
    for i,cfg in enumerate(configs):
        try:
            row=score(cfg)
            if result['best'] is None or row['saving']>result['best']['saving']:
                result['best']=row
                print('BEST',cfg,'saving',row['saving'],'roles',row['roles'],flush=True)
            result['search'].append({k:v for k,v in row.items() if k!='network'})
        except (AssertionError,ValueError,KeyError) as exc:
            result['search'].append(dict(config=cfg,rejected=repr(exc)))
            print('rejected',cfg,repr(exc),flush=True)
        if i%10==9:print('progress',i+1,'/',len(configs),flush=True)
        (search.ROOT/('experimental-'+phase+'.json')).write_text(json.dumps(result,indent=2))
    sidegen.Excl=BASE


if __name__=='__main__':main()
