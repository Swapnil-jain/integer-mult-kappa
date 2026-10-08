"""Explore valid disjoint-sum constructions. All results are conditional candidates."""
import json
import math
from pathlib import Path
import search
import sidegen

BASE = sidegen.Excl


def factory(cutoff, split):
    class Variant(BASE):
        def total(self, vals):
            if not vals: return 0
            if len(vals)==1: return vals[0]
            if split=='left':
                result=0
                for value in vals: result=self.add(result,value)
                return result
            if split=='right':
                result=0
                for value in reversed(vals): result=self.add(value,result)
                return result
            if split=='power':
                mid=1 << ((len(vals)-1).bit_length()-1)
            else:
                mid=len(vals)//2
            return self.add(self.total(vals[:mid]),self.total(vals[mid:]))

        def block(self, pts, edges, wts):
            if len(pts)<=cutoff:
                def total(omit):
                    return self.total([v for p,v in edges.items() if not set(p)&set(omit)]
                        +[v for p,v in wts.items() if p not in omit])
                from itertools import combinations
                return total(()),{a:total((a,)) for a in pts},{(a,b):total((a,b)) for a,b in combinations(pts,2)}
            return super().block(pts,edges,wts)
    return Variant


def main():
    data=[];best=None
    for cutoff in [4,5,6,7,8]:
        for split in ['balanced','power','left','right']:
            sidegen.Excl=factory(cutoff,split)
            for h in [21,22,23,24,25,26,27]:
                try:
                    c=search.build(h,search.orders.gm)
                    net=search.simple(c);m,W,hist=search.verify.normalize_network(net)
                    lo,hi=0.,.001
                    for _ in range(60):
                        a=(lo+hi)/2
                        if sum(n*(w/m)**(1-a) for w,n in hist.items())<W:lo=a
                        else:hi=a
                    v=(lo+hi)/2
                    row=dict(h=h,cutoff=cutoff,split=split,saving=v,roles=c['R'])
                    data.append(row)
                    if best is None or v>best['saving']:
                        best={**row,'network':net}
                        print('BEST',row,flush=True)
                except AssertionError as exc:
                    data.append(dict(h=h,cutoff=cutoff,split=split,rejected=str(exc)))
            print('finished',cutoff,split,flush=True)
            (search.ROOT/'variant-results.json').write_text(json.dumps(dict(search=data,best=best),indent=2))
    sidegen.Excl=BASE


if __name__=='__main__':main()
