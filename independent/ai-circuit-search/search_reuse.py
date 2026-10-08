"""Search existing partial sums for shorter exact, disjoint covers."""
import json
import sys
import experimental_search as exp
import search
import sidegen


def make_reuse(mode='pair',priority='large',minimum=3):
    Parent=exp.make_class(group_size=2,cutoff=4,split='power',loo='prefix',join='original')
    class Reuse(Parent):
        def total(self,vals):
            vals=[v for v in vals if v]
            if len(vals)<minimum:return super().total(vals)
            target=0
            for n in vals:
                assert target&self.support[n]==0
                target|=self.support[n]
            if target in self.lookup:return self.lookup[target]
            candidates=[(s,n) for s,n in self.lookup.items() if s and s!=target and s&~target==0]
            candidates.sort(key=lambda p:((-p[0].bit_count() if priority=='large' else p[0].bit_count()),p[1]))
            for s,n in candidates:
                other=self.lookup.get(target^s)
                if other:return self.add(n,other)
            if mode=='triple' and len(vals)>=4:
                for s,n in candidates[:32]:
                    rest=target^s
                    for t,j in candidates[:64]:
                        if t&~rest==0:
                            k=self.lookup.get(rest^t)
                            if k:return self.add(n,self.add(j,k))
            if mode=='greedy':
                remainder=target;cover=[]
                for s,n in candidates:
                    if s&~remainder==0:
                        cover.append(n);remainder^=s
                        if not remainder:break
                if not remainder and len(cover)<len(vals):
                    out=0
                    for n in cover:out=self.add(out,n)
                    return out
            return super().total(vals)
    return Reuse


def main():
    data=dict(search=[],best=None)
    for mode in ['pair','triple','greedy']:
        for priority in ['large','small']:
            for minimum in [3,4,6]:
                for h in [21,23,25]:
                    cls=make_reuse(mode,priority,minimum)
                    loc=cls(h-1);exp.exact_outputs(loc)
                    sidegen.Excl=cls
                    c=search.build(h,search.orders.gm);net=search.simple(c)
                    root=float(search.verify.approximate_root(net))
                    row=dict(mode=mode,priority=priority,minimum=minimum,h=h,saving=root,roles=c['R'])
                    data['search'].append(row)
                    if data['best'] is None or root>data['best']['saving']:
                        data['best']={**row,'network':net};print('BEST',row,flush=True)
                    print('tested',len(data['search']),row,flush=True)
                    (search.ROOT/'reuse-results.json').write_text(json.dumps(data,indent=2))


if __name__=='__main__':main()
