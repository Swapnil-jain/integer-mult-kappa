"""Deterministic stochastic search over addition-tree parenthesizations."""
import collections
import json
import math
import random
import search
import sidegen

BASE=sidegen.Excl


def producer(cuts):
    class Evolved(BASE):
        sizes=collections.Counter()
        def total(self,vals):
            if not vals:return 0
            if len(vals)==1:return vals[0]
            type(self).sizes[len(vals)]+=1
            mid=cuts.get(len(vals),1<<((len(vals)-1).bit_length()-1))
            assert 0<mid<len(vals)
            return self.add(self.total(vals[:mid]),self.total(vals[mid:]))
    return Evolved


def evaluate(cuts,h=23):
    cls=producer(cuts);sidegen.Excl=cls
    c=search.build(h,search.orders.gm);net=search.simple(c)
    m,W,hist=search.verify.normalize_network(net)
    items=[(w*n/(m*W),math.log(m/w)) for w,n in hist.items()]
    deficit=(m*W-c['s'])/(m*W)
    lo,hi=0.,.001
    for _ in range(60):
        mid=(lo+hi)/2
        if sum(weight*math.expm1(mid*l) for weight,l in items)<deficit:lo=mid
        else:hi=mid
    return dict(cuts={str(k):v for k,v in cuts.items()},saving=(lo+hi)/2,roles=c['R'],
                sizes=dict(cls.sizes),network=net,h=h)


def main():
    rng=random.Random(817362)
    best=evaluate({});data=dict(seed=817362,search=[],best=best)
    active=sorted(int(k) for k in best['sizes'] if int(k)>2)
    print('active lengths',active,'baseline',best['saving'],flush=True)
    seen=set()
    for trial in range(400):
        cuts={int(k):v for k,v in best['cuts'].items()}
        # Mostly hill climbing, with periodic broader restarts and coupled changes.
        if trial%37==0:cuts={}
        count=1 if trial%4 else rng.randrange(2,min(5,len(active))+1)
        for length in rng.sample(active,count):cuts[length]=rng.randrange(1,length)
        key=tuple(sorted(cuts.items()))
        if key in seen:continue
        seen.add(key)
        row=evaluate(cuts)
        data['search'].append({k:v for k,v in row.items() if k not in ('network','sizes')})
        if row['saving']>best['saving']+1e-15:
            best=row;data['best']=row
            active=sorted(set(active)|{int(k) for k in row['sizes'] if int(k)>2})
            print('BEST',trial,row['cuts'],row['saving'],'roles',row['roles'],flush=True)
        if trial%25==24:print('progress',trial+1,'evaluated',len(seen),flush=True)
        (search.ROOT/'evolution-results.json').write_text(json.dumps(data,indent=2))
    print('DONE',len(seen),best['saving'],best['cuts'],flush=True)


if __name__=='__main__':main()
