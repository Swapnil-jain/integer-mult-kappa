"""Explore published histograms and alternative point orders; no theorem claim."""
import ast
import json
import math
import re
import sys
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT
REPO = ROOT.parents[1]
sys.path.insert(0,str(OUT))
import verify
sys.path.insert(0,str(REPO/'scripts'))
sys.path.insert(0,str(REPO/'independent'/'two-stage-bit'))
import certificate_round6 as c6
import certificate_round3 as r3
import cert_rt, rtgm, orders


def simple(c):
    return dict(m=c['m'],W=c['W'],hist={str(k):v for k,v in c['hist'].items()},rank_sum=c['s'])


def build(h, order):
    G=rtgm.build(h,order=order,retain='search')
    C=rtgm.compile_(G)
    rk=Counter()
    for s in range(C['roles']):
        d=rtgm.chain_dims(G,C,s)
        for a,b in zip(d,d[1:]): rk[b-a]+=1
    return cert_rt.rt_histogram(h,C['roles'],rk)


def main():
    nets={}
    _,c=c6.retained_saving()
    nets['jain_bit']=simple(c)
    lean=(REPO/'lean'/'Round6.lean').read_text()
    for label,name,m,W,s in [('bitHist','jain_lean_bit',529,148225616,78410006675),
                            ('cxHist','jain_complex',576,207387136,119453132304)]:
        hist=ast.literal_eval(re.search(r'def '+label+r'.*?:= (\[.*\])',lean).group(1))
        nets[name]=dict(m=m,W=W,hist={str(w):n for w,n in hist},rank_sum=s)
    assert nets['jain_bit']==nets['jain_lean_bit']
    del nets['jain_lean_bit']
    croc=json.loads((ROOT/'croc-certificate.json').read_text())
    for k in ['bit','complex']:
        c=croc[k]['counts']
        nets['community_'+k]=dict(m=c['m'],W=c['W'],rank_sum=c['total_rank'],hist=c['child_multiplicities'])
    data=dict(networks=nets,root_certificates={},jain_assemblies={},search=[])
    for name,net in nets.items():
        cert=verify.certify_root(net)
        data['root_certificates'][name]=cert
        print(name,cert['approximate_root'],flush=True)
    (ROOT/'base-roots.json').write_text(json.dumps(data,indent=2))
    modes=['gm','gmf','gmm','gmrev','natural','rot']
    best=0
    for h in range(18,33):
        for mode in modes:
            try:
                c=build(h,getattr(orders,mode));net=simple(c)
                # Float screening only. Certify the winner independently afterwards.
                m,W,hist=verify.normalize_network(net)
                low,high=0.,.001
                for _ in range(60):
                    a=(low+high)/2
                    v=sum(n*(w/m)**(1-a) for w,n in hist.items())/W
                    if v<1:low=a
                    else:high=a
                value=(low+high)/2
                data['search'].append(dict(h=h,order=mode,approximate_saving=value,roles=c['R']))
                if value>best:
                    best=value;data['search_best']=dict(h=h,order=mode,network=net)
                    print('new best',h,mode,'saving',value,'R',c['R'],flush=True)
            except AssertionError as e:
                data['search'].append(dict(h=h,order=mode,rejected=str(e)))
        print('finished h',h,flush=True)
        (ROOT/'search-results.json').write_text(json.dumps(data,indent=2))
    print('DONE',data['search_best']['h'],data['search_best']['order'],best,flush=True)


if __name__=='__main__':main()
