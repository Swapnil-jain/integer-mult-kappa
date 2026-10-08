"""Regenerate and independently audit the changed side circuit (standard library).

Uses attributed, unmodified upstream producer modules in this repository. The only
algorithm change is grouping a sum at the largest power of two below its length,
instead of splitting it at floor(length/2). This is a conditional finite network
result, not a full formal verification of integer multiplication.
"""
import argparse
import hashlib
import json
import sys
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
sys.path[:0]=[str(REPO/'scripts'),str(REPO/'independent'/'two-stage-bit')]
import sidegen, rtgm, cert_rt
from verify import require, moment_bounds


class PowerSplit(sidegen.Excl):
    def total(self, vals):
        if not vals:return 0
        if len(vals)==1:return vals[0]
        cut=1 << ((len(vals)-1).bit_length()-1)
        return self.add(self.total(vals[:cut]),self.total(vals[cut:]))


def graph_rank(pairs):
    """Independent signless-incidence rank: vertices minus bipartite components."""
    neighbors={}
    for a,b in pairs:
        neighbors.setdefault(a,set()).add(b)
        neighbors.setdefault(b,set()).add(a)
    color={};bipartite=0
    for start in neighbors:
        if start in color:continue
        color[start]=1;queue=[start];ok=True
        for node in queue:
            for nxt in neighbors[node]:
                if nxt not in color:
                    color[nxt]=-color[node];queue.append(nxt)
                elif color[nxt]==color[node]:ok=False
        bipartite+=ok
    return len(neighbors)-bipartite


def set_bits(mask):
    while mask:
        bit=mask & -mask
        yield bit.bit_length()-1
        mask-=bit


def audit_graph(G,C):
    h=G['h'];trip=list(combinations(range(h),3))
    require(G['trip']==trip,'triple enumeration mismatch')
    support=G['sup'];args=G['args'];rank_cache={};rank_checks=0
    for node in sorted(G['active']):
        if args[node]:
            a,b=args[node]
            require(a<node and b<node,'graph not acyclic')
            require(support[a]&support[b]==0,'overlapping sum')
            require(support[node]==support[a]|support[b],'incorrect sum support')
        else:
            require(1<=node<=len(trip) and support[node]==1<<(node-1),'incorrect leaf')
        indices=list(set_bits(support[node]))
        common=set(trip[indices[0]])
        for i in indices[1:]:common.intersection_update(trip[i])
        require(bool(common),'node has no common centre')
        centre=min(common)
        pairs=tuple(sorted(tuple(v for v in trip[i] if v!=centre) for i in indices))
        if pairs not in rank_cache:rank_cache[pairs]=graph_rank(pairs)
        require(G['dn'][node]==rank_cache[pairs],'span-rank mismatch')
        rank_checks+=1
    stars={c:sum(1<<i for i,t in enumerate(trip) if c in t) for c in range(h)}
    for (c,t),node in G['outputs'].items():
        expected=stars[c]
        for v in t:
            if v!=c:expected &= ~stars[v]
        require(support[node]==expected,'incorrect exclusion output')
    for c,node in G['retained'].items():
        require(support[node]==stars[c],'incorrect retained total')
    rk=Counter()
    for slot,chain in enumerate(C['hold']):
        dims=[0]
        last_support=0
        for node in chain:
            require(last_support & ~support[node]==0,'non-nested role chain')
            last_support=support[node]
            if dims[-1]!=G['dn'][node]:dims.append(G['dn'][node])
        if slot in C['out']:
            target=C['out'][slot]
            require(last_support==support[G['outputs'][target]],'wrong role endpoint')
            dims.append(h-1)
        if slot in C['ret']:
            require(last_support==stars[C['ret'][slot]],'wrong retained role endpoint')
        dims.append(h)
        require(all(a<b for a,b in zip(dims,dims[1:])),'non-increasing role ranks')
        for a,b in zip(dims,dims[1:]):rk[b-a]+=1
    require(sum(r*n for r,n in rk.items())==h*C['roles'],'role-rank sum mismatch')
    return rk,dict(active_nodes=len(G['active']),rank_checks=rank_checks,roles=C['roles'],
                   exclusion_outputs=len(G['outputs']),retained_outputs=len(G['retained']))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log',type=Path)
    args=parser.parse_args()
    sidegen.Excl=PowerSplit
    h=23
    G=rtgm.build(h,retain='search');C=rtgm.compile_(G)
    rk,log=audit_graph(G,C)
    c=cert_rt.rt_histogram(h,C['roles'],rk)
    generated=dict(m=c['m'],W=c['W'],hist={str(k):v for k,v in c['hist'].items()},rank_sum=c['s'])
    data=json.loads((ROOT/'certificate.json').read_text())
    require(generated==data['networks']['power_split_bit'],'saved candidate histogram differs')
    a=Q(data['root_certificates']['power_split_bit']['lower'])
    require(moment_bounds(generated,a)[1]<1,'candidate moment fails')
    log.update(status='PASS',saving=str(a),scope='Side sums, span ranks, role nesting and moments; common-basis/tape/analytic transfer remains conditional')
    # A deliberately corrupted sum must be rejected by the circuit audit.
    corrupt=next(n for n in sorted(G['active']) if G['args'][n])
    old=G['sup'][corrupt];G['sup'][corrupt]^=1
    try:audit_graph(G,C)
    except ValueError:log['corrupted_support_rejected']=True
    else:raise ValueError('negative control was not rejected')
    finally:G['sup'][corrupt]=old
    if args.log:args.log.write_text(json.dumps(log,indent=2))
    print(json.dumps(log,indent=2))


if __name__=='__main__':main()
