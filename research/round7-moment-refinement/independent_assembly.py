#!/usr/bin/env python3
"""Rational-only reproduction of unchanged kappa7 assembly and reciprocal ablation."""
import sys
if sys.flags.optimize:raise ValueError('Assertions must remain enabled')
sys.dont_write_bytecode=True
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import independent_moment as proof

HERE=Path(__file__).resolve().parent

def ceil(x):return -((-x.numerator)//x.denominator)
def floor(x):return x.numerator//x.denominator
def js(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {str(k):js(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [js(v) for v in x]
    return x

def assembly(ab):
    ac=Q(36926111,500000000000);beta=Q(1,1000);gap=Q(1,10**16)
    mc,sc=576,119453132304
    tau,sigma=1-ab,1-ac
    chi=tau+(1-beta)*max(sigma-tau,0);leaf=sigma+beta*(1-sigma)
    lam=max(tau,sigma,chi)+gap;lamp=max(lam,leaf)+gap
    eps=Q(floor(Q(10**16)/(2-lamp))-1,10**16)
    lo,hi,_=proof.logarithm(Q(sc))
    assert ceil(mc*lo)==ceil(mc*hi)
    c1=ceil(mc*lo)+2;x=max(3,ceil(eps*c1)+1)
    delta=Q(1,10**22)/max(1,x);poly=(1+x)*delta
    top=min(eps*(1-lamp),1-eps-poly,max(eps,1-eps)*ab)
    kappa=Q(floor(top*10**17)-1,10**17)
    conditions={
        'lambda above tau,sigma,chi':lam-max(tau,sigma,chi),
        'lambda_prime above lambda and leaf':lamp-max(lam,leaf),
        'lambda_prime below one':1-lamp,
        'beta in (0,1)':min(beta,1-beta),
        'guard':1+x-eps*c1,
        'precision':1+x-3*eps,
        'reservations':lamp-(2*eps-1)/eps,
        'record regime':1-eps,
        'delta in (0,1/8)':min(delta,Q(1,8)-delta),
        'independent crude logarithmic bound':c1-2-mc*hi,
    }
    margins={
        'prefix moves and top individual round':1-eps,
        'simultaneous butterflies':eps*(1-lamp),
        'fine-bit exposures (A)':1-eps-poly,
        'Gaussian maps (B)':1-eps-poly,
        'chirps, twists, scalar products':1-eps-poly,
        'packed polynomial products':eps,
        'individually processed reserved axes':1-eps-poly,
        'CRT axis reversal (E)':max(eps,1-eps)*ab,
    }
    assert all(v>0 for v in conditions.values())
    assert all(v>kappa for v in margins.values())
    return dict(bit_saving=ab,complex_saving=ac,beta=beta,eps=eps,c1=c1,x=x,
        lambda_value=lam,lambda_prime=lamp,delta=delta,kappa=kappa,
        binding=min(margins,key=margins.get),gap=min(margins.values())-kappa,
        strict_conditions=conditions,cost_margins=margins,
        crude_log_lower=lo,crude_log_upper=hi,crude_log_ceil=ceil(mc*lo))

def original_reciprocal_bounds(c,a):
    lower=upper=Q(0);terms=[]
    for t,count in c['hist']:
        ratio=Q(c['m'],t);l,u,_=proof.logarithm(ratio)
        k=0;y=ratio
        while y>2:y/=2;k+=1
        # Upstream30-term log upper differs from true log by no more than
        # this explicit tail bound, using z<=1/3 for each scaled term.
        tail=2*(k+1)*Q(1,3)**61/(61*(1-Q(1,3)**2))
        weight=Q(count*t,c['m']*c['W']);assert a*(u+tail)<1
        termlo=proof.down(weight/(1-a*l),proof.SUM_SCALE)
        termhi=proof.up(weight/(1-a*(u+tail)),proof.SUM_SCALE)
        lower+=termlo;upper+=termhi
        terms.append(dict(width=t,base2_steps=k,original_log_tail_bound=tail,lower=termlo,upper=termhi))
    return dict(a=a,lower=lower,upper=upper,terms=terms)

def main():
    import argparse
    from check_sources import check_sources
    check_sources()
    parser=argparse.ArgumentParser()
    parser.add_argument('--claim',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    claim=json.loads(args.claim.read_text());input_path=proof.ROOT/'lean/round7-histograms.json'
    data=json.loads(input_path.read_text());c=data['bit']
    assert claim['histogram']==c
    results={}
    for name,variant in claim['variants'].items():
        got=assembly(Q(variant['bit_saving']));expected=variant['assembly']
        for key in ('kappa','eps','c1','x','gap'):assert Q(got[key])==Q(expected[key]),(name,key)
        assert got['binding']==expected['binding'] and expected['ok'] and not expected['bad']
        results[name]=got;print(name,'PASS',got['kappa'],got['binding'],flush=True)
    a=Q(claim['reciprocal_fine']['saving']);b=Q(claim['reciprocal_fine']['next_saving'])
    assert b-a==Q(1,10**18)
    accepted=original_reciprocal_bounds(c,a);rejected=original_reciprocal_bounds(c,b)
    assert accepted['upper']<1<rejected['lower']
    out=dict(status='PASS independent original assembly and literal30-term reciprocal-majorant thresholds',
        claim_sha256=sha256(args.claim.read_bytes()).hexdigest(),histogram_sha256=sha256(input_path.read_bytes()).hexdigest(),
        checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        moment_helper_sha256=sha256((HERE/'independent_moment.py').read_bytes()).hexdigest(),
        original_assembly_source_sha256=sha256((proof.ROOT/'scripts/certificate_round3.py').read_bytes()).hexdigest(),
        variants=results,original_reciprocal_accepted=accepted,original_reciprocal_rejected=rejected,
        scope='Unchanged rational assembly verified independently; physical, genericity, stage-two and all-size hypotheses remain dependencies.')
    path=args.out;assert not path.exists()
    path.write_text(json.dumps(js(out),indent=2,sort_keys=True)+'\n')
    print('PASS original reciprocal accepted/next, including rigorous30-term log-tail allowance',flush=True)

if __name__=='__main__':main()
