import hashlib
import copy
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OUT=ROOT
sys.path.insert(0,str(OUT))
import verify
from community_assembly import assembly


def serial(value):
    if isinstance(value,Q):return str(value)
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v) for v in value]
    return value


def jain_parameters(ab,ac,network,improved):
    gap=Q(1,10**30) if improved else Q(1,10**16)
    beta=Q(1,1000)
    lp=max(1-ab+gap,1-ac+beta*ac)+gap
    if improved:
        eps=verify.floor_grid(1/(2-lp),10**30)-Q(1,10**30)
        delta=Q(1,10**35)
    else:
        eps=verify.floor_grid(1/(2-lp),10**16)-Q(1,10**16)
        delta=Q(1,10**22*14695)
    p=dict(a_bit=ab,a_complex=ac,beta=beta,epsilon=eps,x=14695,C1=14694,
           lambda_gap=gap,delta=delta,complex_m=576,complex_rank_sum=119453132304,
           bit_network=network,complex_network='jain_complex',kappa=Q(1,10**10))
    slacks,margins=verify.check_jain_assembly(p)
    p['kappa']=verify.floor_grid(min(margins.values()),10**24)-Q(1,10**24)
    slacks,margins=verify.check_jain_assembly(p)
    return p


def main():
    data=json.loads((ROOT/'base-roots.json').read_text())
    variants=json.loads((ROOT/'variant-results.json').read_text())
    candidates=json.loads((ROOT/'search-results.json').read_text())
    data['networks']['power_split_bit']=variants['best']['network']
    data['root_certificates']['power_split_bit']=verify.certify_root(variants['best']['network'])
    ac=Q(36926111,5*10**11)
    p=jain_parameters(Q(36667,10**9),ac,'jain_bit',False)
    p['kappa']=Q(3666565558019,10**17)
    verify.check_jain_assembly(p)
    data['jain_assemblies']['published']=p
    for name in ['jain_bit','power_split_bit']:
        ab=Q(data['root_certificates'][name]['lower'])
        data['jain_assemblies'][name+'_sharpened']=jain_parameters(ab,ac,name,True)
    croc=json.loads((ROOT/'croc-certificate.json').read_text())
    data['community_finite_bridge']=croc['finite_bridge']
    old=croc['assembly']['parameters']
    data['community_assemblies']={'published':{k:old[k] for k in ('a_bit','a_complex','kappa','h')}}
    a=Q(data['root_certificates']['community_bit']['lower'])
    # Prove all original strict inequalities at a smaller positive backoff.
    trial=assembly(croc['finite_bridge'],a,Q(old['kappa']),h=Q(1,10**24))
    kappa=verify.floor_grid(trial['minimum_margin'],10**24)-Q(1,10**24)
    trial=assembly(croc['finite_bridge'],a,kappa,h=Q(1,10**24))
    data['community_assemblies']['sharpened']={k:trial['parameters'][k] for k in ('a_bit','a_complex','kappa','h')}
    data['community_new_slacks']=trial['constraints']
    community_search=json.loads((ROOT/'community-search-results.json').read_text())
    net=community_search['best']['network']
    data['networks']['community_candidate_bit']=net
    data['root_certificates']['community_candidate_bit']=verify.certify_root(net)
    bridge=copy.deepcopy(croc['finite_bridge'])
    bridge['bit']['W']=net['W']
    assert net['W'].bit_length()==bridge['bit']['wire_bits']
    assert max(map(int,net['hist']))==bridge['bit']['maxchild']
    ab=Q(39554359,10**12)
    require_candidate=verify.moment_bounds(net,ab)[1]<1
    assert require_candidate
    kappa=Q(3955279,10**11)
    trial=assembly(bridge,ab,kappa,h=Q(1,10**12))
    params={k:trial['parameters'][k] for k in ('a_bit','a_complex','kappa','h')}
    params.update(bit_network='community_candidate_bit',finite_bridge=bridge)
    data['community_assemblies']['changed_circuit']=params
    data['community_candidate_slacks']=trial['constraints']
    data['community_candidate_configuration']={k:community_search['best'][k] for k in ('split','cutoff','vector','producer')}
    data['profile_limits']={}
    for name,cert in data['root_certificates'].items():
        if name.endswith('_bit'):
            upper=Q(cert['upper'])
            data['profile_limits'][name]=upper/(1+upper)
    counts={}
    for filename in ['search-results.json','variant-results.json','experimental-broad.json',
                     'experimental-smallbase.json','experimental-refine.json','reuse-results.json',
                     'evolution-results.json','community-search-results.json']:
        rows=json.loads((ROOT/filename).read_text())['search']
        counts[filename]=dict(evaluations=len(rows),rejected=sum('rejected' in r for r in rows))
    data['search_summary']=dict(batches=counts,total_evaluations=sum(r['evaluations'] for r in counts.values()),
        best_variant={k:v for k,v in variants['best'].items() if k!='network'},
        note='Exploratory floating scores; only separately enclosed moments are certificates.')
    data['scope']='Finite moment and assembly arithmetic plus audited candidate side circuit. Conditional on upstream fixed-tape, common-basis, analytic and exact-recovery arguments; no full formal proof, runtime or priority claim.'
    data['provenance']=dict(jain_commit='f2176bc1124821bf17eb63725bd366d7bdc020a3',
        community_commit='0605a24a28836168ad29d6239b46064b892298fc',
        community_certificate_sha256=hashlib.sha256((ROOT/'croc-certificate.json').read_bytes()).hexdigest(),
        jain_url='https://github.com/Swapnil-jain/integer-mult-kappa',
        community_url='https://github.com/CrocSwap/integer-mult-bounds')
    (OUT/'certificate.json').write_text(json.dumps(serial(data),indent=2),encoding='utf-8',newline='\n')
    for family in ('jain_assemblies','community_assemblies'):
        for name,row in data[family].items():
            k=Q(row['kappa']);base=Q(data[family]['published']['kappa'])
            print(family,name,str(k),float(k),'improvement_percent',float(100*(k/base-1)))
    print('limits',serial(data['profile_limits']))


if __name__=='__main__':main()
