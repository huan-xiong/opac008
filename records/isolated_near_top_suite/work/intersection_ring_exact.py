"""Independent top intersection checks from defining ring relations.

Lagrangian arithmetic: recursive squarefree reduction of the quadratic
presentation; normalize e_1...e_n to integral 1. No Kostka, Schur or Hilbert
formula is used by the reduction routine.
Ordinary arithmetic: SymPy Groebner reduction of the complete-function
relations; normalize e_k**ell to integral 1.
Comparison: integer Pieri recursion for Kostka numbers, and a separate
hook-length formula for the e_1**D normalization.
"""
from pathlib import Path
from functools import lru_cache
from itertools import product
from math import factorial
from datetime import datetime, timezone
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
START = datetime.now(timezone.utc).isoformat()

@lru_cache(None)
def weighted_monomials(n, D):
    if n == 0:
        return ((),) if D == 0 else ()
    result=[]
    for a in range(D//n+1):
        result.extend(v+(a,) for v in weighted_monomials(n-1,D-n*a))
    return tuple(result)

@lru_cache(None)
def kostka(shape,content):
    if sum(shape)!=sum(content): return 0
    if not content: return int(not shape)
    if not shape: return 0
    w=content[0]
    lower=shape[1:]+(0,)
    total=0
    for mu in product(*(range(lo,hi+1) for lo,hi in zip(lower,shape))):
        if sum(shape)-sum(mu)==w:
            # Interlacing guarantees a partition and horizontal strip.
            newshape=tuple(x for x in mu if x)
            total+=kostka(newshape,content[1:])
    return total

def content_of(exps):
    return tuple(i for i in range(len(exps),0,-1) for _ in range(exps[i-1]))

def tableau_count(shape):
    hooks=[]
    for i,length in enumerate(shape):
        for j in range(length):
            below=sum(row>j for row in shape[i+1:])
            hooks.append(length-j+below)
    denominator=1
    for h in hooks: denominator*=h
    return factorial(sum(shape))//denominator

def lg_reducer(n):
    @lru_cache(None)
    def reduce(exp):
        repeated=next((i for i in range(n,0,-1) if exp[i-1]>=2),None)
        if repeated is None:
            return {exp:1}
        i=repeated
        base=list(exp);base[i-1]-=2
        result={}
        for a in range(1,min(i,n-i)+1):
            v=base.copy();v[i+a-1]+=1
            if i-a: v[i-a-1]+=1
            factor=-2*((-1)**a)
            for mon,c in reduce(tuple(v)).items():
                result[mon]=result.get(mon,0)+factor*c
        return {mon:c for mon,c in result.items() if c}
    return reduce

records=[]
for n in range(1,7):
    D=n*(n+1)//2
    shape=tuple(range(n,0,-1))
    top=(1,)*n
    reduce=lg_reducer(n)
    monomials=weighted_monomials(n,D)
    for exp in monomials:
        nf=reduce(exp)
        assert not(set(nf)-{top})
        actual=nf.get(top,0)
        predicted=kostka(shape,content_of(exp))
        assert actual==predicted,(n,exp,actual,predicted)
    e1=(D,)+(0,)*(n-1)
    assert reduce(e1).get(top,0)==tableau_count(shape)
    rec={'kind':'Lagrangian','n':n,'top_weight':D,
         'top_monomials':len(monomials),'all_integrals':'PASS',
         'e1_top_integral':reduce(e1).get(top,0),
         'hook_tableaux':tableau_count(shape),
         'reduction_cache':reduce.cache_info()._asdict()}
    records.append(rec)
    print(rec,flush=True)

for k,ell in [(1,4),(2,2),(2,5),(3,3),(3,4),(4,4),(4,5)]:
    D=k*ell
    es=s.symbols('e1:'+str(k+1))
    h=[s.Integer(1)]
    for j in range(1,ell+k+1):
        h.append(s.expand(sum((-1)**(i+1)*es[i-1]*h[j-i]
                              for i in range(1,min(k,j)+1))))
    gb=s.groebner(h[ell+1:ell+k+1],*reversed(es),order='grevlex')
    top_poly=gb.reduce(es[-1]**ell)[1]
    assert top_poly!=0
    topP=s.Poly(top_poly,*es)
    pivot=topP.monoms()[0]
    pivot_coefficient=topP.coeff_monomial(pivot)
    shape=(k,)*ell
    monomials=weighted_monomials(k,D)
    for exp in monomials:
        f=s.prod(e**a for e,a in zip(es,exp))
        nf=gb.reduce(f)[1]
        actual=s.Poly(nf,*es).coeff_monomial(pivot)/pivot_coefficient
        assert s.expand(nf-actual*top_poly)==0,(k,ell,exp,'top component not scalar')
        predicted=kostka(shape,content_of(exp))
        assert actual==predicted,(k,ell,exp,actual,predicted)
    e1_nf=gb.reduce(es[0]**D)[1]
    e1_integral=s.Poly(e1_nf,*es).coeff_monomial(pivot)/pivot_coefficient
    assert e1_integral==tableau_count((ell,)*k)
    rec={'kind':'ordinary','k':k,'ell':ell,'top_weight':D,
         'top_monomials':len(monomials),'all_integrals':'PASS',
         'e1_top_integral':int(e1_integral),
         'hook_tableaux':tableau_count((ell,)*k),
         'groebner_basis_size':len(gb.polys)}
    records.append(rec)
    print(rec,flush=True)

out={'status':'PASS','started_at':START,'ended_at':datetime.now(timezone.utc).isoformat(),
     'sympy':s.__version__,
     'Lagrangian_top_monomials':sum(r['top_monomials'] for r in records if r['kind']=='Lagrangian'),
     'ordinary_top_monomials':sum(r['top_monomials'] for r in records if r['kind']=='ordinary'),
     'hook_checks':len(records),
     'scope':'Every elementary-generator monomial in top weight for listed rings; no full filtration ranks',
     'independence':'Defining-relation reduction compared with separately implemented integer horizontal-strip recursion; hook formula is a third calculation.',
     'records':records,'kostka_cache':kostka.cache_info()._asdict()}
(ROOT/'intersection_ring_exact_results.json').write_text(json.dumps(out,indent=2),encoding='utf8')
print('ALL PASS',out['Lagrangian_top_monomials'],out['ordinary_top_monomials'],flush=True)
