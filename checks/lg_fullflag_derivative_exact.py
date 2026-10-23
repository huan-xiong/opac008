"""Independent full transformed-Q flags and derivative contractions in R_n.
Uses raw quadratic squarefree reduction and a local Pfaffian recurrence.
No old helper, predicted ideal dimension, or initial largest-part cutoff is used.
"""
from pathlib import Path
from functools import lru_cache
from itertools import product
from datetime import datetime,timezone
import json
import sympy as s
ROOT=Path(__file__).resolve().parent
START=datetime.now(timezone.utc).isoformat()

@lru_cache(None)
def monomials(weights,t):
    if not weights: return ((),) if t==0 else ()
    w=weights[-1]
    return tuple(v+(a,) for a in range(t//w+1) for v in monomials(weights[:-1],t-w*a))

@lru_cache(None)
def strict_parts(t,maximum):
    if t==0: return ((),)
    return tuple((a,)+rest for a in range(min(t,maximum),0,-1)
                 for rest in strict_parts(t-a,a-1))

@lru_cache(None)
def gaussian(n,r):
    if not 0<=r<=n:return ()
    if r==0 or r==n:return (1,)
    left=gaussian(n-1,r-1);right=gaussian(n-1,r)
    v=[0]*max(len(left),len(right)+r)
    for i,a in enumerate(left):v[i]+=a
    for i,a in enumerate(right):v[i+r]+=a
    return tuple(v)

def target_increment(n,m,t):
    if m%2==0:return 0
    value=0
    for j in range(n-m+1):
        shift=m+j*(j+1)//2
        v=gaussian(m+j,m)
        if 0<=t-shift<len(v):value+=v[t-shift]
    return value

records=[]
for n in range(1,6):
    N=n*(n+1)//2
    zero=(0,)*n
    unit={zero:1}
    @lru_cache(None)
    def reduce(exp):
        i=next((i for i in range(n,0,-1) if exp[i-1]>=2),None)
        if i is None:return {exp:1}
        base=list(exp);base[i-1]-=2
        ans={}
        for a in range(1,min(i,n-i)+1):
            new=base.copy();new[i+a-1]+=1
            if i-a:new[i-a-1]+=1
            for mon,c in reduce(tuple(new)).items():
                ans[mon]=ans.get(mon,0)-2*((-1)**a)*c
        return {mon:c for mon,c in ans.items() if c}
    def add(*terms):
        ans={}
        for scalar,v in terms:
            for mon,c in v.items():ans[mon]=ans.get(mon,0)+scalar*c
        return {mon:c for mon,c in ans.items() if c}
    def multiply(a,b):
        ans={}
        for ea,ca in a.items():
            for eb,cb in b.items():
                for mon,c in reduce(tuple(x+y for x,y in zip(ea,eb))).items():
                    ans[mon]=ans.get(mon,0)+ca*cb*c
        return {mon:c for mon,c in ans.items() if c}
    e=[unit]
    for i in range(n):
        v=list(zero);v[i]=1;e.append({tuple(v):1})
    e.extend([{}]*n)
    q=[add(*((1,multiply(e[i],e[r-i])) for i in range(r+1)))
       for r in range(2*n+1)]
    def q_at(r):return q[r] if 0<=r<=2*n else {}
    @lru_cache(None)
    def pair(a,b):
        return add((1,multiply(q_at(a),q_at(b))),
                   *((2*((-1)**j),multiply(q_at(a+j),q_at(b-j))) for j in range(1,b+1)))
    @lru_cache(None)
    def pf(parts):
        if not parts:return unit
        assert len(parts)%2==0
        return add(*(((-1)**(j+1),multiply(pair(parts[0],parts[j]),pf(parts[1:j]+parts[j+1:])))
                     for j in range(1,len(parts))))
    def Z(parts):return pf(parts+(0,) if len(parts)%2 else parts)
    basis={t:[x for x in product([0,1],repeat=n) if sum((i+1)*a for i,a in enumerate(x))==t]
           for t in range(N+1)}
    def matrix(vs,t):
        return s.Matrix(len(basis.get(t,[])),len(vs),
                        lambda i,j:vs[j].get(basis[t][i],0))
    fullparts={t:strict_parts(t,t) for t in range(N+1)}
    generators={t:[(parts,Z(parts)) for parts in fullparts[t]] for t in range(N+1)}
    transformed_checks=0;flags=0
    ideal_spaces={}
    candidate_total=sum(len(v) for v in fullparts.values())
    outside_retained=sum(1 for vs in fullparts.values() for parts in vs if parts and parts[0]>n)
    for t in range(N+1):
        retained=[v for parts,v in generators[t] if not parts or parts[0]<=n]
        assert matrix(retained,t).rank()==len(basis[t])
        transformed_checks+=1
        for d in range(n+1):
            vectors=[]
            for wt in range(t+1):
                for parts,zval in generators[wt]:
                    if len(parts)>d:
                        vectors.extend(multiply(zval,{b:1}) for b in basis[t-wt])
            I=matrix(vectors,t)
            rank=I.rank()
            expected=sum(1 for b in basis[t] if sum(b)>d)
            assert rank==expected,(n,d,t,rank,expected)
            retained_long=[v for parts,v in generators[t] if len(parts)>d and parts[0]<=n]
            assert matrix(retained_long,t).rank()==expected
            if rank: I=s.Matrix.hstack(*I.columnspace())
            else: I=s.zeros(len(basis[t]),0)
            ideal_spaces[d,t]=I
            flags+=1
    subalgebras=0;derivatives=0
    for m in range(1,n+1):
        weights=tuple(range(1,m+1))
        for t in range(N+1):
            vs=[]
            for exp in monomials(weights,t):
                full=exp+(0,)*(n-m);vs.append(reduce(full))
            dim=matrix(vs,t).rank()
            expected=int(t==0)+sum(target_increment(n,i,t) for i in range(1,m+1))
            assert dim==expected,(n,m,t,dim,expected)
            subalgebras+=1
        if m%2==0:continue
        odd=tuple(range(1,m+1,2))
        for t in range(N+1):
            outputs=monomials(odd,t);inputs=monomials(odd,t+m)
            def lift(exp):
                v=[0]*n
                for i,a in zip(odd,exp):v[i-1]=a
                return tuple(v)
            # All kernel source weights are t+m, including those above the socle.
            source_vectors=[reduce(lift(exp)) for exp in inputs]
            src=matrix(source_vectors,t+m)
            kernel=src.nullspace()
            D=s.zeros(len(outputs),len(inputs))
            output_index={a:i for i,a in enumerate(outputs)}
            for j,a in enumerate(inputs):
                if a[-1]:
                    b=a[:-1]+(a[-1]-1,)
                    D[output_index[b],j]=a[-1]
            W=D*s.Matrix.hstack(*kernel) if kernel else s.zeros(len(outputs),0)
            values=matrix([reduce(lift(exp)) for exp in outputs],t)
            I=ideal_spaces[n-m,t]
            actualQdim=len(outputs)-(s.Matrix.hstack(I,values).rank()-I.cols)
            assert W.rank()==actualQdim,(n,m,t,'dimension',W.rank(),actualQdim)
            assert s.Matrix.hstack(I,values*W).rank()==I.cols,(n,m,t,'inclusion')
            derivatives+=1
    records.append({'n':n,'socle_weight':N,'strict_candidates_through_socle':candidate_total,
                    'candidates_largest_part_gt_n':outside_retained,
                    'full_transformed_basis_components':transformed_checks,
                    'full_flag_quotient_components':flags,
                    'subalgebra_components':subalgebras,
                    'derivative_components':derivatives,'status':'PASS'})
    print(records[-1],flush=True)

out={'status':'PASS','started_at':START,'ended_at':datetime.now(timezone.utc).isoformat(),
     'sympy':s.__version__,'records':records,
     'counts':{key:sum(r[key] for r in records) for key in ['strict_candidates_through_socle',
         'candidates_largest_part_gt_n','full_transformed_basis_components',
         'full_flag_quotient_components','subalgebra_components','derivative_components']},
     'scope':'n<=5, all flag indices and output weights through proven socle; derivative sources extend by newest odd weight.',
     'independence':'Self-contained raw quadratic quotient reduction and explicit E(z)^2 coefficient Pfaffians; no imported old helper or formula-defined ideal.'}
(ROOT/'lg_fullflag_derivative_exact_results.json').write_text(json.dumps(out,indent=2),encoding='utf8')
print(json.dumps(out,indent=2),flush=True)
