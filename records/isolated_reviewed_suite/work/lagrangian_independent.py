"""Independent exact-rational finite-ring verification of Lagrangian claims.
Does not invoke the manuscript's finite-image proof or its rigidity argument.
"""
from functools import lru_cache
from collections import defaultdict
import sympy as s
import json, time, argparse

@lru_cache(None)
def gaussian(n,k):
    if k<0 or k>n: return {}
    if k==0: return {0:1}
    out=dict(gaussian(n-1,k-1))
    for a,b in gaussian(n-1,k).items(): out[a+k]=out.get(a+k,0)+b
    return out

def target(n,m):
    out={0:1}
    for i in range(1,min(m,n)+1,2):
        for j in range(n-i+1):
            for a,b in gaussian(i+j,i).items():
                t=i+j*(j+1)//2+a
                out[t]=out.get(t,0)+b
    return out

def ambient(n):
    @lru_cache(None)
    def reduce(exp):
        pivot=next((i for i,a in enumerate(exp,1) if a>=2),None)
        if pivot is None:
            return {sum((1<<i) for i,a in enumerate(exp) if a):1}
        i=pivot; out=defaultdict(int); base=list(exp);base[i-1]-=2
        for a in range(1,min(i,n-i)+1):
            y=base[:];y[i+a-1]+=1
            if i>a:y[i-a-1]+=1
            for mask,v in reduce(tuple(y)).items():out[mask]+=2*(-1)**(a+1)*v
        return {k:v for k,v in out.items() if v}
    return reduce

def monomials(weights,total):
    if not weights:
        if total==0:yield ()
        return
    w=weights[-1]
    for a in range(total//w+1):
        for rest in monomials(weights[:-1],total-a*w):yield rest+(a,)

def filtration(n):
    reduce=ambient(n);N=n*(n+1)//2
    degree_masks=defaultdict(list)
    for mask in range(1<<n):degree_masks[sum(i+1 for i in range(n) if mask>>i&1)].append(mask)
    result=[]
    for m in range(1,n+1,2):
        actual={}
        for t in range(N+1):
            rows=degree_masks[t]; ridx={x:i for i,x in enumerate(rows)}; cols=[]
            for p in monomials(list(range(1,m+1,2)),t):
                exp=[0]*n
                for i,a in zip(range(1,m+1,2),p):exp[i-1]=a
                col=[0]*len(rows)
                for mask,v in reduce(tuple(exp)).items():col[ridx[mask]]=v
                cols.append(col)
            rank=s.polys.matrices.DomainMatrix.from_Matrix(s.Matrix.hstack(*(s.Matrix(c) for c in cols))).rank() if cols else 0
            if rank:actual[t]=rank
        assert actual==target(n,m),(n,m,actual,target(n,m))
        result.append({'n':n,'m':m,'dimension':sum(actual.values()),'hilbert':actual})
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--max-n',type=int,default=7);args=parser.parse_args()
    start=time.time()
    for n in range(1,args.max_n+1):
        for row in filtration(n):print(json.dumps(row),flush=True)
        print(json.dumps({'finished_n':n,'seconds':round(time.time()-start,2)}),flush=True)
