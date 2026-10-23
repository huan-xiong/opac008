"""Full transformed-Q flag and all potentially nonzero strict generators.
Uses the independently implemented quadratic ambient reducer; sparse exact
integer row elimination is performed without a numerical rank tolerance.
"""
from lagrangian_independent import ambient
from functools import lru_cache
from itertools import combinations
from math import gcd
import json,time,argparse

def primitive(v):
    v={k:int(a) for k,a in v.items() if a}
    if not v:return v
    g=0
    for a in v.values():g=gcd(g,abs(a))
    if v[min(v)]<0:g=-g
    return {k:a//g for k,a in v.items()}

def reduce_row(v,rows):
    v=primitive(v)
    while v:
        p=min(v)
        if p not in rows:return v
        b=rows[p];a=v[p];c=b[p]
        v=primitive({k:c*v.get(k,0)-a*b.get(k,0) for k in v.keys()|b.keys()})
    return v

def add(a,b,c=1):
    o=dict(a)
    for k,v in b.items():o[k]=o.get(k,0)+c*v
    return {k:v for k,v in o.items() if v}

def strict_partitions(maxpart,total):
    yield ()
    def rec(prefix,nextmax,remaining):
        for a in range(1,min(nextmax,remaining)+1):
            p=prefix+(a,)
            yield p
            yield from rec(p,a-1,remaining-a)
    yield from rec((),maxpart,total)

def check(n):
    start=time.time();N=n*(n+1)//2;red=ambient(n)
    weights=[sum(i+1 for i in range(n) if m>>i&1) for m in range(1<<n)]
    e=[{0:1}]+[{1<<i:1} for i in range(n)]
    @lru_cache(None)
    def mm(a,b):
        if weights[a]+weights[b]>N:return {}
        return red(tuple(((a>>i)&1)+((b>>i)&1) for i in range(n)))
    def mul(a,b):
        out={}
        for x,u in a.items():
            for y,v in b.items():out=add(out,mm(x,y),u*v)
        return out
    q=[]
    for a in range(2*n+1):
        o={}
        for i in range(max(0,a-n),min(n,a)+1):o=add(o,mul(e[i],e[a-i]))
        q.append(o)
    @lru_cache(None)
    def pair(a,b):
        out=mul(q[a],q[b])
        for j in range(1,b+1):
            if a+j<=2*n:out=add(out,mul(q[a+j],q[b-j]),2*(-1)**j)
        return out
    @lru_cache(None)
    def pf(lam):
        if not lam:return {0:1}
        if sum(lam)>N:return {}
        if len(lam)%2:return pf(lam+(0,))
        out={}
        for j in range(1,len(lam)):
            out=add(out,mul(pair(lam[0],lam[j]),pf(lam[1:j]+lam[j+1:])),(-1)**(j+1))
        return out
    retained=list(strict_partitions(n,N))
    assert len(retained)==1<<n
    flag=[[{} for t in range(N+1)] for d in range(n+1)]
    allrows=[{} for t in range(N+1)]
    for lam in retained:
        vec=pf(lam);t=sum(lam)
        r=reduce_row(vec,allrows[t]);assert r,('basis',n,lam)
        allrows[t][min(r)]=r
        for d in range(len(lam)):
            r=reduce_row(vec,flag[d][t]);assert r
            flag[d][t][min(r)]=r
    closure=0
    for d in range(n+1):
        for t,rows in enumerate(flag[d]):
            for v in rows.values():
                for i in range(1,min(n,N-t)+1):
                    assert not reduce_row(mul(v,e[i]),flag[d][t+i]),('closure',n,d,t,i)
                    closure+=1
    extended=0
    for lam in strict_partitions(2*n,N):
        vec=pf(lam);t=sum(lam);extended+=1
        for d in range(min(len(lam),n+1)):
            assert not reduce_row(vec,flag[d][t]),('extended generator',n,d,lam)
    return {'n':n,'retained_basis':len(retained),'all_strict_candidates':extended,'closure_comparisons':closure,'flag_indices':n+1,'seconds':round(time.time()-start,3),'status':'PASS'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--max-n',type=int,default=8);args=p.parse_args()
    for n in range(1,args.max_n+1):print(json.dumps(check(n)),flush=True)
