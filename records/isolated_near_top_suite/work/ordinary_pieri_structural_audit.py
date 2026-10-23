"""Independent exact audit using Pieri multiplication, no claimed filtration formula in construction."""
from functools import lru_cache
from itertools import combinations
from math import gcd
from functools import reduce
import json, time, sys


def partitions(k,l):
    out=[]
    def rec(pre,b):
        if len(pre)==k:out.append(tuple(pre));return
        for a in range(b+1):rec(pre+[a],a)
    rec([],l)
    return out

def gaussian(a,b):
    if b<0 or b>a:return []
    rows=[[1]]
    def add(x,y,s=0):
        z=[0]*max(len(x),len(y)+s)
        for i,v in enumerate(x):z[i]+=v
        for i,v in enumerate(y):z[i+s]+=v
        return z
    for n in range(1,a+1):
        prev=rows;rows=[]
        for r in range(n+1):
            if r==0 or r==n:rows.append([1]);continue
            rows.append(add(prev[r-1],prev[r],r))
    return rows[b]
def conv(a,b):
    z=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return z
def accum(a,b,s=0):
    if len(a)<len(b)+s:a.extend([0]*(len(b)+s-len(a)))
    for i,x in enumerate(b):a[i+s]+=x

def addvec(dst,src,m=1):
    for p,x in src.items():
        y=dst.get(p,0)+m*x
        if y:dst[p]=y
        else:dst.pop(p,None)

def primitive(v):
    if not v:return v
    g=reduce(gcd,(abs(x) for x in v.values()))
    if v[min(v)]<0:g=-g
    if g!=1:v={p:x//g for p,x in v.items()}
    return v

def insert(basis,v):
    v=dict(v)
    while v:
        p=min(v)
        if p not in basis:
            basis[p]=primitive(v);return True
        b=basis[p];a=v[p];bp=b[p]
        v={q:bp*x for q,x in v.items()}
        addvec(v,b,-a)
        v=primitive(v)
    return False

class Ring:
    def __init__(self,k,l):
        self.k=k;self.l=l;self.N=k*l
        self.ps=partitions(k,l);self.id={p:i for i,p in enumerate(self.ps)}
        self.w=[sum(p) for p in self.ps]
        self.by=[[] for _ in range(self.N+1)]
        for i,w in enumerate(self.w):self.by[w].append(i)
        self.unit={self.id[(0,)*k]:1}
        self.algebra_bases={};self.ideal_bases={}
    @lru_cache(None)
    def pieri(self,i,p,kind):
        if i==0:return (p,)
        w=self.w[p]+i
        if w>self.N:return ()
        a=self.ps[p]
        if kind=='e':
            return tuple(t for t in self.by[w] if all(0<=b-c<=1 for b,c in zip(self.ps[t],a)))
        return tuple(t for t in self.by[w] if all(self.ps[t][j]>=a[j] for j in range(self.k)) and all(a[j]>=self.ps[t][j+1] for j in range(self.k-1)))
    def apply(self,i,v,kind='e'):
        out={}
        for p,c in v.items():
            for t in self.pieri(i,p,kind):out[t]=out.get(t,0)+c
        return {p:c for p,c in out.items() if c}
    @lru_cache(None)
    def a_basis(self,n,p):
        out={}
        for i in range(max(0,n-self.l),min(n,self.l)+1):
            for q in self.pieri(n-i,p,'h'):
                for t in self.pieri(i,q,'h'):out[t]=out.get(t,0)+1
        return tuple(out.items())
    def apply_a(self,n,v):
        out={}
        for p,c in v.items():addvec(out,dict(self.a_basis(n,p)),c)
        return out
    def algebra(self,m):
        bs=[{} for _ in range(self.N+1)];insert(bs[0],self.unit)
        for n in range(1,self.N+1):
            for i in range(1,min(m,n)+1):
                for v in bs[n-i].values():insert(bs[n],self.apply(i,v))
        self.algebra_bases[m]=bs
        return [len(b) for b in bs]
    def minor(self,rows,s):
        L=self.l-self.k+1
        @lru_cache(None)
        def rec(cols):
            r=s-len(cols)
            if not cols:return tuple(self.unit.items())
            v={}
            for j,c in enumerate(cols):
                t=self.apply_a(L+rows[r]+c,dict(rec(cols[:j]+cols[j+1:])))
                addvec(v,t,(-1)**j)
            return tuple(v.items())
        return dict(rec(tuple(range(s))))
    def ideal(self,s):
        bs=[{} for _ in range(self.N+1)]
        for rows in combinations(range(self.k),s):
            v=self.minor(rows,s)
            if v:insert(bs[self.w[next(iter(v))]],v)
        for n in range(self.N+1):
            for i in range(1,min(self.k,n)+1):
                for v in bs[n-i].values():insert(bs[n],self.apply(i,v))
        self.ideal_bases[s]=bs
        return [len(b) for b in bs]
    def transformed_schur(self,lam):
        k=self.k
        @lru_cache(None)
        def rec(cols):
            r=k-len(cols)
            if not cols:return tuple(self.unit.items())
            v={}
            for j,c in enumerate(cols):
                n=lam[r]-r+c
                if n<0:continue
                t=self.apply_a(n,dict(rec(cols[:j]+cols[j+1:])))
                addvec(v,t,(-1)**j)
            return tuple(v.items())
        return dict(rec(tuple(range(k))))
    def target_increment(self,m):
        t=[0]*(self.N+1)
        for j in range(self.k-m+1):accum(t,conv(gaussian(self.l,m),gaussian(m+j-1,j)),m+j*(self.l-m+1))
        return t
    def target_polar(self,s):
        t=[0]*(self.N+1);L=self.l-self.k+1
        for r in range(s):accum(t,conv(gaussian(self.k,r),gaussian(self.l,self.k-r)),r*(r+L-1))
        return t

def check(k,l,ideals=True):
    st=time.time();R=Ring(k,l);actual=[[1]+[0]*R.N]
    for m in range(1,k+1):
        h=R.algebra(m);inc=[x-y for x,y in zip(h,actual[-1])]
        expected=R.target_increment(m)
        assert inc==expected,('FILTRATION MISMATCH',k,l,m,inc,expected)
        actual.append(h)
        print(json.dumps({'case':[k,l],'m':m,'dim_A':sum(h),'status':'increment exact'}),flush=True)
    polar=[]
    if ideals:
        for s in range(1,k+1):
            ranks=R.ideal(s);quot=[len(R.by[n])-ranks[n] for n in range(R.N+1)]
            expected=R.target_polar(s)
            assert quot==expected,('POLAR MISMATCH',k,l,s,quot,expected)
            polar.append({'s':s,'dimension':sum(quot),'hilbert':quot})
            print(json.dumps({'case':[k,l],'s':s,'dim_polar':sum(quot),'status':'Jacobian minors exact'}),flush=True)
    # Full-space checks: independence of all transformed Schurs, outside-hook
    # containment and bases of actual ideals, and ideal containment in A_m.
    if ideals:
        transformed=[{} for _ in range(R.N+1)]
        for lam in R.ps:
            v=R.transformed_schur(lam);n=sum(lam)
            assert insert(transformed[n],v),('TRANSFORMED BASIS DEPENDENCE',k,l,lam)
            durfee=max([0]+[i+1 for i,a in enumerate(lam) if a>=i+l-k+1])
            for s in range(1,durfee+1):
                probe=dict(R.ideal_bases[s][n])
                assert not insert(probe,v),('OUTSIDE HOOK NOT IN IDEAL',k,l,s,lam)
        for s in range(1,k+1):
            m=k-s+1
            for n,basis in enumerate(R.ideal_bases[s]):
                for v in basis.values():
                    probe=dict(R.algebra_bases[m][n])
                    assert not insert(probe,v),('IDEAL NOT IN SUBALGEBRA',k,l,s,m,n)
        print(json.dumps({'case':[k,l],'status':'transformed basis and all subspace containments exact'}),flush=True)
    out={'k':k,'ell':l,'ring_dimension':len(R.ps),'filtration_hilbert':actual,'polar':polar,'seconds':time.time()-st}
    print(json.dumps({'case':[k,l],'seconds':out['seconds'],'status':'PASS'}),flush=True)
    return out

if __name__=='__main__':
    cases=[tuple(map(int,a.split(','))) for a in sys.argv[1:]] or [(4,4),(4,5),(5,5),(5,6)]
    results=[]
    for k,l in cases:
        results.append(check(k,l))
        with open(__file__.replace('.py','_results.json'),'w') as f:json.dump(results,f,indent=2)
