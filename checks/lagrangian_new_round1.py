"""Independent exact checks of Lagrangian formulas; no old audit code imported."""
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations
from collections import defaultdict
import sympy as s
import json, time
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def add(a,b,c=1):
    z=a.copy()
    for k,v in b.items():
        z[k]=z.get(k,0)+c*v
        if z[k]==0: del z[k]
    return z

class Space:
    def __init__(self): self.piv={}
    def insert(self,v):
        v={k:F(x) for k,x in v.items() if x}
        while v:
            k=min(v)
            if k not in self.piv:
                c=v[k]; self.piv[k]={i:x/c for i,x in v.items()}; return True
            v=add(v,self.piv[k],-v[k])
        return False
    def contains(self,v):
        v=v.copy()
        while v:
            k=min(v)
            if k not in self.piv:return False
            v=add(v,self.piv[k],-v[k])
        return True
    def __len__(self):return len(self.piv)

class LG:
    def __init__(self,n): self.n=n; self.N=n*(n+1)//2
    @lru_cache(None)
    def nf(self,p):
        for k,a in enumerate(p):
            if a>=2:break
        else:return {sum((1<<k) for k,a in enumerate(p) if a):1}
        i=k+1; p=list(p);p[k]-=2;ans={}
        for a in range(1,min(i,self.n-i)+1):
            t=p.copy();t[i+a-1]+=1
            if i>a:t[i-a-1]+=1
            ans=add(ans,self.nf(tuple(t)),2 if a%2 else -2)
        return ans
    @lru_cache(None)
    def mulmask(self,a,b):
        return self.nf(tuple(((a>>i)&1)+((b>>i)&1) for i in range(self.n)))
    def mul(self,a,b):
        z={}
        for i,x in a.items():
            for j,y in b.items():z=add(z,self.mulmask(i,j),x*y)
        return z
    def e(self,i):return {0:1} if i==0 else ({1<<(i-1):1} if i<=self.n else {})
    def weight(self,mask):return sum(i+1 for i in range(self.n) if (mask>>i)&1)
    def coeff_square(self,r):
        z={}
        for i in range(r+1):z=add(z,self.mul(self.e(i),self.e(r-i)))
        return z
    @lru_cache(None)
    def qpair(self,a,b):
        z=self.mul(self.coeff_square(a),self.coeff_square(b))
        for j in range(1,b+1):z=add(z,self.mul(self.coeff_square(a+j),self.coeff_square(b-j)),2*(-1)**j)
        return z
    @lru_cache(None)
    def zlambda(self,lam):
        if not lam:return {0:1}
        if len(lam)%2:lam=lam+(0,)
        ans={}
        for j in range(1,len(lam)):
            ans=add(ans,self.mul(self.qpair(lam[0],lam[j]),self.zlambda(lam[1:j]+lam[j+1:])),(-1)**(j+1))
        return ans
    def algebra(self,gens):
        spaces=defaultdict(Space); spaces[0].insert({0:1}); queue=[(0,{0:1})]
        for wt,v in queue:
            for gw,g in gens:
                if wt+gw>self.N:continue
                z=self.mul(v,g)
                if spaces[wt+gw].insert(z):queue.append((wt+gw,z))
        return {wt:len(sp) for wt,sp in spaces.items() if len(sp)}

def gbin(n,k):
    if k<0 or k>n:return [0]
    row=[[1]]
    for a in range(1,n+1):
        new=[]
        for b in range(a+1):
            if b==0 or b==a:new.append([1]);continue
            x=row[b-1];y=[0]*b+row[b];z=[0]*max(len(x),len(y))
            for i,t in enumerate(x):z[i]+=t
            for i,t in enumerate(y):z[i]+=t
            new.append(z)
        row=new
    return row[k]

def targets(n,m):
    ans={0:1}
    for i in range(1,m+1,2):
        for j in range(n-i+1):
            sh=i+j*(j+1)//2
            for k,x in enumerate(gbin(i+j,i)):
                if x:ans[k+sh]=ans.get(k+sh,0)+x
    return ans

def stricts(total,hi=None):
    if total==0:yield ();return
    if hi is None:hi=total
    for x in range(min(hi,total),0,-1):
        for tail in stricts(total-x,x-1):yield (x,)+tail

def lg_checks(n):
    A=LG(n);out={'n':n,'filtrations':[]}
    for m in range(1,n+1):
        got=A.algebra([(i,A.e(i)) for i in range(1,m+1,2)])
        assert got==targets(n,m),(n,m,got,targets(n,m))
        out['filtrations'].append({'m':m,'dimension':sum(got.values())})
    if n<=6:
        flag=[defaultdict(Space) for _ in range(n+1)]
        retained=[]
        for mask in range(1<<n):
            lam=tuple(i+1 for i in range(n-1,-1,-1) if (mask>>i)&1)
            wt=sum(lam);z=A.zlambda(lam);retained.append((lam,wt,z))
            for d in range(len(lam)):
                assert flag[d][wt].insert(z),(n,d,lam,'dependence')
        count=0
        for wt in range(A.N+1):
            for lam in stricts(wt):
                z=A.zlambda(lam);count+=1
                for d in range(min(len(lam),n+1)):
                    assert flag[d][wt].contains(z),(n,d,lam,'outside generator')
        for d in range(n+1):
            for lam,wt,z in retained:
                if len(lam)<=d:continue
                for i in range(1,n+1):
                    assert flag[d][wt+i].contains(A.mul(z,A.e(i))),(n,d,lam,i,'not ideal')
        out['all_strict_partitions_checked']=count
        out['flag_checked']=True
    return out

def transition_check(d,odd):
    M=s.Symbol('M');g=(s.Integer(1),)+s.symbols('g1:'+str(d+1));u={t:s.eye(d)[t-1,:] for t in range(1,d+1)};v={}
    max_t=2*d-1 if odd else max(d,2*d-2)
    for t in range(1,max_t+1):
        if t>d:
            u[t]=s.zeros(1,d)
            for i in range(1,d+1):u[t]+=-g[i]*u[t-i]+i*g[i]*v[t-i]
        v[t]=-u[t]/(M+t) if t%2==int(odd) else s.zeros(1,d)
    if odd:rows=[v[t] for t in range(1,2*d,2)];den=s.prod(M+2*j-1 for j in range(1,d+1));sign=(-1)**(d*(d+1)//2)
    else:
        b=u[d].copy()
        for i in range(1,d):b+=g[i]*u[d-i]-i*g[i]*v[d-i]
        rows=[v[t] for t in range(2,2*d,2)]+[b/(M+d)];den=(M+d)*s.prod(M+2*j for j in range(1,d));sign=(-1)**((d-1)*(d-2)//2)
    mat=s.Matrix.vstack(*rows);hur=s.det(s.Matrix(d-1,d-1,lambda a,j:g[d+a+1-2*(j+1)] if 0<=d+a+1-2*(j+1)<=d else 0))
    # Clear each row denominator then Bareiss avoids expensive rational expansion.
    scaled=mat.copy()
    if odd:
        for a in range(d):scaled[a,:]=scaled[a,:]*(M+2*a+1)
    else:
        for a in range(d-1):scaled[a,:]=scaled[a,:]*(M+2*a+2)
        scaled[d-1,:]=scaled[d-1,:]*(M+d)
    assert s.cancel(mat.det(method='domain-ge')*den-sign*hur)==0,(d,odd)
    return {'d':d,'M_parity':'odd' if odd else 'even','determinant':'verified symbolically in M'}

if __name__=='__main__':
    out={'started':time.time(),'filtration_and_flag':[],'transition':[]}
    for n in range(1,9):
        z=lg_checks(n);out['filtration_and_flag'].append(z);print(json.dumps(z),flush=True)
    for d in range(1,6):
        for odd in (False,True):
            z=transition_check(d,odd);out['transition'].append(z);print(json.dumps(z),flush=True)
    out['finished']=time.time();(ROOT/'lagrangian_new_round1_results.json').write_text(json.dumps(out,indent=2));print('ALL CHECKS PASSED',flush=True)
