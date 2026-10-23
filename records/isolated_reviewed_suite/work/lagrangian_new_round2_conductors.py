"""Compute full finite conductors from ALL ambient multiplication constraints.
The positive cases and negative boundary controls are both intentional.
"""
import sympy as s
from functools import lru_cache
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json
from lagrangian_new_round1 import Space,add
from lagrangian_new_round2_ordinary import conv,power_one
ROOT=Path(__file__).resolve().parent

def conductor_calc(vv,weights,relations,E,rho):
    gb=s.groebner(relations,*vv,order='grevlex',domain=s.QQ)
    lm=[p.LM(order='grevlex').exponents for p in gb.polys]
    zero=(0,)*len(vv);basis=[zero];seen={zero}
    for a in basis:
        for i in range(len(vv)):
            b=list(a);b[i]+=1;b=tuple(b)
            if b in seen:continue
            seen.add(b)
            if not any(all(x>=y for x,y in zip(b,p)) for p in lm):basis.append(b)
    wt=lambda a:sum(x*y for x,y in zip(a,weights))
    def expr(a):return s.prod(x**k for x,k in zip(vv,a))
    def vector(p):return {a:F(c) for a,c in s.Poly(gb.reduce(s.expand(p))[1],*vv).terms() if c}
    @lru_cache(None)
    def mon_nf(a):return vector(expr(a))
    def multiply(a,b):
        out={}
        for p,x in a.items():
            for q,y in b.items():out=add(out,mon_nf(tuple(i+j for i,j in zip(p,q))),x*y)
        return out
    Bspace=defaultdict(Space);Bbasis=defaultdict(list);queue=[(0,{zero:F(1)})];Bspace[0].insert(queue[0][1]);Bbasis[0].append(queue[0][1])
    gens=[(i,vector(p)) for i,p in enumerate(E) if i>0]
    for w,v in queue:
        for iw,g in gens:
            z=multiply(v,g)
            if Bspace[w+iw].insert(z):queue.append((w+iw,z));Bbasis[w+iw].append(z)
    def remainder(v,w):
        out={}
        while v:
            p=min(v)
            if p in Bspace[w].piv:v=add(v,Bspace[w].piv[p],-v[p])
            else:
                out[p]=v[p];v=v.copy();del v[p]
        return out
    conductor={}
    for w,cols in Bbasis.items():
        constraints={}
        for j,c in enumerate(cols):
            for bi,b in enumerate(basis):
                v=remainder(multiply(c,{b:F(1)}),w+wt(b))
                for a,co in v.items():
                    key=(bi,a)
                    if key not in constraints:constraints[key]=[F(0)]*len(cols)
                    constraints[key][j]=co
        rank=s.Matrix(list(constraints.values())).rank() if constraints else 0
        if len(cols)-rank:conductor[w]=len(cols)-rank
    rv=vector(rho);rhospace=defaultdict(Space)
    for b in basis:
        v=multiply(rv,{b:F(1)});w=wt(b)+(wt(next(iter(rv))) if rv else 0)
        if v:
            assert Bspace[w].contains(v),'rhoT not in coefficient algebra'
            rhospace[w].insert(v)
    rhilb={w:len(p) for w,p in rhospace.items() if len(p)}
    return {'dim_T':len(basis),'dim_B':sum(len(p) for p in Bspace.values()),'conductor_hilbert':conductor,'rho_hilbert':rhilb,'full_conductor_equals_rhoT':conductor==rhilb}

def ordinary(d,e,m):
    g=(s.Integer(1),)+s.symbols('g1:'+str(d+1));b=(s.Integer(1),)+s.symbols('b1:'+str(e+1));h=d+e;N=m+h;D=power_one(conv(g,b,N),-s.Rational(1,2),N)
    C=s.Matrix(h,h,lambda r,j:(g[d+r-j] if 0<=d+r-j<=d else 0) if r<e else (b[e+(r-e)-j] if 0<=e+(r-e)-j<=e else 0))
    rho=(-1)**(d*e)*C.det();E=conv(g,D,m+d)
    out=conductor_calc(g[1:]+b[1:],list(range(1,d+1))+list(range(1,e+1)),D[m+1:m+h+1],E,rho);out.update({'case':'ordinary','d':d,'e':e,'m':m})
    assert out['full_conductor_equals_rhoT']==(m>=1)
    return out

def lagrangian(d,M):
    g=(s.Integer(1),)+s.symbols('g1:'+str(d+1));r=M//2;n=M+d-1;N=M+2*d+2
    Gminus=[(-1)**i*x for i,x in enumerate(g)];D=power_one(conv(g,Gminus,N),-s.Rational(1,2),N)
    rel=[D[i] for i in range(M+1,M+2*d,2)] if M%2 else [D[i] for i in range(M+2,M+2*d,2)]+[g[d]*D[M]]
    rho=s.det(s.Matrix(d-1,d-1,lambda a,j:g[d+a+1-2*(j+1)] if 0<=d+a+1-2*(j+1)<=d else 0))
    out=conductor_calc(g[1:],list(range(1,d+1)),rel,conv(g,D,n),rho);out.update({'case':'lagrangian','d':d,'M':M})
    assert out['full_conductor_equals_rhoT']==(M>=2)
    return out
if __name__=='__main__':
    out=[]
    for d,e,m in [(1,1,0),(1,1,1),(1,2,1),(1,2,2),(2,2,1),(2,3,1)]:
        z=ordinary(d,e,m);out.append(z);print(json.dumps(z),flush=True);(ROOT/'lagrangian_new_round2_conductors_results.json').write_text(json.dumps(out,indent=2))
    for d,M in [(2,1),(2,2),(2,3),(3,2)]:
        z=lagrangian(d,M);out.append(z);print(json.dumps(z),flush=True);(ROOT/'lagrangian_new_round2_conductors_results.json').write_text(json.dumps(out,indent=2))
    print('ALL FULL-CONDUCTOR AND NEGATIVE-BOUNDARY CHECKS PASSED',flush=True)
