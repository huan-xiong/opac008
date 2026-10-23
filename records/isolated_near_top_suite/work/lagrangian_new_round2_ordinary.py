"""Targeted independent tests: asymmetric resultant matrices and monic cancellation.
No ordinary audit code is imported. Series powers use direct binomial expansion,
not the differential recurrence being audited.
"""
import sympy as s
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent

def conv(a,b,N):
    return [s.expand(sum(a[i]*b[n-i] for i in range(max(0,n-len(b)+1),min(n,len(a)-1)+1))) for n in range(N+1)]

def power_one(A,p,N):
    a=list(A)+[s.Integer(0)]*(N+1-len(A));a=a[:N+1];a[0]-=1
    v=[s.Integer(1)]+[s.Integer(0)]*N;D=[s.Integer(1)]+[s.Integer(0)]*N
    for j in range(1,N+1):
        v=conv(v,a,N);c=s.binomial(p,j)
        D=[s.expand(x+c*y) for x,y in zip(D,v)]
    return D

def trans(d,e):
    M=s.Symbol('M');h=d+e;g=(s.Integer(1),)+s.symbols('g1:'+str(d+1));b=(s.Integer(1),)+s.symbols('b1:'+str(e+1));uu={};vv={};dd={}
    for t in range(1,h+1):
        uu[t]=s.eye(h)[t-1,:] if t<=d else sum((-g[i]*uu[t-i]+i*g[i]*dd[t-i] for i in range(1,d+1)),s.zeros(1,h))
        vv[t]=s.eye(h)[d+t-1,:] if t<=e else sum((-b[i]*vv[t-i]+i*b[i]*dd[t-i] for i in range(1,e+1)),s.zeros(1,h))
        dd[t]=-(uu[t]+vv[t])/(2*(M+t))
    A=s.Matrix.vstack(*(dd[t] for t in range(1,h+1)))
    C=s.Matrix(h,h,lambda r,j:(g[d+r-j] if 0<=d+r-j<=d else 0) if r<e else (b[e+(r-e)-j] if 0<=e+(r-e)-j<=e else 0))
    rho=(-1)**(d*e)*C.det();scaled=s.cancel(A.det(method='domain-ge')*2**h*s.prod(M+t for t in range(1,h+1))/rho)
    assert scaled in (1,-1),(d,e,scaled)
    return g,b,M,A,rho,scaled

def identities(d,e,m):
    g,b,M,A,rho,sign=trans(d,e);N=m+d+e;Q=conv(g,b,N);D=power_one(Q,-s.Rational(1,2),N)
    invG=power_one(g,-1,N);invB=power_one(b,-1,N)
    zg=[s.Integer(0)]+[i*g[i] for i in range(1,d+1)];zb=[s.Integer(0)]+[i*b[i] for i in range(1,e+1)]
    U=conv(conv(zg,D,N),invG,N);V=conv(conv(zb,D,N),invB,N)
    jj=s.Matrix(U[m+1:m+d+1]+V[m+1:m+e+1]);ff=s.Matrix(D[m+1:N+1]);assert all(s.expand(x)==0 for x in ff-A.subs(M,m)*jj),(d,e,m)
    # The binomial-generated coefficients obey the independently constructed matrix.
    variables=g[1:]+b[1:];gb=s.groebner(list(ff),*variables,order='grevlex',domain=s.QQ)
    for j in jj:assert gb.reduce(s.expand(rho*j))[1]==0,(d,e,m,'colon containment')
    return {'d':d,'e':e,'m':m,'symbolic_det_ratio_sign':int(sign),'coefficient_identity':'pass','rho_J_subset_I':'pass'}

def nonreduced_witness():
    g,b,c,be,t=s.symbols('g b c beta t');N=4;G=[1,g];B=[1,b,c];D=power_one(conv(G,B,N),-s.Rational(1,2),N)
    source=s.groebner(D[2:5],g,b,c,domain=s.QQ,order='grevlex')
    f=t*t-be*t/2+3*be*be/8;ext=s.groebner([be**3,f],t,be,domain=s.QQ,order='grevlex')
    sub={g:-t,b:be-t,c:-be*t}
    for a in D[2:5]:assert ext.reduce(s.expand(a.subs(sub,simultaneous=True)))[1]==0
    assert s.expand((g*g-b*g+c).subs(sub,simultaneous=True))==0
    assert ext.reduce(be**2)[1]!=0,'nilpotent lost'
    E=conv(G,D,4);Dp=[1,-be/2,3*be*be/8,0,0]
    for i in range(5):assert ext.reduce(s.expand(E[i].subs(sub,simultaneous=True)-Dp[i]))[1]==0
    # Establish that the resulting source factorization is finite in the ring itself.
    Ep=E[:3];Hp=conv(B,D,4)[:4];prod=conv(Ep,Hp,5)
    assert source.reduce(s.expand(prod[0]-1))[1]==0
    for a in prod[1:]:assert source.reduce(a)[1]==0
    return {'d':1,'e':2,'m':1,'smaller_ring':'Q[beta]/(beta^3)','monic_polynomial':'t^2-beta*t/2+3*beta^2/8','beta_squared_nonzero_in_extension':True,'status':'pass'}
if __name__=='__main__':
    out={'asymmetric_cases':[]}
    for d,e,m in [(1,1,1),(1,3,1),(1,3,2),(2,2,1),(2,3,1),(3,1,2)]:
        z=identities(d,e,m);out['asymmetric_cases'].append(z);print(json.dumps(z),flush=True)
    out['nonreduced_witness']=nonreduced_witness();print(json.dumps(out['nonreduced_witness']),flush=True)
    (ROOT/'lagrangian_new_round2_ordinary_results.json').write_text(json.dumps(out,indent=2));print('TARGETED ORDINARY CHECKS PASSED',flush=True)
