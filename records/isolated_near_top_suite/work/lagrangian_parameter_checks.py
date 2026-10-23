"""Independent Groebner-basis tests of finite parameter images and colon ranks.
The image is generated from the actual coefficient polynomials by exact linear
algebra in S/a; colon dimensions are ranks of multiplication by rho in S/a.
"""
import sympy as s
from functools import lru_cache
from collections import defaultdict
import json, time
from lagrangian_independent import gaussian, monomials

def expected_H(n,d):
    out=defaultdict(int)
    for a in range(d+1):
        for t,c in gaussian(n,a).items():out[t+a*(a+1)//2]+=c
    return dict(out)

def coeffs(M,d):
    g=(s.Integer(1),)+s.symbols('g1:'+str(d+1)); v=g[1:]
    cutoff=M+2*d+2
    w=[sum((-1)**j*g[i]*g[j] for i in range(d+1) for j in range(d+1) if i+j==2*a) for a in range(d+1)]
    D=[s.Integer(0)]*(cutoff+1);D[0]=s.Integer(1)
    for a in range(1,cutoff//2+1):
        D[2*a]=s.expand(-sum((s.Rational(a)-s.Rational(i,2))*w[i]*D[2*(a-i)] for i in range(1,min(d,a)+1))/a)
    U=[s.Integer(0)]*(cutoff+1)
    for t in range(1,cutoff+1):
        U[t]=s.expand(sum(i*g[i]*D[t-i]-g[i]*U[t-i] for i in range(1,min(t,d)+1)))
    a=[D[j] for j in (range(M+2,M+2*d,2) if M%2==0 else range(M+1,M+2*d,2))]
    if M%2==0:a.append(s.expand(g[d]*D[M]))
    b=U[M+1:M+d+1]
    rho=s.det(s.Matrix(d-1,d-1,lambda a,j:g[d+(a+1)-2*(j+1)] if 0<=d+(a+1)-2*(j+1)<=d else 0)) if d>1 else s.Integer(1)
    E=[s.expand(sum(g[i]*D[t-i] for i in range(min(d,t)+1))) for t in range(cutoff+1)]
    return v,a,b,s.expand(rho),E


def all_standard_monomials(gb):
    """Enumerate the complete monomial basis, with no predicted degree cutoff."""
    leads=[p.LM(order=gb.order).exponents for p in gb.polys]
    d=len(gb.gens)
    assert all(any(l[i]>0 and all(l[j]==0 for j in range(d) if j!=i) for l in leads) for i in range(d)), 'leading ideal is not zero-dimensional'
    zero=(0,)*d;seen={zero};pending=[zero]
    while pending:
        a=pending.pop()
        for i in range(d):
            b=list(a);b[i]+=1;b=tuple(b)
            if b in seen or any(all(x>=y for x,y in zip(b,l)) for l in leads):continue
            seen.add(b);pending.append(b)
    return seen

def run_case(M,d):
    begin=time.time();n=M+d-1;v,a,b,rho,E=coeffs(M,d)
    gb=s.groebner(a,*v,order='grevlex',domain=s.QQ)
    leads=[p.LM(order=gb.order).exponents for p in gb.polys]
    all_mons=all_standard_monomials(gb)
    top=max(sum((i+1)*a for i,a in enumerate(e)) for e in all_mons)
    assert top==d*M+d*(d-1)//2, ('socle degree',M,d,top)
    mons={t:sorted(e for e in all_mons if sum((i+1)*a for i,a in enumerate(e))==t) for t in range(top+1)}
    @lru_cache(None)
    def nf(p):return gb.reduce(p)[1]
    def asvec(p,t):
        co=s.Poly(p,*v).as_dict();assert not (set(co)-set(mons.get(t,[]))), ('normal-form weight',M,d,t,co)
        return s.Matrix([co.get(e,0) for e in mons.get(t,[])])
    def pol(e):return s.prod(x**a for x,a in zip(v,e))
    def basis(polys,t):
        polys=[nf(s.expand(p)) for p in polys];polys=[p for p in polys if p!=0]
        if not polys:return []
        mat=s.Matrix.hstack(*(asvec(p,t) for p in polys));pivots=s.polys.matrices.DomainMatrix.from_Matrix(mat).rref()[1]
        return [polys[j] for j in pivots]
    B={0:[s.Integer(1)]};actual={0:1}
    for t in range(1,top+1):
        B[t]=basis([E[i]*p for i in range(1,min(n,t)+1) for p in B[t-i]],t)
        if B[t]:actual[t]=len(B[t])
    assert actual==expected_H(n,d),(M,d,'image',actual,expected_H(n,d))
    shift=d*(d-1)//2;colon={}
    for t in range(top-shift+1):
        ps=[nf(s.expand(rho*pol(e))) for e in mons[t]]
        rank=len(basis(ps,t+shift))
        if rank:colon[t]=rank
    assert colon==gaussian(M+d,d),(M,d,'colon',colon,gaussian(M+d,d))
    for f in b:assert nf(s.expand(rho*f))==0,(M,d,'annihilator',f)
    for i in range(n+1,len(E)):assert nf(E[i])==0
    return {'M':M,'d':d,'n':n,'parameter_dimension':sum(map(len,mons.values())),'image_dimension':sum(actual.values()),'colon_quotient_dimension':sum(colon.values()),'seconds':round(time.time()-begin,3)}

if __name__=='__main__':
    for d in range(1,5):
        for M in range(1,6):
            print(json.dumps(run_case(M,d)),flush=True)
