import sympy as s
from itertools import product
from math import comb
from time import time

def invpower(poly, coeff, alpha, nmax):
    out=[s.S(1)]
    # q D' = alpha q' D
    for n in range(1,nmax+1):
        out.append(s.expand(sum((alpha*i-(n-i))*coeff[i]*out[n-i] for i in range(1,min(n,len(coeff)-1)+1))/n))
    return out

def multiply(a,b,nmax):
    return [s.expand(sum(a[j]*b[n-j] for j in range(max(0,n-len(b)+1),min(n,len(a)-1)+1))) for n in range(nmax+1)]

def mons(weights,n,prefix=()):
    if not weights:
        if n==0: yield prefix
        return
    for j in range(n//weights[0]+1):
        yield from mons(weights[1:],n-j*weights[0],prefix+(j,))

def polynomial(exp,vars):
    return s.prod(v**a for a,v in zip(exp,vars))

def direct_component(rel, relweights,vars,weights,n):
    basis=list(mons(weights,n)); loc={x:i for i,x in enumerate(basis)}
    cols=[]
    for f,w in zip(rel,relweights):
        for t in mons(weights,n-w):
            p=s.Poly(f*polynomial(t,vars),*vars)
            col=[s.S(0)]*len(basis)
            for exp,c in p.terms(): col[loc[exp]]=c
            cols.append(col)
    return basis, s.Matrix(len(basis), len(cols),lambda i,j:cols[j][i])

def rank_polys(polys,vars,weights,n):
    basis=list(mons(weights,n)); loc={x:i for i,x in enumerate(basis)}
    return s.Matrix(len(basis),len(polys),lambda i,j:s.S(0)) if not polys else s.polys.matrices.DomainMatrix.from_Matrix(s.Matrix([[s.Poly(f,*vars).coeff_monomial(polynomial(exp,vars)) for f in polys] for exp in basis])).rank()

def rank(M):
    return s.polys.matrices.DomainMatrix.from_Matrix(M).rank()

def gbin(n,k):
    if k<0 or k>n:return [0]
    if not k or k==n:return [1]
    a=gbin(n-1,k-1); b=[0]*k+gbin(n-1,k)
    return [ (a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b))) ]

def convolution(a,b):
    ans=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):ans[i+j]+=x*y
    return ans

def run(d,e,m):
    start=time(); h=d+e; k=m+d; ell=m+e
    var=s.symbols('g1:'+str(d+1))+s.symbols('b1:'+str(e+1))
    g=(s.S(1),)+var[:d]; b=(s.S(1),)+var[d:]
    weights=list(range(1,d+1))+list(range(1,e+1)); top=m*h+d*e
    q=multiply(g,b,h); D=invpower(None,q,s.Rational(-1,2),m+h)
    # logarithmic derivatives independently use ordinary power inverse
    Ginv=invpower(None,g,-1,m+h); Binv=invpower(None,b,-1,m+h)
    U=multiply(multiply([i*x for i,x in enumerate(g)],Ginv,m+h),D,m+h)
    V=multiply(multiply([i*x for i,x in enumerate(b)],Binv,m+h),D,m+h)
    I=D[m+1:m+h+1]; J=U[m+1:m+d+1]+V[m+1:m+e+1]
    iw=list(range(m+1,m+h+1)); jw=list(range(m+1,m+d+1))+list(range(m+1,m+e+1))
    w=s.symbols('w'); cg=sum(g[i]*w**(d-i) for i in range(d+1)); cb=sum(b[i]*w**(e-i) for i in range(e+1))
    rho=s.resultant(cg,cb,w)
    E=multiply(g,D,k); Ec=E[1:]
    expected=[0]*(k*ell+1)
    for r in range(d+1):
        pol=convolution(gbin(k,r),gbin(ell,k-r)); shift=r*(r+ell-k)
        for n,x in enumerate(pol):expected[n+shift]+=x
    data=[]
    for n in range(top+1):
        basis,matI=direct_component(I,iw,var,weights,n); ri=rank(matI)
        # image by all monomials in R generators, quotient by I, no assumed polar equality
        ecmon=[s.expand(s.prod(Ec[i]**a for i,a in enumerate(exp))) for exp in mons(list(range(1,k+1)),n)]
        ecos=s.Matrix([[s.Poly(f,*var).coeff_monomial(polynomial(exp,var)) for f in ecmon] for exp in basis])
        img=rank(matI.row_join(ecos))-ri
        assert img==expected[n],('image',d,e,m,n,img,expected[n])
        if n>=d*e:
            baselow,matJ=direct_component(J,jw,var,weights,n-d*e)
            rhopol=[s.expand(rho*polynomial(exp,var)) for exp in baselow]
            rhomat=s.Matrix([[s.Poly(f,*var).coeff_monomial(polynomial(exp,var)) for f in rhopol] for exp in basis])
            rr=rank(matI.row_join(rhomat))-ri
            cj=len(baselow)-rank(matJ)
            assert rr==cj,('colon dimension',d,e,m,n,rr,cj)
            # whole rhoJ in I is checked, dimensions then establish colon equality degreewise
            if matJ.cols:
                assert rank(matI.row_join(rhomat*matJ))==ri,('rhoJ containment',d,e,m,n)
        data.append({'n':n,'Tdim':len(basis)-ri,'image':img})
        print('degree',d,e,m,n,'passed',flush=True)
    return {'d':d,'e':e,'m':m,'dimensions':data,'seconds':round(time()-start,3),'method':'Exact QQ homogeneous Macaulay ranks; independently enumerate all generator monomial images; direct resultant multiplication and J containment'}

if __name__=='__main__':
 import json,sys
 records=[]
 for a in [(1,1,1),(1,2,1),(1,1,2),(1,2,2),(2,2,1)]:
    records.append(run(*a))
    with open(sys.argv[1],'w',encoding='utf8') as f:json.dump(records,f,indent=2)
 print('ALL CHECKS PASSED',flush=True)
