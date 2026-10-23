import sympy as S
from itertools import product
from math import comb
import json,time

def expand(x):return S.expand(x)
def coeffs(d,e,m):
    gs=S.symbols('g1:'+str(d+1));bs=S.symbols('b1:'+str(e+1));xs=gs+bs;g=(S.Integer(1),)+gs;b=(S.Integer(1),)+bs;h=d+e
    q=[sum(g[i]*b[n-i] for i in range(max(0,n-e),min(d,n)+1)) for n in range(h+1)]
    D=[S.Integer(1)]
    for n in range(1,m+h+1):D.append(expand(-sum((n-S.Rational(i,2))*q[i]*D[n-i] for i in range(1,min(h,n)+1))/n))
    U=[S.Integer(0)];V=[S.Integer(0)]
    for n in range(1,m+h+1):
        U.append(expand(sum(i*g[i]*D[n-i]-g[i]*U[n-i] for i in range(1,min(d,n)+1))))
        V.append(expand(sum(i*b[i]*D[n-i]-b[i]*V[n-i] for i in range(1,min(e,n)+1))))
    ur=[None];vr=[None];dr=[None]
    for t in range(1,h+1):
        ur.append(S.eye(h)[t-1,:] if t<=d else S.zeros(1,h));vr.append(S.eye(h)[d+t-1,:] if t<=e else S.zeros(1,h))
        if t>d:
            for i in range(1,d+1):ur[t]+=g[i]*(-ur[t-i]+i*dr[t-i])
        if t>e:
            for i in range(1,e+1):vr[t]+=b[i]*(-vr[t-i]+i*dr[t-i])
        dr.append((-(ur[t]+vr[t])/(2*(m+t))).applyfunc(expand))
    M=S.Matrix.vstack(*dr[1:]);f=S.Matrix(D[m+1:]);j=S.Matrix(U[m+1:m+d+1]+V[m+1:m+e+1]);w=S.Symbol('w')
    C=S.Matrix(h,h,lambda r,j: (g[d+r+1-(j+1)] if 0<=d+r+1-(j+1)<=d else 0) if r<e else (b[e+(r-e)+1-(j+1)] if 0<=e+(r-e)+1-(j+1)<=e else 0))
    rho=S.expand((-1)**(d*e)*C.det(method='domain-ge'))
    assert all(expand(x)==0 for x in f-M*j)
    det=M.det(method='domain-ge');ratio=S.cancel(det/rho)
    expected=S.Rational(1,2**h*product_int(range(m+1,m+h+1)))
    assert ratio in [expected,-expected],(d,e,m,ratio,expected)
    return xs,list(f),list(j),rho,ratio

def product_int(xs):
    a=1
    for x in xs:a*=x
    return a

def standard_monomials(G,h):
    lms=[p.LM(order=G.order).exponents for p in G.polys]
    bounds=[min(lm[i] for lm in lms if lm[i]>0 and all(lm[j]==0 for j in range(h) if j!=i)) for i in range(h)]
    return [es for es in product(*(range(b) for b in bounds)) if not any(all(e>=a for e,a in zip(es,lm)) for lm in lms)]

def colon_check(d,e,m):
    st=time.time();xs,I,J,rho,ratio=coeffs(d,e,m);h=d+e;weights=tuple(range(1,d+1))+tuple(range(1,e+1));top=m*h+d*e
    G=S.groebner(I,*xs,order='grevlex');H=S.groebner(J,*xs,order='grevlex')
    assert all(G.reduce(expand(rho*x))[1]==0 for x in J)
    assert all(H.reduce(x)[1]==0 for x in I)
    mons=standard_monomials(G,h)
    assert max(sum(e*w for e,w in zip(es,weights)) for es in mons)==top
    idx={es:i for i,es in enumerate(mons)}
    cols=[]
    for es in mons:
        r=G.reduce(expand(rho*product_int(x**a for x,a in zip(xs,es))))[1]
        col=[0]*len(mons)
        for powers,c in S.Poly(r,*xs).terms():
            if c:col[idx[powers]]=c
        cols.append(col)
    rank=S.Matrix.hstack(*(S.Matrix(c) for c in cols)).rank()
    jmons=standard_monomials(H,h)
    assert len(jmons)==rank,('DIRECT COLON DIMENSION MISMATCH',d,e,m,len(jmons),rank)
    expected=comb(m+d,d)*comb(m+e,e)
    assert rank==expected,(d,e,m,rank,expected)
    expected_len=product_int(range(m+1,m+h+1))//(product_int(range(1,d+1))*product_int(range(1,e+1)))
    assert len(mons)==expected_len
    print(json.dumps({'d':d,'e':e,'m':m,'det_ratio':str(ratio),'parameter_length':len(mons),'rho_multiplication_rank':rank,'seconds':time.time()-st,'status':'exact polynomial identities and colon PASS'}),flush=True)

if __name__=='__main__':
    for case in [(1,1,1),(1,1,2),(1,2,1),(1,2,2),(2,2,1),(2,2,2),(1,3,1)]:colon_check(*case)
