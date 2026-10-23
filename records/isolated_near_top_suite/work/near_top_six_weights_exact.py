"""Independent exact near-top verification from defining Q[u,v] relations.
Finite checks validate the new analytic proof; they do not prove a classification.
"""
import sympy as s
import json,time
from pathlib import Path
u,v=s.symbols('u v');start=time.monotonic()
def h_list(n):
    h=[s.S.One,u]
    while len(h)<=n:h.append(s.expand(u*h[-1]+(v-u*u)*h[-2]))
    return h

def coords(f,N):
    f=s.Poly(f,u,v)
    return [f.coeff_monomial(u**(N-2*a)*v**a) for a in range(N//2+1)]
def moment(j,shift=0):return sum(s.binomial(j,i)*(-1)**(j-i)*s.catalan(i+shift) for i in range(j+1))
records=[]
for k in range(2,21):
    h=h_list(k+2)
    for N in range(max(0,2*k-5),2*k+1):
        relations=[]
        for n in (k+1,k+2):
            excess=N-n
            for a in range(excess//2+1):relations.append(coords(u**(excess-2*a)*v**a*h[n],N))
        relation_matrix=s.Matrix(relations) if relations else s.zeros(0,N//2+1)
        relation_rank=relation_matrix.rank();quotient_dim=N//2+1-relation_rank
        d=N-k
        indices=range(N//2+1) if N<=k else [0]+list(range(d+1,N//2+1))
        candidates=[coords(u**(N-2*a)*v**a,N) for a in indices]
        image_rank=s.Matrix(relations+candidates).rank()-relation_rank
        codim=quotient_dim-image_rank
        expected=1 if (k,N) in ((5,8),(8,11)) else 0
        assert codim==expected
        row={'k':k,'N':N,'excess':d,'quotient_dimension':quotient_dim,'candidate_count':len(candidates),'actual_span_rank':image_rank,'codimension':codim}
        if d>=1 and d%2 and (k-d)//2==1:
            m=d+1;sigma=k-d-2
            def gamma(b):return sum(s.binomial(b,i)*(-1)**i*s.catalan(k-i) for i in range(b+1))
            actual=s.Matrix([[gamma(0),gamma(1)],[gamma(m),gamma(m+1)]]).det()
            shifted=moment(m,sigma);prev=moment(m-1,sigma)
            formula=s.Rational(m,m+2+sigma)*s.catalan(m+sigma)*(3*prev-shifted)
            assert actual==formula
            row.update({'m':m,'sigma':sigma,'pairing_determinant':str(actual),'Pearson_formula_passed':True})
        if d>=1 and d%2 and (k-d)//2==2:
            m=d+1;sigma=k-d-4
            def gamma(b):return sum(s.binomial(b,i)*(-1)**i*s.catalan(k-i) for i in range(b+1))
            actual=s.Matrix([[gamma(a+j) for j in range(3)] for a in (0,m,m+1)]).det()
            T=moment(m,sigma);ratio=s.cancel(moment(m-1,sigma)/T)
            if sigma==0:
                Q=9*m*(3*m+7)*ratio**2+3*(6*m*m+13*m-6)*ratio-(9*m*m+38*m+48)
                formula=s.cancel(T*T*s.catalan(m+1)*2*m*(m+1)*Q/((m+2)**2*(m+3)**2*(m+4)))
            else:
                Q=9*m*(2*m+7)*ratio**2+6*(2*m*m+9*m+6)*ratio-(6*m*m+37*m+60)
                formula=s.cancel(T*T*s.catalan(m+2)*3*m*(m+1)*Q/((m+3)**2*(m+4)**2*(m+5)))
            assert actual==formula
            row.update({'m':m,'sigma':sigma,'pairing_determinant':str(actual),'p2_quadratic_formula_passed':True})
        records.append(row);print(json.dumps(row),flush=True)
result={'research_only':True,'method':'exact weighted defining-ideal quotient ranks over Q[u,v], then an independent Catalan-moment/Pearson determinant check','cases':len(records),'all_passed':True,'seconds':time.monotonic()-start,'records':records}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2),encoding='utf8')
