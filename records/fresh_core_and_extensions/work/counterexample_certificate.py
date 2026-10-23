"""Exact polynomial certificate and Schur coordinates, independent of filtration proofs."""
import json
import sympy as s
x,y=s.symbols('x y')
h=[s.Integer(1),x]
for j in range(2,9): h.append(s.expand(x*h[-1]-y*h[-2]))
certificate=s.expand(3*x**8-14*h[2]**4-(14*y-24*x*x)*h[6]-13*x*h[7])
assert certificate==0
G=s.groebner([h[6],h[7]],x,y,order='grevlex')
S53=s.expand(h[5]*h[3]-h[6]*h[2])
S44=s.expand(h[4]*h[4]-h[5]*h[3])
columns=[[28,14],[6,3],[19,9]]
for p,(a,b) in zip([h[1]**8,h[2]**4,h[2]*h[1]**6],columns):
    assert G.reduce(s.expand(p-a*S53-b*S44))[1]==0
assert s.Matrix([[28,19],[14,9]]).det()==-14
def quotient_rank(polys):
    normals=[G.reduce(s.expand(p))[1] for p in polys]
    supports=sorted(set().union(*(set(s.Poly(p,x,y).monoms()) for p in normals)))
    cols=[]
    for p in normals:
        dic=dict(s.Poly(p,x,y).terms())
        cols.append(s.Matrix([dic.get(e,0) for e in supports]))
    return s.Matrix.hstack(*cols).rank()
assert quotient_rank([S53,S44])==2
assert quotient_rank([h[1]**8,h[2]**4])==1
assert quotient_rank([h[1]**8,h[2]*h[1]**6])==2
print(json.dumps({'schur_basis_rank':2,'candidate_rank':1,'replacement_rank':2,'certificate_zero':True,'schur_columns':columns,'replacement_determinant':-14,'sympy':s.__version__}))
