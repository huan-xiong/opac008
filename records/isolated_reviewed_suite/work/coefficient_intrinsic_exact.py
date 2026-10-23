"""Exact checks of the two intrinsic coefficient-algebra obstructions.
No intersection data or predicted basis matrix is used to form multiplication.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import sympy as s
ROOT=Path(__file__).resolve().parent
start=datetime.now(timezone.utc).isoformat()
x,y=s.symbols('x y')
basis=[s.Integer(1),x,y,x**2,x*y,y**2]
relations={'C':[3*x**3-8*x*y,x**2*y-2*y**2],
           'H':[x**3-2*x*y,x**2*y-y**2]}
records={}
for name,rels in relations.items():
    gb=s.groebner(rels,x,y,order='lex',domain=s.QQ)
    def vector(f):
        nf=gb.reduce(s.expand(f))[1]
        p=s.Poly(nf,x,y)
        result=s.Matrix([p.coeff_monomial(b) for b in basis])
        assert s.expand(nf-sum(result[i]*basis[i] for i in range(6)))==0
        return result
    multiplications=[s.Matrix.hstack(*(vector(a*b) for b in basis)) for a in basis]
    Mx,My=multiplications[1],multiplications[2]
    mspace=s.eye(6)[:,1:]
    m2=s.Matrix.hstack(*(T*mspace for T in multiplications[1:]))
    m2=s.Matrix.hstack(*m2.columnspace())
    m3=s.Matrix.hstack(*(T*m2 for T in multiplications[1:]))
    m3=s.Matrix.hstack(*m3.columnspace())
    desired=s.Matrix.hstack(vector(x*y),vector(y**2))
    assert m3.rank()==2 and s.Matrix.hstack(m3,desired).rank()==2
    soc=s.Matrix.vstack(Mx,My).nullspace()
    assert len(soc)==1 and s.Matrix.hstack(soc[0],vector(y**2)).rank()==1
    I=s.Matrix.vstack(*(sum((m3[i,j]*multiplications[i] for i in range(6)),s.zeros(6))
                       for j in range(m3.cols))).nullspace()
    assert len(I)==4
    expectedI=s.Matrix.hstack(*(vector(f) for f in [y,x**2,x*y,y**2]))
    assert s.Matrix.hstack(*I,expectedI).rank()==4
    # Compute the intrinsic pair from multiplication in the chosen quotient basis.
    pair=[]
    for a in [y,x**2]:
        row=[]
        for b in [y,x**2]:
            v=vector(a*b)
            assert all(v[i]==0 for i in range(5))
            row.append(v[5])
        pair.append(row)
    pair=s.Matrix(pair)
    assert pair==(s.Matrix([[1,2],[2,s.Rational(16,3)]]) if name=='C' else s.Matrix([[1,1],[1,2]]))
    records[name]={'dimension':6,'m2_dimension':m2.cols,'m3_dimension':m3.cols,
                   'socle_dimension':len(soc),'ann_m3_dimension':len(I),
                   'pairing':[[str(a) for a in row] for row in pair.tolist()],
                   'determinant':str(pair.det()),'status':'PASS'}

# The rational field restriction is essential: exhibit a graded isomorphism over Q(sqrt(3)).
t=s.sqrt(3)
gbH=s.groebner(relations['H'],x,y,order='lex',extension=t)
phi_y=-t*y/4+(3+t)*x**2/8
assert all(gbH.reduce(s.expand(f.subs(y,phi_y)))[1]==0 for f in relations['C'])
assert -t/4!=0
records['quadratic_extension']={'field':'Q(sqrt(3))','map_x':'x',
    'map_y':'-sqrt(3)*y/4 + (3+sqrt(3))*x^2/8',
    'relations_vanish':'PASS','surjective':'nonzero y coefficient, x unchanged',
    'scope':'Additional boundary check: no all-field nonisomorphism claim for the Lagrangian example.'}

# The ordinary characteristic-zero obstruction: independent symbolic reduction and resultant.
g,b,a,c=s.symbols('g b alpha beta')
gbC=s.groebner([5*g**3+2*g**2*b+g*b**2,g**2*b+2*g*b**2+5*b**3],g,b,order='lex',domain=s.QQ)
F=-2*a**3+15*a**2*c-c**3
G=-a**3+15*a*c**2-2*c**3
# Work coefficientwise in alpha,beta so the reduction field remains Q.
for i in range(4):
    scalar=s.expand(F/5).coeff(a,3-i).coeff(c,i)
    scalar2=s.expand(G/5).coeff(a,3-i).coeff(c,i)
    expected=scalar*g**2*b+scalar2*g*b**2
    actual=s.binomial(3,i)*g**(3-i)*b**i
    assert gbC.reduce(s.expand(actual-expected))[1]==0
v1=gbC.reduce(g**2*b)[1]
v2=gbC.reduce(g*b**2)[1]
terms=list(set(s.Poly(v1,g,b).monoms()+s.Poly(v2,g,b).monoms()))
rank=s.Matrix([[s.Poly(v1,g,b).coeff_monomial(m),s.Poly(v2,g,b).coeff_monomial(m)] for m in terms]).rank()
assert rank==2
z=s.symbols('z')
resultant=s.resultant(F.subs({a:z,c:1}),G.subs({a:z,c:1}),z)
assert resultant!=0
assert F.subs({a:1,c:0})==-2 and G.subs({a:1,c:0})==-1
assert s.expand(F-G+(a-c)*(a**2-14*a*c+c**2))==0
assert s.expand(F+G+3*(a+c)*(a**2-6*a*c+c**2))==0
records['ordinary']={'cube_coordinate_identities':'PASS','degree3_coordinate_rank':rank,
    'affine_resultant':str(resultant),'point_at_infinity':'not a zero',
    'field_scope':'Every characteristic-zero extension of Q; resultant is a nonzero rational integer.',
    'ungraded_scope':'Local algebra isomorphism induces an isomorphism of gr_m; both original algebras are standard graded.'}

out={'status':'PASS','started_at':start,'ended_at':datetime.now(timezone.utc).isoformat(),
     'sympy':s.__version__,'records':records}
(ROOT/'coefficient_intrinsic_exact_results.json').write_text(json.dumps(out,indent=2),encoding='utf8')
print(json.dumps(out,indent=2),flush=True)
