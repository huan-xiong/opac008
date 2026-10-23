"""Exact negative controls for two intentionally false parity extensions."""
import sympy as s
import json
g1,g2=s.symbols('g1 g2')
D2=(g1*g1-2*g2)/2
D4=s.Rational(3,8)*(g1*g1-2*g2)**2-g2*g2/2
G=s.groebner([D4],g1,g2,order='grevlex')
missing_top=G.reduce(s.expand(g2*D2))[1]
assert missing_top != 0
x,y,z=s.symbols('x y z')
R=s.groebner([x*x-2*y,y*y-2*x*z,z*z],x,y,z,order='grevlex')
a=R.reduce(x**3)[1]
b=R.reduce(2*x**3-4*z)[1]
vecs=[]
terms=sorted(set(s.Poly(a,x,y,z).monoms()+s.Poly(b,x,y,z).monoms()))
for p in [a,b]:
    dic=dict(s.Poly(p,x,y,z).terms())
    vecs.append([dic.get(e,0) for e in terms])
assert s.Matrix.hstack(*(s.Matrix(v) for v in vecs)).rank()==2
print(json.dumps({'omitted_top_relation_remainder':str(missing_top),'even_containment_rank':2,'status':'PASS'}))
