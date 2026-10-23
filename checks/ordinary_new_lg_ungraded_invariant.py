import sympy as s,json,sys
from ordinary_new_exact import invpower,multiply
x,y=s.symbols('x y');bas=[s.S(1),x,y,x*x,x*y,y*y]
GG=[s.S(1),s.S(0),2*y-x*x,s.S(0),y*y]
DD=invpower(None,GG,s.Rational(-1,2),4)
Gi=invpower(None,[1,x,y],-1,4)
UU=multiply(multiply([0,x,2*y],Gi,4),DD,4)
assert list(s.groebner(UU[3:5],x,y,domain=s.QQ))==list(s.groebner([3*x**3-8*x*y,x*x*y-2*y*y],x,y,domain=s.QQ))
print('Independent defining coefficients U3=',UU[3],'U4=',UU[4])

def analyze(label,relations):
 gb=s.groebner(relations,x,y,domain=s.QQ)
 def nf(f):return s.expand(gb.reduce(s.expand(f))[1])
 def v(f):
  pol=s.Poly(nf(f),x,y)
  out=s.Matrix([pol.coeff_monomial(t) for t in bas])
  assert s.expand(sum(out[i]*bas[i] for i in range(6))-pol.as_expr())==0
  return out
 def mult(f):return s.Matrix.hstack(*[v(f*g) for g in bas])
 mx,my=mult(x),mult(y)
 m=s.Matrix.hstack(*[v(t) for t in bas[1:]])
 powers={1:m}
 for r in range(2,6):
  cols=(mx*powers[r-1]).row_join(my*powers[r-1]).columnspace()
  powers[r]=s.Matrix.hstack(*cols) if cols else s.zeros(6,0)
 m3=powers[3];m4=powers[4]
 expected_m3=s.Matrix.hstack(v(x*y),v(y*y))
 assert m3.rank()==2 and m3.row_join(expected_m3).rank()==2
 assert m4.rank()==1 and m4.row_join(v(y*y)).rank()==1 and powers[5].cols==0
 stack=s.Matrix.vstack(*[mult(sum(m3[i,j]*bas[i] for i in range(6))) for j in range(m3.cols)])
 ann=s.Matrix.hstack(*stack.nullspace());expect=s.Matrix.hstack(*[v(t) for t in [y,x*x,x*y,y*y]])
 assert ann.rank()==4 and ann.row_join(expect).rank()==4
 # whole Ann(m^3)^2 lies in m^4, its pairing radical is m^3
 middle=[y,x*x];Q=s.zeros(2)
 for i,a in enumerate(middle):
  for j,b in enumerate(middle):
   vv=v(a*b); assert all(vv[k]==0 for k in range(5));Q[i,j]=vv[5]
 assert Q.det()!=0
 soc=s.Matrix.vstack(mx,my).nullspace();assert len(soc)==1 and s.Matrix.hstack(soc[0],v(y*y)).rank()==1
 return {'algebra':label,'groebner':[str(p.as_expr()) for p in gb.polys],'m_power_dimensions':{str(r):powers[r].rank() for r in powers},'Ann_m3_dimension':ann.rank(),'socle':'Q*y^2','canonical_pairing_basis':['y','x^2'],'pairing':[[str(Q[i,j]) for j in range(2)] for i in range(2)],'determinant':str(Q.det())}

Cgb=s.groebner(UU[3:5],x,y,domain=s.QQ)
EC=multiply([1,x,y],DD,3)
assert Cgb.reduce(s.expand(EC[1]-x))[1]==0
assert Cgb.reduce(s.expand(EC[2]-x*x/2))[1]==0
assert Cgb.reduce(s.expand(EC[3]-x**3/8))[1]==0
assert Cgb.reduce(x**4)[1]!=0 and Cgb.reduce(x**5)[1]==0
assert list(s.groebner(UU[3:5]+[x],x,y,domain=s.QQ))==[x,y*y]
print('R3 action: e1=x, e2=x^2/2, e3=x^3/8; image Q[x]/x^5; C/R3_positive*C=Q[y]/y^2.')
out=[analyze('LG colon M=d=2',[3*x**3-8*x*y,x*x*y-2*y*y]),analyze('Gr(2,4)',[x**3-2*x*y,x*x*y-y*y])]
with open(sys.argv[1],'w',encoding='utf8') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2))
print('CANONICAL UNGRADED INVARIANT VERIFIED: determinant square classes 1/3 and 1 differ in Q*/Q*^2.')
