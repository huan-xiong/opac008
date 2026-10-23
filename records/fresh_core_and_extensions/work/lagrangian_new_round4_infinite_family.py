"""Independent fourth-round symbolic checks for the infinite two-row family.
No ordinary_new program or report imported/read. Proof is in the review note.
"""
from fractions import Fraction
from pathlib import Path
import json
import sympy as s
r,t=s.symbols('r t')
A=1-2*r*t+s.Rational(3,2)*r*(r-1)*t**2+s.Rational(1,2)*r*(r-1)*(3-r)*t**3
B=1-r*t-r*t**2+s.Rational(1,2)*r*(r-1)*(r+2)*t**3
for X,x0,x1 in [(A,1,1-2*t),(B,1,1-t-t*t)]:
 assert s.expand(X.subs(r,0)-x0)==0 and s.expand(X.subs(r,1)-x1)==0
 residual=s.expand(X.subs(r,r+2)-(2-3*t)*X.subs(r,r+1)+(1-t)**3*X)
 assert all(residual.coeff(t,j)==0 for j in range(4))
M=s.Matrix([[-2*r,1,-r],[s.Rational(3,2)*r*(r-1),-2*r,-r],[s.Rational(1,2)*r*(r-1)*(3-r),s.Rational(3,2)*r*(r-1),s.Rational(1,2)*r*(r-1)*(r+2)]])
k=s.Matrix([2*(2*r+1),r*(3*r+1),-(5*r+3)])
assert all(s.expand(x)==0 for x in M*k)
minor=s.factor(M[:2,:2].det());assert minor==r*(5*r+3)/2
assert s.expand(k[0]+k[2]+r+1)==0

def transpose(lam):return tuple(sum(row>=j for row in lam) for j in range(1,lam[0]+1))
def hooks(lam):
 col=transpose(lam)
 return [[lam[i]-j+col[j-1]-i for j in range(1,lam[i]+1)] for i in range(len(lam))]
def bounded(lam,k):return tuple(sum(h<=k for h in row) for row in hooks(lam))
def add(p,q):
 n=max(len(p),len(q));return [(p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0) for i in range(n)]
def scalar(p,c):return [c*x for x in p]
def shift(p):return [0]+p
H=[[1],[1]]
for n in range(2,27):H.append(add(add(H[-1],scalar(H[-2],-1)),shift(H[-2])))
checks=[]
for rr in range(2,9):
 ell=3*rr-1;N=3*rr+2;F=N//2
 core=(ell+3,3)
 assert all(ell+1 not in row for row in hooks(core))
 assert bounded(core,ell)==(ell,3)
 assert bounded(transpose(core),ell)==(1,)*N
 for a in range(4,F+1):
  nu=(N-a,a);assert max(map(max,hooks(nu)))<=ell
  assert bounded(transpose(nu),ell)==(2,)*a+(1,)*(N-2*a)
 P=add(add(scalar(H[3*rr],2*(2*rr+1)),scalar(shift(H[3*rr]),rr*(3*rr+1))),scalar(H[3*rr+1],-(5*rr+3)))
 assert P[0]==(-1)**(rr+1)*(rr+1) and P[1:4]==[0,0,0]
 assert all(x==0 for x in P[F+1:])
 count=1+len(range(4,F+1));assert count==F-2
 checks.append({'r':rr,'ell':ell,'weight':N,'candidate_count':count,'certificate_coefficients':P[:F+1],'full_piece_dimension':F-2,'candidate_span_dimension':F-3,'status':'pass'})

# Post-review addition, 2026-09-14 06:40 UTC or later.
# Manuscript SHA256 f3f8c99352e342c28179f0d53e840ad6550e101be09bb40d8caf78944b60ff1c.
# This block was NOT part of the initial round-4 run.
tau=s.symbols('tau')
L=s.Matrix([[r*(r-1)/2,r-1,1]])
assert all(s.expand(x)==0 for x in L*M)
h_value=r*(r-1)/2
schur_value=s.factor((L*s.Matrix([1,-(tau-1),(tau-1)*(tau-2)/2]))[0])
assert s.expand(schur_value-(r-tau)*(r-tau+1)/2)==0
# Actual r=2 ring and Schur coordinates, independently reconstructed.
u,v=s.symbols('u v');hh=[s.Integer(1),u]
for n in range(2,8):hh.append(s.expand(u*hh[-1]-(u*u-v)*hh[-2]))
s53=s.expand(hh[5]*hh[3]-hh[6]*hh[2])
s44=s.expand(hh[4]*hh[4]-hh[5]*hh[3])
new_level2=s.expand(v*(u*u-v)**3)
assert s.expand(new_level2-s53)==0
GB=s.groebner([hh[6],hh[7]],u,v,domain=s.QQ,order='grevlex')
assert GB.reduce(s.expand(v**4-6*s53-3*s44))[1]==0
assert GB.reduce(s.expand(u**6*v-19*s53-9*s44))[1]==0
normal=[s.Poly(GB.reduce(x)[1],u,v) for x in [s53,s44]]
assert s.Matrix([[p.coeff_monomial(u**(8-2*a)*v**a) for p in normal] for a in range(5)]).rank()==2
r2_h=s.Matrix([[19,6],[9,3]])
r2_level2=s.Matrix([[1,6],[0,3]])
assert r2_h.det()==3 and r2_level2.det()==3
replacement={'added_after_initial_round4':True,'manuscript_sha256':'f3f8c99352e342c28179f0d53e840ad6550e101be09bb40d8caf78944b60ff1c','symbolic_left_kernel':'pass','h_new_column_value':str(h_value),'level2_new_column_value':str(schur_value),'r2_schur_basis':['s_(5,3)','s_(4,4)'],'r2_complete_columns':[list(r2_h[:,j]) for j in range(2)],'r2_level2_columns':[list(r2_level2[:,j]) for j in range(2)],'r2_two_determinants':[3,3],'r2_schur_classes_independent_by_QQ_normal_forms':True,'status':'pass'}
# Convert SymPy integers only for JSON serialization.
replacement['r2_complete_columns']=[[int(x) for x in col] for col in replacement['r2_complete_columns']]
replacement['r2_level2_columns']=[[int(x) for x in col] for col in replacement['r2_level2_columns']]

out={'symbolic_jet_recurrences':'pass','symbolic_missing_coordinate_matrix_kernel':'pass','symbolic_minor':str(minor),'cases':checks,'uniform_replacement_addition':replacement,'status':'pass'}
Path(__file__).with_name('lagrangian_new_round4_infinite_family_results.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2),flush=True)
