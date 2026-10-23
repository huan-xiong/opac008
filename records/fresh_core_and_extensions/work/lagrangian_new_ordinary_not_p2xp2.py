"""New targeted exact check: ordinary m=2, d=e=1 colon algebra.
Direct binomial-series expansion and Groebner basis over QQ; no audit helpers.
"""
import json
from pathlib import Path
import sympy as s
g,b,z,a,c,q=s.symbols('g b z alpha beta q')
# Coefficients are computed directly from the defining generating functions.
D=s.series((1+g*z)**(-s.Rational(1,2))*(1+b*z)**(-s.Rational(1,2)),z,0,4).removeO().expand()
U=s.series(z*g*D/(1+g*z),z,0,4).removeO().expand().coeff(z,3)
V=s.series(z*b*D/(1+b*z),z,0,4).removeO().expand().coeff(z,3)
r1=5*g**3+2*g*g*b+g*b*b;r2=g*g*b+2*g*b*b+5*b**3
assert s.expand(U-s.Rational(3,8)*r1)==0
assert s.expand(V-s.Rational(3,8)*r2)==0
GB=s.groebner([r1,r2],g,b,domain=s.QQ,order='grevlex')
lms=[p.LM(order=GB.order).exponents for p in GB.polys]
# BFS produces every standard monomial, with no a priori top cutoff.
std={(0,0)};front=[(0,0)]
while front:
    ii,jj=front.pop()
    for ex in [(ii+1,jj),(ii,jj+1)]:
        if ex in std or any(all(ex[k]>=lm[k] for k in range(2)) for lm in lms):continue
        std.add(ex);front.append(ex)
hilbert=[sum(sum(ex)==n for ex in std) for n in range(max(map(sum,std))+1)]
assert hilbert==[1,2,3,2,1]
degree3=sorted(ex for ex in std if sum(ex)==3)
base_columns=[s.Poly(GB.reduce(v)[1],g,b) for v in [g*g*b,g*b*b]]
assert s.Matrix([[p.coeff_monomial(g**ex[0]*b**ex[1]) for p in base_columns] for ex in degree3]).rank()==2
F=-2*a**3+15*a*a*c-c**3;H=-a**3+15*a*c*c-2*c**3
cube=s.expand((a*g+c*b)**3-F*g*g*b/5-H*g*b*b/5)
assert all(GB.reduce(co)[1]==0 for co in s.Poly(cube,a,c).coeffs())
assert s.expand(F-H+(a-c)*(a*a-14*a*c+c*c))==0
assert s.expand(F+H+3*(a+c)*(a*a-6*a*c+c*c))==0
# A projective zero would lie in one of the two standard affine charts.
assert s.groebner([F.subs(c,1),H.subs(c,1)],a,domain=s.QQ).is_zero_dimensional is False
assert list(s.groebner([F.subs(c,1),H.subs(c,1)],a,domain=s.QQ))==[1]
assert list(s.groebner([F.subs(a,1),H.subs(a,1)],c,domain=s.QQ))==[1]
out={'U3':str(s.factor(U)),'V3':str(s.factor(V)),'groebner_basis':[str(p.as_expr()) for p in GB.polys],'hilbert_actual_standard_monomials':hilbert,'degree_three_basis':['g^2*b','g*b^2'],'five_times_cube_coefficients':[str(F),str(H)],'projective_linear_cube_zero_locus':'empty over every characteristic-zero extension','status':'pass'}
Path(__file__).with_name('lagrangian_new_ordinary_not_p2xp2_results.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
