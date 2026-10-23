import sympy as s,json,sys
from math import comb
from ordinary_new_exact import invpower,multiply,mons,polynomial,direct_component,rank,gbin
from time import time

def setup(M,d):
 var=s.symbols('g1:'+str(d+1)); g=[s.S(1)]+list(var)
 q=multiply(g,[(-1)**i*a for i,a in enumerate(g)],2*d)
 D=invpower(None,q,s.Rational(-1,2),M+2*d)
 Gi=invpower(None,g,-1,M+2*d)
 U=multiply(multiply([i*a for i,a in enumerate(g)],Gi,M+2*d),D,M+2*d)
 if M%2:ids=list(range(M+1,M+2*d,2));A=[D[i] for i in ids];aw=ids
 else:ids=list(range(M+2,M+2*d,2));A=[D[i] for i in ids]+[s.expand(g[-1]*D[M])];aw=ids+[M+d]
 B=U[M+1:M+d+1];bw=list(range(M+1,M+d+1))
 rho=s.det(s.Matrix(d-1,d-1,lambda a,j:g[d+(a+1)-2*(j+1)] if 0<=d+(a+1)-2*(j+1)<=d else 0)) if d>1 else s.S(1)
 return var,g,D,U,A,aw,B,bw,s.expand(rho)

from functools import lru_cache

def conductor(M,d):
 var,g,D,U,A,aw,B,bw,rho=setup(M,d);weights=list(range(1,d+1));n=M+d-1;top=M*d+d*(d-1)//2
 gb=s.groebner(A,*var,domain=s.QQ);leading=[p.LM(order=gb.order).exponents for p in gb.polys]
 std={}
 for wt in range(top+1):
  std[wt]=[ex for ex in mons(weights,wt) if not any(all(a>=b for a,b in zip(ex,lm)) for lm in leading)]
 @lru_cache(None)
 def normal(ex):
  wt=sum(a*w for a,w in zip(ex,weights))
  if wt>top:return s.zeros(0,1)
  rem=s.Poly(gb.reduce(polynomial(ex,var))[1],*var)
  return s.Matrix([rem.coeff_monomial(polynomial(e,var)) for e in std[wt]])
 def nfpoly(f,wt):
  rem=s.Poly(gb.reduce(s.expand(f))[1],*var)
  return s.Matrix([rem.coeff_monomial(polynomial(e,var)) for e in std[wt]])
 def annihilating_rows(mat):
  # A row basis for functionals vanishing on the column space.
  vec=mat.T.nullspace()
  return s.Matrix.vstack(*[v.T for v in vec]) if vec else s.zeros(0,mat.rows)
 E=multiply(g,D,n)[1:];Bcols={};Bperp={}
 for wt in range(top+1):
  cols=[nfpoly(s.prod(E[i]**a for i,a in enumerate(ex)),wt) for ex in mons(list(range(1,n+1)),wt)]
  bm=s.Matrix.hstack(*cols) if cols else s.zeros(len(std[wt]),0)
  Bcols[wt]=bm;Bperp[wt]=annihilating_rows(bm)
 output=[]
 for wt in range(top+1):
  equations=[]
  for tw in range(top-wt+1):
   perp=Bperp[wt+tw]
   if perp.rows==0:continue
   for mul_exp in std[tw]:
    mat=s.Matrix.hstack(*[normal(tuple(a+b for a,b in zip(ex,mul_exp))) for ex in std[wt]]) if std[wt] else s.zeros(len(std[wt+tw]),0)
    equations.append(perp*mat)
  allmat=s.Matrix.vstack(*equations) if equations else s.zeros(0,len(std[wt]))
  cdim=len(std[wt])-rank(allmat)
  shift=d*(d-1)//2
  if wt>=shift:
   rc=[nfpoly(rho*polynomial(ex,var),wt) for ex in std[wt-shift]]
   rr=rank(s.Matrix.hstack(*rc)) if rc else 0
  else:rr=0
  expect=len(std[wt]) if M==1 else rr
  assert cdim==expect,('full conductor',M,d,wt,cdim,expect)
  output.append({'weight':wt,'T':len(std[wt]),'B':rank(Bcols[wt]),'conductor':cdim,'rhoT':rr})
 print('full conductor',M,d,'Tdim',sum(x['T'] for x in output),'cond_dim',sum(x['conductor'] for x in output),'rhoT_dim',sum(x['rhoT'] for x in output),flush=True)
 return {'M':M,'d':d,'data':output,'method':'Solve t*T subset B degreewise using multiplication by EVERY Groebner standard basis element of T; no a priori conductor formula used.'}

out=[]
for M,d in [(1,1),(1,2),(1,3),(1,4),(2,1),(2,2),(2,3),(3,2),(3,3),(4,2)]:
 out.append(conductor(M,d))
 with open(sys.argv[1],'w',encoding='utf8') as f:json.dump(out,f,indent=2)
print('ALL FULL CONDUCTOR CHECKS PASSED',flush=True)
