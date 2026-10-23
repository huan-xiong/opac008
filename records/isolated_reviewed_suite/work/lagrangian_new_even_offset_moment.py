"""Independent exact checks of the new even-offset moment argument.
No ordinary_new source/report is imported. 'actual' entries are original-ring
Schur top coefficients computed from the two-root complete polynomials.
"""
from math import comb,prod
from pathlib import Path
import json

def conv(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out

def actual_top(ell,k):
 # Top Schur coefficient in h1^(2ell-2k)*h2^k, roots (1,z).
 f=[comb(2*ell-2*k,j) for j in range(2*ell-2*k+1)]
 for _ in range(k):f=conv(f,[1,1,1])
 return f[ell]-f[ell-1]

def cat(n):return comb(2*n,n)//(n+1)
def gamma(ell,k):return sum((-1)**i*comb(k,i)*cat(ell-i) for i in range(k+1))

def det(mat):
 a=[row[:] for row in mat];n=len(a)
 if not n:return 1
 prev=1;sign=1
 for k in range(n-1):
  if a[k][k]==0:
   rr=next((j for j in range(k+1,n) if a[j][k]),None)
   if rr is None:return 0
   a[k],a[rr]=a[rr],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    num=a[i][j]*pivot-a[i][k]*a[k][j]
    assert num%prev==0
    a[i][j]=num//prev
  for i in range(k+1,n):a[i][k]=0
  prev=pivot
 return sign*a[-1][-1]

def hs(values,maxn):
 out=[1]+[0]*maxn
 for x in values:
  for n in range(1,maxn+1):out[n]+=x*out[n-1]
 return out

def rectangular_schur(values,d):
 p=len(values)-1;h=hs(values,d+p)
 return det([[h[d-i+j] if d-i+j>=0 else 0 for j in range(p)] for i in range(p)])

cases=[];entries=0
for ell in range(2,15):
 moments=[actual_top(ell,k) for k in range(ell+1)]
 for k,v in enumerate(moments):assert v==gamma(ell,k),(ell,k,v);entries+=1
 for d in range(0,ell+1,2):
  p=(ell-d)//2;a=[0]+list(range(d+1,d+p+1))
  assert len(a)==p+1 and max(a)+p<=ell
  matrix=[[moments[ai+j] for j in range(p+1)] for ai in a]
  value=det(matrix);assert value>0,(ell,d,value)
  cases.append({'ell':ell,'d':d,'weight':ell+d,'p':p,'actual_pairing_determinant':value})
rect_checks=[]
for p in range(4):
 for d in [0,2,4]:
  for z in [list(range(-2,-2+p+1)),[0]*(p+1),[0]+list(range(1,p+1)),[-1]*(p+1),list(range(1,p+2))]:
   lhs=rectangular_schur(z,d)
   cof=[prod(z[j] for j in range(p+1) if j!=i) for i in range(p+1)]
   rhs=hs(cof,d)[d]
   assert lhs==rhs and lhs>=0,(p,d,z,lhs,rhs)
   if all(z):assert lhs>0
   rect_checks.append({'p':p,'d':d,'z':z,'value':lhs})

# The known odd-offset failure checks that the matrix really detects the old
# candidate dependence, rather than a different positive Gram problem.
odd_matrix=[[actual_top(5,ai+j) for j in range(2)] for ai in [0,4]]
assert odd_matrix==[[42,28],[9,6]] and det(odd_matrix)==0
assert all(3*odd_matrix[0][j]==14*odd_matrix[1][j] for j in range(2))

out={'actual_top_entries_compared_to_Catalan_moments':entries,'positive_pairing_cases':cases,'rectangle_zero_negative_endpoint_cases':len(rect_checks),'known_odd_offset_negative_control':{'ell':5,'d':3,'matrix':odd_matrix,'determinant':0},'status':'pass'}
Path(__file__).with_name('lagrangian_new_even_offset_moment_results.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'actual_top_entries':entries,'positive_pairing_cases':len(cases),'rectangle_checks':len(rect_checks),'p0_cases':sum(c['p']==0 for c in cases),'status':'pass'},indent=2),flush=True)
