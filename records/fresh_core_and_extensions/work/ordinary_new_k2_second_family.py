import sympy as s,json,sys
from functools import lru_cache
from math import comb
from time import time

@lru_cache(None)
def tri(a):
 if a==0:return (1,)
 prev=tri(a-1);out=[0]*(2*a+1)
 for j,x in enumerate(prev):
  for k in range(3):out[j+k]+=x
 return tuple(out)

def from_roots(a,b,c=0):
 # h2^a h1^b e2^c in roots X,Y: set X=1, Y=z.
 n=2*a+b+2*c;out=[0]*(n+1)
 for i,x in enumerate(tri(a)):
  for j in range(b+1):out[i+j+c]+=x*comb(b,j)
 return [out[j]-(out[j-1] if j else 0) for j in range(n//2+1)]

def pieri(vec,r):
 out={}
 for (p,q),coef in vec.items():
  for dv in range(r+1):
   if q+dv<=p:
    key=(p+r-dv,q+dv);out[key]=out.get(key,0)+coef
 return out

@lru_cache(None)
def actual_h(a,b):
 if b:return pieri(actual_h(a,b-1),1)
 if a:return pieri(actual_h(a-1,0),2)
 return {(0,0):1}

def actual_second(a,b):
 if a==0:return actual_h(0,b)
 eps=b%2;c=b//2
 # e2 multiplication adds one box to each row; partitions with >2 rows vanish.
 return {(p+c,q+c):coef for (p,q),coef in actual_h(a,eps).items()}

def hooks(part):
 return [[part[i]-j+sum(x>j for x in part)-i-1 for j in range(part[i])] for i in range(len(part))]
def transpose(part):return tuple(sum(x>=j for x in part) for j in range(1,max(part,default=0)+1))
def bounded(part,k):return tuple(sum(x<=k for x in row) for row in hooks(part) if any(x<=k for x in row))

start=time();bound=int(sys.argv[1]);count=0;hook_count=0;fail=[];stats=[]
for ell in range(2,bound+1):
 for r in range(ell+1):
  for ss in range(r+1):
   mu=tuple(x for x in (r,ss) if x)
   core=tuple(x for x in ((r,ss) if r<ell else (ell+ss,ss)) if x)
   expected=tuple([2]*ss+[1]*(r-ss)) if r<ell else tuple([1]*(ell+ss))
   assert all((ell+1) not in row for row in hooks(core))
   assert bounded(core,ell)==mu,(ell,mu,core)
   assert bounded(transpose(core),ell)==expected,(ell,mu,core,expected)
   hook_count+=1
 failed=[]
 for n in range(2*ell+1):
  t=n//2;aa=[0]+list(range(max(1,n-ell+1),t+1));rows=list(range(max(0,n-ell),t+1))
  first=[from_roots(a,n-2*a) for a in aa]
  second=[from_roots(0,n) if a==0 else from_roots(a,n%2,t-a) for a in aa]
  fm=s.Matrix([[c[j] for c in first] for j in rows]);sm=s.Matrix([[c[j] for c in second] for j in rows])
  change=s.zeros(len(aa));change[0,0]=1
  for ci,a in enumerate(aa[1:],1):
   for ri,j in enumerate(aa[1:],1):
    if j>=a:change[ri,ci]=(-1)**(j-a)*comb(t-a,j-a)
  assert change.det()==1 and sm==fm*change,(ell,n,'triangular')
  for a,col in zip(aa,second):
   pv=actual_second(a,n-2*a)
   assert all(col[j]==pv.get((n-j,j),0) for j in range(t+1)),(ell,n,a,'actual Pieri')
   count+=1
  det=sm.det(method='domain-ge')
  if det==0:
   rank=sm.rank();assert rank==fm.rank();failed.append(n)
   fail.append({'ell':ell,'weight':n,'dimension':len(rows),'rank':rank})
 stats.append({'ell':ell,'failed_weights':failed})
 if ell%10==0:print('second family completed ell',ell,'seconds',round(time()-start,2),flush=True)
# Independently compute Riordan numbers as actual h2 powers in equal rows.
rr=[from_roots(a,0)[a] for a in range(201)]
assert rr[:7]==[1,0,1,1,3,6,15]
assert all((n+1)*rr[n]==(n-1)*(2*rr[n-1]+3*rr[n-2]) for n in range(2,201))
endpoint=[a+1 for a in range(1,200) if (a+2)*rr[a+1]==3*a*rr[a]]
assert endpoint==[5],endpoint
for a in range(1,200):
 h=from_roots(a,0)
 assert h[a-1]==rr[a+1],a
res={'bound':bound,'core_hook_checks':hook_count,'second_family_actual_Pieri_column_checks':count,'triangular_span_identity_verified_all_weights':True,'failures':fail,'stats':stats,'endpoint_degenerate_ell_through_200':endpoint,'riordan_numbers_0_to_20':rr[:21],'seconds':time()-start}
with open(sys.argv[2],'w',encoding='utf8') as f:json.dump(res,f,indent=2)
print(json.dumps({k:v for k,v in res.items() if k not in ['stats']},indent=2),flush=True)
