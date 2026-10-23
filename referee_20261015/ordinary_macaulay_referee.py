import itertools, json, pathlib, sympy as sp, hashlib, time, sys
from functools import lru_cache
OUT=pathlib.Path(__file__).with_suffix('.json')
SOURCE=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else None
start=time.time()

@lru_cache(None)
def mons(weights,n):
    if n<0: return ()
    if not weights: return ((),) if n==0 else ()
    return tuple((a,)+b for a in range(n//weights[0]+1) for b in mons(weights[1:],n-a*weights[0]))
def mul(vars,exp):
    return sp.prod(v**a for v,a in zip(vars,exp))
def vec(poly,vars,weights,n):
    basis=mons(weights,n); pos={a:i for i,a in enumerate(basis)}
    result=sp.zeros(len(basis),1)
    if poly==0:return result
    if not vars:
        result[0]=poly;return result
    for exp,coef in sp.Poly(sp.expand(poly),*vars).terms():
        if coef:
            assert sum(w*a for w,a in zip(weights,exp))==n,(poly,n)
            result[pos[exp]]=coef
    return result
class Ring:
    def __init__(self,vars,weights,gens):
        self.vars=tuple(vars);self.weights=tuple(weights);self.gens=tuple((sp.expand(f),int(w)) for f,w in gens if f!=0)
        self.cache={}; self.qcache={}
    def mat(self,n):
        if n not in self.cache:
            cols=[vec(f*mul(self.vars,a),self.vars,self.weights,n) for f,w in self.gens for a in mons(self.weights,n-w)]
            self.cache[n]=sp.Matrix.hstack(*cols) if cols else sp.zeros(len(mons(self.weights,n)),0)
        return self.cache[n]
    def quot(self,n):
        if n not in self.qcache:
            ns=self.mat(n).T.nullspace()
            self.qcache[n]=sp.Matrix.vstack(*(v.T for v in ns)) if ns else sp.zeros(0,len(mons(self.weights,n)))
        return self.qcache[n]
    def dim(self,n):return self.quot(n).rows

def inverse_coeff(coeff,N):
    out=[sp.Integer(1)]
    for n in range(1,N+1):out.append(sp.expand(-sum(coeff[i]*out[n-i] for i in range(1,min(n,len(coeff)-1)+1))))
    return out

def make_parameter(d,e,m):
    gv=sp.symbols('g1:'+str(d+1));bv=sp.symbols('b1:'+str(e+1));vars=gv+bv;weights=tuple(range(1,d+1))+tuple(range(1,e+1));h=d+e
    G=[sp.Integer(1)]+list(gv);B=[sp.Integer(1)]+list(bv);N=m+h
    q=[sp.expand(sum(G[i]*B[n-i] for i in range(len(G)) if 0<=n-i<len(B))) for n in range(h+1)]
    # Solve [z^n](GB)D^2=0 successively. No published window recurrence or transition matrix is used.
    D=[sp.Integer(1)]
    for n in range(1,N+1):
        residual=sum(q[i]*D[a]*D[n-i-a] for i in range(min(h,n)+1) for a in range(n-i+1) if a<len(D) and n-i-a<len(D))
        D.append(sp.expand(-residual/2))
    Gi=inverse_coeff(G,N);Bi=inverse_coeff(B,N)
    U=[sp.Integer(0)]+[sp.expand(sum(i*G[i]*D[a]*Gi[n-i-a] for i in range(1,min(d,n)+1) for a in range(n-i+1))) for n in range(1,N+1)]
    V=[sp.Integer(0)]+[sp.expand(sum(i*B[i]*D[a]*Bi[n-i-a] for i in range(1,min(e,n)+1) for a in range(n-i+1))) for n in range(1,N+1)]
    Ir=Ring(vars,weights,[(D[n],n) for n in range(m+1,N+1)])
    Jr=Ring(vars,weights,[(U[n],n) for n in range(m+1,m+d+1)]+[(V[n],n) for n in range(m+1,m+e+1)])
    w=sp.Symbol('w')
    rho=sp.resultant(sum(G[i]*w**(d-i) for i in range(d+1)),sum(B[i]*w**(e-i) for i in range(e+1)),w) if d and e else sp.Integer(1)
    E=[sp.expand(sum(G[i]*D[n-i] for i in range(min(n,d)+1))) for n in range(m+d+1)]
    return vars,weights,Ir,Jr,rho,E

def cols(polys,vars,weights,n):
    vv=[vec(f,vars,weights,n) for f in polys]
    return sp.Matrix.hstack(*vv) if vv else sp.zeros(len(mons(weights,n)),0)

def multiplication(vars,weights,n,exp):
    outn=n+sum(a*w for a,w in zip(exp,weights)); pos={a:i for i,a in enumerate(mons(weights,outn))}; M=sp.zeros(len(pos),len(mons(weights,n)))
    for j,a in enumerate(mons(weights,n)):M[pos[tuple(x+y for x,y in zip(a,exp))],j]=1
    return M

def audit(k,ell,m,derivatives=False):
    d=k-m;e=ell-m; pv,pw,Ir,Jr,rho,E=make_parameter(d,e,m)
    c=sp.symbols('c1:'+str(k+1)); cw=tuple(range(1,k+1));N=k*ell
    H=inverse_coeff([sp.Integer(1)]+list(c),ell+k)
    R=Ring(c,cw,[(H[n],n) for n in range(ell+1,ell+k+1)])
    row=sp.Matrix([[sp.diff(H[ell+i],c[j-1]) for j in range(m,k+1)] for i in range(1,k+1)])
    s=k-m+1; minors=[sp.expand(row[list(rr),:].det()) for rr in itertools.combinations(range(k),s)]
    RJ=Ring(c,cw,list(R.gens)+[(f,sum(ell+i+1 for i in rr)-sum(range(m,k+1))) for f,rr in zip(minors,itertools.combinations(range(k),s))])
    ttop=m*(d+e)+d*e; bound=max(ttop,N)
    substitution=dict(zip(c,E[1:])); image={}; results=[]
    for t in range(bound+1):
        mapped=[mul(c,a).subs(substitution,simultaneous=True) for a in mons(cw,t)]
        image[t]=cols(mapped,pv,pw,t)
    for t in range(bound+1):
        imdim=(Ir.quot(t)*image[t]).rank()
        assert imdim==RJ.dim(t),('parameter-kernel',k,ell,m,t,imdim,RJ.dim(t))
        # Compare actual colon multiplication kernel with the actual ideal J/I, degree by degree.
        product=cols([rho*mul(pv,a) for a in mons(pw,t)],pv,pw,t+d*e)
        kernel=Ir.dim(t)-(Ir.quot(t+d*e)*product).rank()
        JmodI=(Ir.quot(t)*Jr.mat(t)).rank()
        assert kernel==JmodI,('colon',k,ell,m,t,kernel,JmodI)
        assert Ir.quot(t+d*e)*cols([rho*f for f,wgt in Jr.gens if wgt==t],pv,pw,t+d*e)==sp.zeros(Ir.dim(t+d*e),len([f for f,wgt in Jr.gens if wgt==t]))
        if t<=ttop:
            constraints=[]
            for j in range(ttop-t+1):
                target=t+j; span=sp.Matrix.hstack(Ir.mat(target),image[target])
                ns=span.T.nullspace(); annihilator=sp.Matrix.vstack(*(v.T for v in ns)) if ns else sp.zeros(0,span.rows)
                for a in mons(pw,j): constraints.append(annihilator*multiplication(pv,pw,t,a))
            C=sp.Matrix.vstack(*constraints) if constraints else sp.zeros(0,len(mons(pw,t)))
            cond=len(mons(pw,t))-C.rank()-Ir.mat(t).rank()
            rhoimage=cols([rho*mul(pv,a) for a in mons(pw,t-d*e)],pv,pw,t)
            rhodim=(Ir.quot(t)*rhoimage).rank()
            assert cond==rhodim,('conductor',k,ell,m,t,cond,rhodim)
        results.append({'weight':t,'T':Ir.dim(t),'RmodJacobian':RJ.dim(t),'parameter_image':imdim,'colon_kernel':kernel,'JmodI':JmodI})
    # Verify generators individually map to zero; no Hilbert formula is used.
    for f,w in RJ.gens:
        assert Ir.quot(w)*vec(f.subs(substitution,simultaneous=True),pv,pw,w)==sp.zeros(Ir.dim(w),1)
    derivative_results=[]
    if derivatives:
        for mm in range(1,k+1):
            js=k-mm+1; jac=sp.Matrix([[sp.diff(H[ell+i],c[j-1]) for j in range(mm,k+1)] for i in range(1,k+1)])
            rr=list(itertools.combinations(range(k),js)); jj=[sp.expand(jac[list(a),:].det()) for a in rr]
            RQ=Ring(c,cw,list(R.gens)+[(f,sum(ell+i+1 for i in a)-sum(range(mm,k+1))) for f,a in zip(jj,rr)])
            sw=cw[:mm]
            for t in range(N+1):
                src=mons(sw,t+mm); dst=mons(sw,t)
                embsrc=cols([mul(c[:mm],a) for a in src],c,cw,t+mm)
                kernels=(R.quot(t+mm)*embsrc).nullspace()
                K=sp.Matrix.hstack(*kernels) if kernels else sp.zeros(len(src),0)
                diff=sp.zeros(len(dst),len(src));pos={a:i for i,a in enumerate(dst)}
                for j,a in enumerate(src):
                    if a[-1]:
                        b=a[:-1]+(a[-1]-1,);diff[pos[b],j]=a[-1]
                W=diff*K
                embdst=cols([mul(c[:mm],a) for a in dst],c,cw,t)
                Qconstraints=RQ.quot(t)*embdst; qdim=len(dst)-Qconstraints.rank()
                assert Qconstraints*W==sp.zeros(Qconstraints.rows,W.cols),('derivative-inclusion',mm,t)
                assert W.rank()==qdim,('derivative-equality',mm,t,W.rank(),qdim)
                derivative_results.append({'m':mm,'target_weight':t,'source_weight':t+mm,'source_above_Rtop':t+mm>N,'Wdim':W.rank(),'Qdim':qdim})
    print('PASSED',k,ell,m,'Tdim',sum(r['T'] for r in results),'image_dim',sum(r['parameter_image'] for r in results),flush=True)
    return {'parameters':{'k':k,'ell':ell,'m':m,'d':d,'e':e},'weighted_components':results,'derivative_components':derivative_results}

report={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest() if SOURCE else None,'method':'Direct weighted Macaulay ideal multiples and exact rational kernels; formal coefficient equations; actual Jacobian minors; full multiplication conductor constraints. No target Hilbert formula, Schur/Pieri multiplication, transition matrix or old program is used.','cases':[]}
for args in [(2,2,1,True),(2,3,2,False),(2,2,2,False)]:
    report['cases'].append(audit(*args)); OUT.write_text(json.dumps(report,indent=2),encoding='utf8')
report['elapsed_seconds']=time.time()-start;OUT.write_text(json.dumps(report,indent=2),encoding='utf8')
print('ALL PASSED in',report['elapsed_seconds'],'seconds',flush=True)
