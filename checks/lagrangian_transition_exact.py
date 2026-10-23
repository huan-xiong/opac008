"""Exact check of the manuscript's two polynomial transition matrices.
This constructs the rows from the recurrence (not from the claimed determinant).
It checks their determinants, signs, and generator identities before quotienting.
"""
from pathlib import Path
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
m = s.symbols('M', positive=True)
results = []

def transition(d, parity, M):
    gs = (s.Integer(1),) + s.symbols('g1:'+str(d+1))
    us = s.symbols('u1:'+str(d+1))
    u = {j: us[j-1] for j in range(1,d+1)}
    def v(j):
        correct = (j % 2 == 0) if parity == 'even' else (j % 2 == 1)
        return -u[j] / (M+j) if correct else s.Integer(0)
    end = 2*d-2 if parity == 'even' else 2*d-1
    for t in range(d+1,end+1):
        u[t] = s.expand(-sum(gs[i]*u[t-i] for i in range(1,d+1)) +
                        sum(i*gs[i]*v(t-i) for i in range(1,d+1)))
    if parity == 'even':
        final = (u[d] + sum(gs[i]*u[d-i] for i in range(1,d)) -
                 sum(i*gs[i]*v(d-i) for i in range(1,d)))/(M+d)
        rows = [v(2*j) for j in range(1,d)] + [final]
    else:
        rows = [v(2*j-1) for j in range(1,d+1)]
    A = s.Matrix([[s.diff(row,x) for x in us] for row in rows])
    return gs, A

for d in range(1,6):
    for parity in ['even','odd']:
        gs,A = transition(d,parity,m)
        def g(i): return gs[i] if 0 <= i <= d else 0
        H = s.Matrix(d-1,d-1,lambda a,j:g(d+(a+1)-2*(j+1)))
        rho = H.det() if d>1 else s.Integer(1)
        if parity == 'even':
            target = (-1)**((d-1)*(d-2)//2)*rho/((m+d)*s.prod(m+2*j for j in range(1,d)))
        else:
            target = (-1)**(d*(d+1)//2)*rho/s.prod(m+2*j-1 for j in range(1,d+1))
        assert s.cancel(A.det()-target) == 0
        results.append({'type':'symbolic_determinant','d':d,'parity':parity,'M':'indeterminate','status':'PASS'})
        print('determinant',d,parity,'PASS',flush=True)

# Actual formal coefficients, computed through logarithmic recurrences.
for d in range(1,5):
    for M in range(1,7):
        parity = 'even' if M % 2 == 0 else 'odd'
        gs,A = transition(d,parity,s.Integer(M))
        z = s.symbols('z')
        G = sum(gs[i]*z**i for i in range(d+1))
        gg = s.Poly(s.expand(G*G.subs(z,-z)),z)
        ws = [gg.nth(2*i) for i in range(d+1)]
        nmax = M + 2*d
        D = [s.Integer(1)] + [s.Integer(0)]*nmax
        for a in range(1,nmax//2+1):
            D[2*a] = s.expand(-sum((s.Integer(a)-s.Rational(i,2))*ws[i]*D[2*(a-i)]
                                  for i in range(1,min(d,a)+1))/a)
        # GU=zG'D solves U without taking square roots in a quotient.
        U = [s.Integer(0)]*(nmax+1)
        for n in range(1,nmax+1):
            U[n] = s.expand(sum(i*gs[i]*D[n-i] for i in range(1,min(d,n)+1))-
                            sum(gs[i]*U[n-i] for i in range(1,min(d,n)+1)))
        seed = s.Matrix([U[M+t] for t in range(1,d+1)])
        if parity == 'even':
            f = [D[M+2*j] for j in range(1,d)] + [gs[d]*D[M]]
        else:
            f = [D[M+2*j-1] for j in range(1,d+1)]
        diff = A*seed-s.Matrix(f)
        assert all(s.expand(v)==0 for v in diff)
        results.append({'type':'polynomial_generator_identity','d':d,'M':M,'parity':parity,'status':'PASS'})
        print('generator identity',d,M,'PASS',flush=True)

out = {'status':'PASS','sympy':s.__version__,'symbolic_determinant_cases':10,
       'polynomial_generator_cases':24,'scope':'Exact polynomial identities over Q[g_1,...,g_d]; no quotient reduction, finite-range evidence only',
       'results':results}
(ROOT/'lagrangian_transition_exact_results.json').write_text(json.dumps(out,indent=2),encoding='utf8')
print('ALL PASS',flush=True)
