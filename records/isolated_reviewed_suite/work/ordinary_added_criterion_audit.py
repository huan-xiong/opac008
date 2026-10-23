"""Exact audit of the added coefficient criterion and odd excess formulas.
Coefficients are generated from h_n = u h_(n-1) - (u^2-v)h_(n-2).
Ranks are independently compared with top intersection pairings computed
from two-row Pieri/Catalan numbers, for both prescribed polynomial families.
This is a finite consistency check, not a proof of all-parameter assertions.
"""
from fractions import Fraction
from math import comb
from pathlib import Path
import json
import hashlib
from datetime import datetime, timezone

def rank(a):
    a = [[Fraction(x) for x in row] for row in a]
    if not a:
        return 0
    r = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(r, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        v = a[r][col]
        a[r] = [x/v for x in a[r]]
        for j in range(r+1, len(a)):
            v = a[j][col]
            if v:
                a[j] = [x-v*y for x,y in zip(a[j],a[r])]
        r += 1
        if r == len(a):
            break
    return r

MAX_K = 32
C = [[0]*(MAX_K+2) for _ in range(MAX_K+3)]
C[0][0] = C[1][0] = 1
for n in range(2, MAX_K+3):
    for a in range(n//2+1):
        C[n][a] = C[n-1][a]-C[n-2][a]+(C[n-2][a-1] if a else 0)
def coeff(n,a):
    return C[n][a] if 0 <= a <= n//2 else 0

def catalan(n):
    assert n >= 0
    return comb(2*n,n)//(n+1)
def moment(k,b):
    assert 0 <= b <= k
    return sum((-1)**i*comb(b,i)*catalan(k-i) for i in range(b+1))

results = []
for k in range(2,MAX_K+1):
    for d in range(1,k+1):
        n = k+d
        t = n//2
        p = t-d
        assert p >= 0
        matrix = [[coeff(k+1,a-j) for j in range((d-1)//2+1)]
                  +[coeff(k+2,a-j) for j in range((d-2)//2+1)]
                  for a in range(1,d+1)]
        assert len(matrix) == d and all(len(row)==d for row in matrix)
        indices = [0]+list(range(d+1,t+1))
        assert len(indices) == p+1
        gram_h = [[moment(k,a+j) for j in range(p+1)] for a in indices]
        gram_s = [gram_h[0]]+[
            [sum((-1)**i*comb(a+j,i)*catalan(k-t+a-i)
                 for i in range(a+j+1)) for j in range(p+1)]
            for a in indices[1:]]
        nullity = d-rank(matrix)
        defect_h, defect_s = p+1-rank(gram_h), p+1-rank(gram_s)
        assert nullity == defect_h == defect_s, (k,d,nullity,defect_h,defect_s)
        if d % 2 == 0 or d == 1:
            assert nullity == 0, (k,d,nullity)
        if d == 3:
            assert nullity == int(k % 3 == 2), (k,d,nullity)
        # Hook tests for the explicit exceptional core (k+d,d).
        core = [k+d,d]
        hooks = [[core[row]-col + sum(other > col for other in core[row+1:])
                  for col in range(core[row])] for row in range(2)]
        assert all(k+1 not in row for row in hooks)
        assert [sum(h <= k for h in row) for row in hooks] == [k,d]
        transpose = [2]*d+[1]*k
        transhooks = [[transpose[row]-col + sum(other > col for other in transpose[row+1:])
                      for col in range(transpose[row])] for row in range(k+d)]
        assert all(sum(h <= k for h in row)==1 for row in transhooks)
        results.append({'k':k,'d':d,'N':n,'dimension':p+1,'defect':nullity})
    print(f'PASS k={k}',flush=True)

# Verify all three residue rows, including n=0,1,2 boundaries.
residue_cases = 0
for r in range(0,11):
    rows = [
        [1,-2*r,Fraction(3*r*(r-1),2),Fraction(-r*(r-3)*(r-1),2)],
        [1,-r,-r,Fraction(r*(r-1)*(r+2),2)],
        [0,r+1,Fraction(-3*r*(r+1),2),r*(r-1)*(r+1)]
    ]
    for c in range(3):
        n=3*r+c
        for a in range(4):
            actual = (-1)**r*coeff(n,a)
            assert actual == rows[c][a], (n,a,actual,rows[c][a])
            residue_cases += 1

# Symbolic checks of the printed determinant and nonzero-minor formulas.
import sympy as sp
r=sp.symbols('r')
A0=sp.Matrix([[-r,1,r+1],[-r,-r,-sp.Rational(3,2)*r*(r+1)],
              [r*(r-1)*(r+2)/2,-r,r*(r-1)*(r+1)]])
A1=sp.Matrix([[r+1,0,2*(r+1)],[-sp.Rational(3,2)*r*(r+1),r+1,-sp.Rational(3,2)*r*(r+1)],
              [r*(r-1)*(r+1),-sp.Rational(3,2)*r*(r+1),r*(r-2)*(r+1)/2]])
A2=sp.Matrix([[2*(r+1),-1,r+1],[-sp.Rational(3,2)*r*(r+1),2*(r+1),r+1],
              [r*(r-2)*(r+1)/2,-sp.Rational(3,2)*r*(r+1),-r*(r+1)*(r+3)/2]])
assert sp.factor(A0.det()) == r**2*(r+1)**2*(3*r+2)/4
assert sp.factor(A1.det()) == r*(r+1)**3*(3*r+4)/4
assert A2.det().expand() == 0
assert sp.factor(A2[:2,:2].det()-(r+1)*(5*r+8)/2)==0
out={'status':'PASS','range':{'k_min':2,'k_max':MAX_K,'d':'1..k'},
     'criterion_pairing_cases':len(results),'coefficient_table_entries':residue_cases,
     'symbolic_determinants':3,'symbolic_minors':1,'results':results,
     'python_sympy':sp.__version__,'ended_utc':datetime.now(timezone.utc).isoformat(),
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('ordinary_added_criterion_audit_results.json').write_text(
    json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps({key:value for key,value in out.items() if key!='results'},ensure_ascii=False),flush=True)
