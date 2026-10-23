"""Independent derivation of the near-top 3x3 Catalan pairing determinant."""
import sympy as s,json,hashlib
from pathlib import Path
m,t=s.symbols('m t')
records=[]
for sigma in (0,1):
    R=[s.S.One,(2*m+3*sigma+3*m*t)/(m+2+sigma)]
    R.append(s.cancel(((2*(m+1)+3*sigma)*R[1]+3*(m+1)*R[0])/(m+3+sigma)))
    R.append(s.cancel(((2*(m+2)+3*sigma)*R[2]+3*(m+2)*R[1])/(m+4+sigma)))
    c1=s.Rational(2)*(2*m+2*sigma+3)/(m+sigma+3)
    c2=c1*2*(2*m+2*sigma+5)/(m+sigma+4)
    b0=R[0]+3*R[1]+3*R[2]+R[3]
    b1=R[1]+2*R[2]+R[3]
    b2=R[2]+R[3]
    b3=R[3]
    K=s.Matrix([[c2,c2-c1,c2-2*c1+1],[b0,b1,b2],[b1,b2,b3]])
    K=K.applyfunc(s.cancel)
    terms=[K[0,0]*K[1,1]*K[2,2],K[0,1]*K[1,2]*K[2,0],K[0,2]*K[1,0]*K[2,1],-K[0,2]*K[1,1]*K[2,0],-K[0,1]*K[1,0]*K[2,2],-K[0,0]*K[1,2]*K[2,1]]
    det=s.factor(s.cancel(sum(s.cancel(term) for term in terms)))
    num,den=s.fraction(det)
    content=s.factor(s.gcd_list(s.Poly(num,t).all_coeffs()))
    primitive=s.factor(num/content)
    disc=s.factor(s.discriminant(primitive,t))
    expected_q = (9*m*(3*m+7)*t**2 + 3*(6*m**2+13*m-6)*t - (9*m**2+38*m+48)) if sigma==0 else (9*m*(2*m+7)*t**2 + 6*(2*m**2+9*m+6)*t - (6*m**2+37*m+60))
    expected_det = (2*m*(m+1)*expected_q/((m+2)**2*(m+3)**2*(m+4))) if sigma==0 else (3*m*(m+1)*expected_q/((m+3)**2*(m+4)**2*(m+5)))
    expected_disc = 81*(m+2)**2*(16*m**2+32*m+1) if sigma==0 else 144*(m+3)**2*(4*m**2+14*m+1)
    assert s.cancel(det-expected_det)==0
    assert s.expand(primitive-expected_q)==0
    assert s.expand(disc-expected_disc)==0
    assert s.degree(primitive,t)==2
    if sigma==0:
        assert s.expand(16*m**2+32*m+1-(4*m+3)**2-8*(m-1))==0
        assert s.expand((4*m+4)**2-(16*m**2+32*m+1)-15)==0
    else:
        assert s.expand(4*m**2+14*m+1-(2*m+3)**2-2*(m-4))==0
        assert s.expand((2*m+4)**2-(4*m**2+14*m+1)-(2*m+15))==0
        assert expected_q.subs({m:4,t:s.Rational(4,9)})==0
    record={'status':'PASS','sigma':sigma,'normalized_determinant':str(det),'content':str(content),'denominator':str(den),'primitive_quadratic':str(primitive),'discriminant':str(disc)}
    records.append(record);print(json.dumps(record),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'independent':True,'records':records,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),encoding='utf8')
