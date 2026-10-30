"""Supplementary rational model/coordinate checks, not a replacement for the proof."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import copy,json
import ternary_model as T
import diagonal_model as D
from certificate import verify

def dot(a,b):return sum(x*y for x,y in zip(a,b))
def phi(a,A,t,h,k,l,w):
    x,y,z=w
    return a*(x+h*y+k*z)**2+A*(y+l*z)**2+t*z*z

def main():
    counts={'ternary_domain_points':0,'ternary_original_Q_identities':0,'diagonal_domain_points':0,'diagonal_original_Q_identities':0,'malformed_proofs_rejected':0}
    points=[]
    for d,ar,Ar,h,k,l in product([F(9,20),F(12,25),F(1,2),F(13,25),F(8,15)],
        [F(1,3),F(3,8),F(5,12)], [F(0),F(1,2),F(1)],
        [F(0),F(1,4),F(1,2)], [F(-1,2),F(-1,4),F(0),F(1,4),F(1,2)], [F(0),F(1,4),F(1,2)]):
        a=d*ar;delta=d+3*a/4
        if delta/2>1-delta:continue
        A=delta/2+Ar*(1-F(3,2)*delta);t=d**5/(2*a*A)
        if a**3>d**5 or t<A*(1-l*l):continue
        if any(phi(a,A,t,h,k,l,z)<a for z in product(range(-2,3),repeat=3) if 0<dot(z,z)<=5):continue
        points.append((d,a,A,t,h,k,l))
    # Deterministic spread through the whole candidate list, both shear signs.
    step=max(1,len(points)//100)
    roots={sign:T.root_rows(sign) for sign in (-1,1)}
    for d,a,A,t,h,k,l in points[::step]:
        beta=(F(-2,7),F(1,5),F(3,8));nu=F(2,5)
        theta=[d,a,A,t,a*h,a*k,a*h*h,a*h*k,a*k*k,A*l,A*l*l,*beta,nu]
        assert all(dot(row,theta)<=rhs for row,rhs in roots[1 if k>=0 else -1])
        counts['ternary_domain_points']+=1
        c=(F(1,3),F(-2,5),F(1,7));c1=F(1,7)
        gram=((a,a*h,a*k),(a*h,A+a*h*h,A*l+a*h*k),(a*k,A*l+a*h*k,t+A*l*l+a*k*k))
        ell=tuple(beta[i]+2*dot(gram[i],c) for i in range(3))
        a2=nu-c1-dot(ell,c)+phi(a,A,t,h,k,l,c)
        for z in [(0,0,0),(1,-1,0),(-1,1,-2),(2,0,1)]:
            j=1;shifted=tuple(F(x)+y for x,y in zip(z,c))
            original=j+c1+a2+dot(ell,shifted)-phi(a,A,t,h,k,l,shifted)
            row,offset=T.witness_affine(z,j)
            assert original==dot(row,theta)+offset
            counts['ternary_original_Q_identities']+=1
    diagroots=D.root_rows()
    for family,d,c,b,alpha,kappa in product(range(3),[F(16,25),F(33,50),F(17,25)],
         [F(0),F(1,6),F(1,3),F(1,2)],[F(-1,2),F(0),F(2,5)],[F(-1,3),F(1,4)],[F(-1,4),F(1,5)]):
        a=2*d**5;theta=[d,a,a*c,a*c*c,b,alpha,kappa]
        assert all(dot(row,theta)<=rhs for row,rhs in diagroots)
        counts['diagonal_domain_points']+=1
        (e,cy),(f,cz)=D.family_parameters(family);c1=kappa-b*c-e*cy-f*cz
        for n,z,j in [(1,(0,1,-1),0),(-1,(-1,0,2),1),(2,(1,-1,0),-1),(-3,(0,2,1),2)]:
            x,y,w=tuple(F(z[i])+v for i,v in enumerate((c,cy,cz)))
            original=(F(j)+c1+alpha*n+b*x+e*y+f*w)*n-a*x*x-(y*y+w*w)/2
            row,offset=D.witness_affine(family,n,z,j)
            assert original==dot(row,theta)+offset
            counts['diagonal_original_Q_identities']+=1
    p=json.loads((Path(__file__).resolve().parent/'certificates'/'low_plus.json').read_text())
    mutations=[]
    q=copy.deepcopy(p);q['complete']=False;mutations.append(q)
    q=copy.deepcopy(p);q['nodes'][0]['children'].pop();mutations.append(q)
    q=copy.deepcopy(p)
    leaf=next(n for n in q['nodes'] if 'certificate' in n);leaf['certificate'][0][1]='-1';mutations.append(q)
    q=copy.deepcopy(p);q['nodes'][0]['children'][0]=0;mutations.append(q)
    for q in mutations:
        try:verify(q)
        except ValueError:counts['malformed_proofs_rejected']+=1
        else:raise AssertionError('Malformed proof accepted')
    assert counts['ternary_domain_points']>0
    print(json.dumps(counts,indent=2))
if __name__=='__main__':main()
