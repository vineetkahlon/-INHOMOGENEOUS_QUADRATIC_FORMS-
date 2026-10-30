"""Check the new binary lifting identity and strict-certificate rejection rules."""
import copy,json
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from binary_resonance_model import NAMES,DIM,root_rows,witness_affine
from verify_resonance import verify,verify_leaf

def rejects(call):
    try:call()
    except ValueError:return
    raise ValueError('Invalid proof was accepted')

def run():
    # Exact determinant-compatible geometric point at t=1/6.
    d=F(21,50);a=F(9,40);A=3*d**5/a
    alpha0=F(2,7);kappa0=F(-3,11);checks=0
    for lam,q,C in product([F(1,4),F(1,3),F(1,2)],[F(-1,2),F(0),F(1,2)],[F(-1,2),F(0),F(1,2)]):
        data=dict(d=d,a=a,A=A,p=A*lam,v=A*lam**2,r=A*q,s=A*lam*q,
                  u=A*C,w=A*lam*C,mu=A*C*C,alpha=alpha0-A*q*q,kappa=kappa0+2*A*C*q)
        theta=[data[name] for name in NAMES]
        for row,rhs in root_rows():
            if sum(x*y for x,y in zip(row,theta))>rhs:raise ValueError('Physical point excluded')
        for family,n,y,z in product(range(2),[-3,-1,1,2,4],[-2,0,1],[-1,0,2]):
            c=F(family,2);e=F(1-family,6);j=-1
            direct=n*j+alpha0*n*n+kappa0*n+e*n*z-A*(y+lam*z+C-q*n)**2-(F(z)+c)**2/6-a/4
            row,offset=witness_affine(family,n,(y,z),j)
            if direct!=offset+sum(x*y for x,y in zip(row,theta)):raise ValueError('Incorrect lifting identity')
            checks+=1
    zero=[F(0)]*DIM
    rejects(lambda:verify_leaf([(zero,F(0),False)],[[0,'1']]))
    verify_leaf([(zero,F(0),True)],[[0,'1']])
    rejects(lambda:verify_leaf([(zero,F(0),True)],[[0,'-1']]))
    proof=json.loads((Path(__file__).resolve().parent/'certificates/resonance_0.json').read_text())
    broken=copy.deepcopy(proof);broken['nodes'][0]['children']=broken['nodes'][0]['children'][:1]
    rejects(lambda:verify(broken))
    broken=copy.deepcopy(proof);broken['complete']=False
    rejects(lambda:verify(broken))
    return {'exact_lifting_checks':checks,'physical_root_points':27,'invalid_certificate_checks':4,'result':'passed'}

if __name__=='__main__':print(json.dumps(run(),indent=2))
