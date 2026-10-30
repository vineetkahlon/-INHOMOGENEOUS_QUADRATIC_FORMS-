"""Supplementary checks of the reduction and certificate rejection behavior."""
from fractions import Fraction as F
from pathlib import Path
import copy,json
from model import FAMILIES,phi,witness_affine
from verify import verify


def rejected(call):
    try:call()
    except ValueError:return
    raise ValueError('Malformed certificate accepted')


def run():
    checks=0
    alpha,c1=F(2,7),F(-3,11)
    for i,f in enumerate(FAMILIES):
        kappa=c1+sum(b*c for b,c in zip(f['b'],f['c']))
        for n in (-4,-1,1,2,3):
            for w in ((-1,0,1),(0,0,0),(1,-1,0)):
                X=tuple(F(x)+c for x,c in zip(w,f['c']))
                j=2
                direct=(F(j)+c1+alpha*n+sum(b*x for b,x in zip(f['b'],X)))*n-phi(f['kind'],X)
                row,offset=witness_affine(i,n,w,j)
                if direct!=row[0]*alpha+row[1]*kappa+offset:raise ValueError('Incorrect original-form identity')
                checks+=1
    p=json.loads((Path(__file__).resolve().parent/'certificates/family_00.json').read_text())
    q=copy.deepcopy(p);q['complete']=False;rejected(lambda:verify(q))
    q=copy.deepcopy(p);q['nodes'][0]['children']=q['nodes'][0]['children'][:1];rejected(lambda:verify(q))
    q=copy.deepcopy(p)
    leaf=next(n for n in q['nodes'] if 'certificate' in n)
    leaf['certificate'][0][1]='-1';rejected(lambda:verify(q))
    return {'lifting_checks':checks,'invalid_certificate_checks':3,'result':'passed'}


if __name__=='__main__':print(json.dumps(run(),indent=2))
