"""Find short exact certificates for the four remaining fixed families."""
import argparse,itertools,json,time,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from four_families_model import DIM,DELTA,root_rows,witness_affine,witness_rows
from verify_four import verify
from discovery_utils import farkas

def discover(family,seconds=100,maxn=4,radius=2):
    root=root_rows()
    pairs=[(n,z) for n in range(-maxn,maxn+1) if n for z in itertools.product(range(-radius,radius+1),repeat=2)]
    affine=[witness_affine(family,n,z) for n,z in pairs]
    coeff=np.array([a for a,b in affine],dtype=float);offset=np.array([b for a,b in affine],dtype=float)
    periods=np.array([abs(n) for n,z in pairs]);signs=np.array([1 if n>0 else -1 for n,z in pairs])
    costs=np.array([1+sum(x*x for x in z)+(n*n-1)/2 for n,z in pairs])
    nodes=[None];stack=[(0,[])];start=time.monotonic();failure=None
    while stack and time.monotonic()-start<seconds:
        idx,path=stack.pop();rows=root+path
        A=np.array([a for a,b in rows],dtype=float);b=np.array([b for a,b in rows],dtype=float)
        r=linprog([0]*DIM+[-1],A_ub=np.column_stack((A,np.max(np.abs(A),axis=1))),b_ub=b,
                  bounds=[(None,None)]*DIM+[(0,1)],method='highs')
        try:
            if not r.success:nodes[idx]={'certificate':farkas(rows,DIM)}
            else:
                vals=coeff@r.x[:DIM]+offset;floors=np.floor(vals/periods).astype(int)
                js=-signs*floors;rem=vals-periods*floors
                margins=np.minimum(rem,float(DELTA)-rem);best=int(np.argmax(margins/costs))
                if margins[best]<1e-7:
                    failure={'point':r.x[:DIM].tolist(),'margin':float(margins[best])};stack.append((idx,path));break
                n,z=pairs[best];j=int(js[best]);left=len(nodes);right=left+1
                nodes.extend([None,None]);nodes[idx]={'n':n,'z':list(z),'j':j,'children':[left,right]}
                cuts=witness_rows(family,n,z,j)
                stack.append((right,path+[cuts[1]]));stack.append((left,path+[cuts[0]]))
        except ValueError as exc:failure={'reason':str(exc)};stack.append((idx,path));break
    proof={'family':family,'delta':str(DELTA),'complete':not stack,'nodes':nodes,'discovery_seconds':time.monotonic()-start,
           'pending':[[idx,[[list(map(str,a)),str(b)] for a,b in path]] for idx,path in stack]}
    if failure:proof['failure']=failure
    if proof['complete']:proof['verification']=verify(proof)
    return proof

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--family',type=int,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();p=discover(args.family);args.output.write_text(json.dumps(p,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in p.items() if k not in ('nodes','pending')},indent=2))
