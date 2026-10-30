"""Discover rational certificates with exact handling of strict failure rows."""
import argparse,itertools,json,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from binary_resonance_model import DIM,NAMES,root_rows,witness_affine,witness_rows
from exact_linear import solve_exact
from verify_resonance import verify,verify_leaf

def farkas(rows):
    for zero in (False,True):
        columns=[[*a,b]+([int(strict)] if zero else []) for a,b,strict in rows]
        target=[F(0)]*DIM+([F(0),F(1)] if zero else [F(-1)])
        r=linprog(np.ones(len(rows)),A_eq=np.array(columns,dtype=float).T,b_eq=np.array(target,dtype=float),bounds=(0,None),method='highs')
        if not r.success:continue
        support=np.flatnonzero(r.x>1e-9).tolist()
        try:
            weights=solve_exact([columns[i] for i in support],target)
            cert=[[i,str(w)] for i,w in zip(support,weights) if w]
            verify_leaf(rows,cert);return cert
        except ValueError:continue
    raise ValueError('No exact infeasibility identity recovered')

def discover(family,seconds,radius,shape_box=None,resume=None):
    if resume:shape_box=resume.get('shape_box')
    root=[(a,b,False) for a,b in root_rows(shape_box)]
    pairs=[(n,z) for n in range(-4,5) if n for z in itertools.product(range(-radius,radius+1),repeat=2)]
    affine=[witness_affine(family,n,z) for n,z in pairs]
    coeff=np.array([a for a,b in affine],dtype=float);offsets=np.array([b for a,b in affine],dtype=float)
    periods=np.array([abs(n) for n,z in pairs]);signs=np.array([1 if n>0 else -1 for n,z in pairs])
    costs=np.array([1+sum(x*x for x in z)+(n*n-1)/2 for n,z in pairs])
    nodes=[None];stack=[(0,[])];previous=0;maxdepth=0
    if resume:
        if resume.get('model_version')!=2 or resume['family']!=family or resume['complete'] or resume['target']!='0<f<=d':raise ValueError('Bad resume')
        nodes=resume['nodes'];previous=resume['discovery_seconds'];maxdepth=resume['max_depth']
        stack=[(idx,[(list(map(F,a)),F(b),s) for a,b,s in path]) for idx,path in resume['pending']]
    start=time.monotonic();last=start;checked=0;failure=None
    while stack and time.monotonic()-start<seconds:
        idx,path=stack.pop();rows=root+path;checked+=1;maxdepth=max(maxdepth,len(path))
        A=np.array([a for a,b,s in rows],dtype=float);b=np.array([b for a,b,s in rows],dtype=float)
        strict=np.array([int(s) for a,b,s in rows]);norms=np.max(np.abs(A),axis=1)
        r=linprog([0]*DIM+[-1],A_ub=np.column_stack((A,norms)),b_ub=b,bounds=[(None,None)]*DIM+[(0,1)],method='highs')
        try:
            if r.success and r.x[-1]<1e-8 and any(strict):
                r=linprog([0]*DIM+[-1],A_ub=np.column_stack((A,strict)),b_ub=b,bounds=[(None,None)]*DIM+[(0,1)],method='highs')
                if r.success and r.x[-1]<1e-8:
                    nodes[idx]={'certificate':farkas(rows)};continue
            if not r.success:nodes[idx]={'certificate':farkas(rows)}
            else:
                point=r.x[:DIM];values=coeff@point+offsets;floors=np.floor(values/periods).astype(int);js=-signs*floors;rem=values-periods*floors
                eligible=(rem>1e-7)&(rem<=point[0]+1e-8)
                score=np.where(eligible,(np.minimum(rem,point[0]-rem)+1e-6)/costs,-np.inf)
                best=int(np.argmax(score))
                if not eligible[best]:
                    failure={'reason':'No witness in finite set at outer-model point','point':dict(zip(NAMES,point.tolist()))}
                    stack.append((idx,path));break
                n,z=pairs[best];z=list(z);j=int(js[best]);left=len(nodes);right=left+1
                nodes.extend([None,None]);nodes[idx]={'n':n,'z':z,'j':j,'children':[left,right]}
                cuts=witness_rows(family,n,z,j)
                stack.append((right,path+[(*cuts[1],True)]));stack.append((left,path+[(*cuts[0],False)]))
        except ValueError as exc:
            failure={'reason':str(exc)};stack.append((idx,path));break
        if time.monotonic()-last>10:
            print(json.dumps({'family':family,'checked':checked,'pending':len(stack),'max_depth':maxdepth}),flush=True);last=time.monotonic()
    proof={'family':family,'model_version':2,'target':'0<f<=d','complete':not stack,'nodes':nodes,'discovery_seconds':previous+time.monotonic()-start,'max_depth':maxdepth,
           'pending':[[idx,[[list(map(str,a)),str(b),s] for a,b,s in path]] for idx,path in stack]}
    if shape_box is not None:proof['shape_box']=[list(map(str,pair)) for pair in shape_box]
    if failure:proof['failure']=failure
    if proof['complete']:proof['verification']=verify(proof)
    return proof

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--family',type=int,choices=[0,1],required=True)
    ap.add_argument('--seconds',type=float,default=120);ap.add_argument('--radius',type=int,default=3)
    ap.add_argument('--shape-box',type=json.loads);ap.add_argument('--resume',type=Path)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    p=discover(args.family,args.seconds,args.radius,args.shape_box,json.loads(args.resume.read_text()) if args.resume else None)
    args.output.write_text(json.dumps(p,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in p.items() if k not in ('nodes','pending')},indent=2))
