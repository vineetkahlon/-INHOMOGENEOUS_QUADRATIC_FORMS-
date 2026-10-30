"""Discover rational certificates with exact handling of strict failure rows."""
import argparse,itertools,json,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from model import DIM,NAMES,root_rows,witness_affine,witness_rows
from exact_linear import solve_exact
from verify_cover import verify,verify_leaf

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

def discover(L,seconds,radius,box=None,shape_box=None,resume=None,region=None):
    if resume:box=resume.get('d_interval');shape_box=resume.get('shape_box')
    if resume:region=resume.get('region')
    root=[(a,b,False) for a,b in root_rows(L,2,box,shape_box,region)]
    zs=sorted(itertools.product(range(-radius,radius+1),repeat=3),key=lambda z:sum(x*x for x in z))
    coeff=np.array([witness_affine(z)[0] for z in zs],dtype=float)
    costs=np.array([1+sum(x*x for x in z) for z in zs])
    nodes=[None];stack=[(0,[])];previous=0;maxdepth=0
    if resume:
        if resume.get('model_version')!=2 or resume['L']!=L or resume['complete'] or resume['target']!='0<f<=d':raise ValueError('Bad resume')
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
                point=r.x[:DIM];values=coeff@point;js=-np.floor(values).astype(int);rem=values+js
                eligible=(rem>1e-7)&(rem<=point[0]+1e-8)
                score=np.where(eligible,(np.minimum(rem,point[0]-rem)+1e-6)/costs,-np.inf)
                best=int(np.argmax(score))
                if not eligible[best]:
                    failure={'reason':'No witness in finite set at outer-model point','point':dict(zip(NAMES,point.tolist()))}
                    stack.append((idx,path));break
                z=list(zs[best]);j=int(js[best]);left=len(nodes);right=left+1
                nodes.extend([None,None]);nodes[idx]={'z':z,'j':j,'children':[left,right]}
                cuts=witness_rows(z,j)
                stack.append((right,path+[(*cuts[1],True)]));stack.append((left,path+[(*cuts[0],False)]))
        except ValueError as exc:
            failure={'reason':str(exc)};stack.append((idx,path));break
        if time.monotonic()-last>10:
            print(json.dumps({'L':L,'checked':checked,'pending':len(stack),'max_depth':maxdepth}),flush=True);last=time.monotonic()
    proof={'L':L,'model_version':2,'target':'0<f<=d','complete':not stack,'nodes':nodes,'discovery_seconds':previous+time.monotonic()-start,'max_depth':maxdepth,
           'pending':[[idx,[[list(map(str,a)),str(b),s] for a,b,s in path]] for idx,path in stack]}
    if box is not None:proof['d_interval']=list(map(str,box))
    if shape_box is not None:proof['shape_box']=[list(map(str,pair)) for pair in shape_box]
    if region is not None:proof['region']=region
    if failure:proof['failure']=failure
    if proof['complete']:proof['verification']=verify(proof)
    return proof

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--L',type=int,choices=[1,2,3,4],default=2);ap.add_argument('--seconds',type=float,default=120)
    ap.add_argument('--radius',type=int,default=3);ap.add_argument('--lo',type=F);ap.add_argument('--hi',type=F)
    ap.add_argument('--shape-box',type=json.loads);ap.add_argument('--resume',type=Path);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--region')
    args=ap.parse_args();box=(args.lo,args.hi) if args.lo is not None else None
    p=discover(args.L,args.seconds,args.radius,box,args.shape_box,json.loads(args.resume.read_text()) if args.resume else None,args.region)
    args.output.write_text(json.dumps(p,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in p.items() if k not in ('nodes','pending')},indent=2))
