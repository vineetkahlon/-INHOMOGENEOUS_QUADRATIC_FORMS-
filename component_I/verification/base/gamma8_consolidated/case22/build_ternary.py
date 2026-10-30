"""Numerical discovery; complete trees are accepted only by verify.py.

Incomplete searches and points without a witness are saved explicitly.
"""
import argparse,itertools,json,time,sys
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from model import DIM,NAMES,root_rows,witness_affine,witness_rows
from verify_ternary import verify,verify_leaf
from exact_linear import solve_exact

def farkas(rows):
    mat=np.array([[*a,b] for a,b in rows],dtype=float).T
    r=linprog(np.ones(len(rows)),A_eq=mat,b_eq=[0]*DIM+[-1],bounds=(0,None),method='highs')
    if not r.success:raise ValueError('No numerical Farkas multiplier')
    support=np.flatnonzero(r.x>1e-9).tolist()
    weights=solve_exact([[*rows[i][0],rows[i][1]] for i in support],[F(0)]*DIM+[F(-1)])
    cert=[[i,str(w)] for i,w in zip(support,weights) if w]
    verify_leaf(rows,cert)
    return cert

def discover(L,seconds,radius,resume=None,box=None,center=False,penalty=0):
    version=resume.get('model_version',1) if resume else 3
    if resume:box=resume.get('d_interval')
    root=root_rows(L,version,box)
    zs=sorted(itertools.product(range(-radius,radius+1),repeat=3),key=lambda z:sum(x*x for x in z))
    coeff=np.array([witness_affine(z)[0] for z in zs],dtype=float)
    costs=np.array([(1+sum(x*x for x in z))**penalty for z in zs])
    nodes=[None];stack=[(0,[])];maxdepth=0
    previous_seconds=0
    if resume:
        if resume['L']!=L or resume['complete']:raise ValueError('Invalid resume state')
        nodes=resume['nodes'];maxdepth=resume['max_depth'];previous_seconds=resume['discovery_seconds']
        stack=[(idx,[(list(map(F,a)),F(b)) for a,b in path]) for idx,path in resume['pending']]
    start=time.monotonic();last=start;checked=0;fail=None
    while stack:
        if time.monotonic()-start>seconds:break
        idx,path=stack.pop();rows=root+path;checked+=1;maxdepth=max(maxdepth,len(path))
        A=np.array([a for a,b in rows],dtype=float);b=np.array([b for a,b in rows],dtype=float)
        if center:
            norms=np.max(np.abs(A),axis=1)
            r=linprog([0]*DIM+[-1],A_ub=np.column_stack((A,norms)),b_ub=b,
                      bounds=[(None,None)]*DIM+[(0,1)],method='highs')
        else:r=linprog(np.zeros(DIM),A_ub=A,b_ub=b,bounds=[(None,None)]*DIM,method='highs')
        try:
            if r.status==2:nodes[idx]={'certificate':farkas(rows)}
            elif r.success:
                point=r.x[:DIM]
                vals=coeff@point;j=-np.floor(vals).astype(int);rem=vals+j
                margins=np.minimum(rem,point[0]-rem);best=int(np.argmax(margins/costs))
                if margins[best]<1e-7:
                    fail={'reason':'no strict witness at outer-model point','point':dict(zip(NAMES,point.tolist())),'margin':float(margins[best])}
                    stack.append((idx,path));break
                z=list(zs[best]);jj=int(j[best]);left=len(nodes);right=left+1
                nodes.extend([None,None]);nodes[idx]={'z':z,'j':jj,'children':[left,right]}
                cuts=witness_rows(z,jj)
                stack.append((right,path+[cuts[1]]));stack.append((left,path+[cuts[0]]))
            else:nodes[idx]={'certificate':farkas(rows)}
        except ValueError as exc:
            fail={'reason':str(exc)};stack.append((idx,path));break
        if time.monotonic()-last>10:
            print(json.dumps({'L':L,'checked':checked,'pending':len(stack),'maxdepth':maxdepth}),flush=True);last=time.monotonic()
    proof={'L':L,'model_version':version,'radius':radius,'complete':not stack,'nodes':nodes,'discovery_seconds':previous_seconds+time.monotonic()-start,'max_depth':maxdepth,
           'pending':[[idx,[[list(map(str,a)),str(b)] for a,b in path]] for idx,path in stack]}
    if box is not None:proof['d_interval']=list(map(str,box))
    if fail:proof['failure']=fail
    if proof['complete']:proof['verification']=verify(proof)
    return proof

def upgrade(proof):
    """Retain old certificates while appending valid geometric root cuts."""
    old=proof.get('model_version',1)
    if old==3:return proof
    n=len(root_rows(proof['L'],old));new=len(root_rows(proof['L'],3))
    if root_rows(proof['L'],old)!=root_rows(proof['L'],3)[:n]:raise ValueError('Root prefix changed')
    for node in proof['nodes']:
        if node and 'certificate' in node:
            node['certificate']=[[i if i<n else i+new-n,w] for i,w in node['certificate']]
    proof['model_version']=3
    return proof

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--L',type=int,default=4);ap.add_argument('--seconds',type=float,default=120)
    ap.add_argument('--radius',type=int,default=2);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--resume',type=Path)
    ap.add_argument('--upgrade',action='store_true');ap.add_argument('--lo',type=F);ap.add_argument('--hi',type=F)
    ap.add_argument('--center',action='store_true');ap.add_argument('--penalty',type=float,default=0)
    args=ap.parse_args();resume=json.loads(args.resume.read_text()) if args.resume else None
    if args.upgrade:resume=upgrade(resume)
    box=(args.lo,args.hi) if args.lo is not None else None
    p=discover(args.L,args.seconds,args.radius,resume,box,args.center,args.penalty)
    args.output.write_text(json.dumps(p,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in p.items() if k not in ('nodes','pending')},indent=2))
