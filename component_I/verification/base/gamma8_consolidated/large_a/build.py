"""Numerical discovery only; finished trees are checked as exact fractions."""
import argparse,itertools,json,time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
from model import DIM,FAMILIES,DELTAS,root_rows,witness_affine,witness_rows
from discovery_utils import farkas
from verify import verify


def candidates(index,maxn,radius):
    out=[]
    for n in range(-maxn,maxn+1):
        if not n:continue
        seen={}
        for w in sorted(itertools.product(range(-radius,radius+1),repeat=3),key=lambda w:sum(x*x for x in w)):
            _,offset=witness_affine(index,n,w)
            residue=offset%abs(n)
            if residue not in seen:seen[residue]=(n,w)
        out.extend(seen.values())
    return out


def discover(index,maxn=4,radius=3,seconds=30):
    root=root_rows();delta=DELTAS[FAMILIES[index]['kind']]
    pairs=candidates(index,maxn,radius)
    affine=[witness_affine(index,n,w) for n,w in pairs]
    coeff=np.array([a for a,b in affine],dtype=float)
    offset=np.array([b for a,b in affine],dtype=float)
    periods=np.array([abs(n) for n,w in pairs])
    signs=np.array([1 if n>0 else -1 for n,w in pairs])
    costs=np.array([1+(n*n-1)/2+sum(x*x for x in w)/10 for n,w in pairs])
    nodes=[None];stack=[(0,[])];start=time.monotonic();failure=None
    while stack and time.monotonic()-start<seconds:
        idx,path=stack.pop();rows=root+path
        A=np.array([a for a,b in rows],dtype=float);b=np.array([b for a,b in rows],dtype=float)
        result=linprog([0,0,-1],A_ub=np.column_stack((A,np.max(np.abs(A),axis=1))),
                       b_ub=b,bounds=[(None,None),(None,None),(0,1)],method='highs')
        try:
            if result.status==2:nodes[idx]={'certificate':farkas(rows,DIM)}
            elif result.success:
                vals=coeff@result.x[:DIM]+offset
                floors=np.floor(vals/periods).astype(int)
                js=-signs*floors;rem=vals-periods*floors
                margins=np.minimum(rem,float(delta)-rem)
                best=int(np.argmax(margins/costs))
                if margins[best]<1e-8:
                    failure={'point':result.x[:DIM].tolist(),'margin':float(margins[best])}
                    stack.append((idx,path));break
                n,w=pairs[best];j=int(js[best]);left=len(nodes);right=left+1
                nodes.extend([None,None]);nodes[idx]={'n':n,'w':list(w),'j':j,'children':[left,right]}
                cuts=witness_rows(index,n,w,j)
                stack.append((right,path+[cuts[1]]));stack.append((left,path+[cuts[0]]))
            else:raise ValueError(result.message)
        except ValueError as exc:
            failure={'reason':str(exc)};stack.append((idx,path));break
    proof={'family':index,'delta':str(delta),'complete':not stack,'nodes':nodes,
           'pending_count':len(stack),'discovery_seconds':time.monotonic()-start}
    if failure:proof['failure']=failure
    if proof['complete']:proof['verification']=verify(proof)
    return proof


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,default=0);ap.add_argument('--end',type=int,default=72)
    ap.add_argument('--maxn',type=int,default=4);ap.add_argument('--radius',type=int,default=3)
    args=ap.parse_args();root=Path(__file__).resolve().parent/'certificates';root.mkdir(exist_ok=True)
    for i in range(args.start,args.end):
        p=root/f'family_{i:02d}.json'
        if p.exists() and json.loads(p.read_text()).get('complete'):continue
        proof=discover(i,args.maxn,args.radius)
        p.write_text(json.dumps(proof,separators=(',',':'))+'\n')
        print(json.dumps({k:v for k,v in proof.items() if k!='nodes'}),flush=True)
