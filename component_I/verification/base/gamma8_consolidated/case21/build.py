"""Optional discovery: numerical LP followed by exact rational reconstruction.

The proof is checked by verify_all.py without NumPy, SciPy, or this file.
"""
import argparse,itertools,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from certificate import model,affine,cuts,verify
from exact_linear import solve_exact,verify_leaf

def farkas(rows,dim):
    columns=[[*a,b] for a,b in rows];target=[0]*dim+[-1]
    fit=linprog(np.ones(len(rows)),A_eq=np.array(columns,dtype=float).T,b_eq=target,bounds=(0,None),method='highs')
    if not fit.success:raise ValueError('No numerical Farkas identity')
    support=np.flatnonzero(fit.x>1e-9).tolist()
    weights=solve_exact([columns[i] for i in support],target)
    cert=[[i,str(w)] for i,w in zip(support,weights) if w]
    verify_leaf(rows,cert,dim)
    return cert

def discover(domain,seconds=120,radius=3):
    module,root=model(domain);dim=module.DIM
    zs=sorted(itertools.product(range(-radius,radius+1),repeat=3),key=lambda z:sum(x*x for x in z))
    candidates=[];periods=[];signs=[];costs=[]
    for n in ([1] if domain['kind']=='ternary' else [i for i in range(-4,5) if i]):
        for z in zs:
            w={'z':list(z),'j':0}
            if domain['kind']=='diagonal':w['n']=n
            candidates.append(w);periods.append(abs(n));signs.append(1 if n>0 else -1)
            costs.append(1+sum(x*x for x in z)+(n*n-1)/2)
    values=[affine(domain,w) for w in candidates]
    coeff=np.array([a for a,b in values],dtype=float);offsets=np.array([b for a,b in values],dtype=float)
    periods=np.array(periods);signs=np.array(signs);costs=np.array(costs)
    nodes=[None];stack=[(0,[])];start=time.monotonic();last=start
    while stack:
        if time.monotonic()-start>seconds:raise ValueError('Discovery time limit; no proof claimed')
        idx,path=stack.pop();rows=root+path
        A=np.array([a for a,b in rows],dtype=float);b=np.array([b for a,b in rows],dtype=float)
        norms=np.max(np.abs(A),axis=1)
        fit=linprog([0]*dim+[-1],A_ub=np.column_stack((A,norms)),b_ub=b,bounds=[(None,None)]*dim+[(0,1)],method='highs')
        if not fit.success:nodes[idx]={'certificate':farkas(rows,dim)}
        else:
            point=fit.x[:dim];raw=coeff@point+offsets;floors=np.floor(raw/periods).astype(int)
            rem=raw-periods*floors;eligible=(rem>1e-7)&(rem<point[0]-1e-7)
            score=np.where(eligible,np.minimum(rem,point[0]-rem)/costs,-np.inf);best=int(np.argmax(score))
            if not eligible[best]:raise ValueError('No finite-set witness at outer-model point; no proof claimed')
            w=dict(candidates[best]);w['j']=int(-signs[best]*floors[best]);left=len(nodes);right=left+1
            nodes.extend([None,None]);nodes[idx]={**w,'children':[left,right]}
            c=cuts(domain,w);stack.append((right,path+[c[1]]));stack.append((left,path+[c[0]]))
        if time.monotonic()-last>10:
            print(json.dumps({'nodes':len(nodes),'pending':len(stack)}),flush=True);last=time.monotonic()
    proof={'model_version':1,'target':'0<f<d','complete':True,'domain':domain,'nodes':nodes}
    verify(proof)
    return proof

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--kind',choices=['ternary','diagonal'],required=True)
    ap.add_argument('--sign',type=int,choices=[-1,1],default=1);ap.add_argument('--family',type=int,choices=[0,1,2],default=0)
    ap.add_argument('--seconds',type=float,default=120);ap.add_argument('--radius',type=int,default=3)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    domain={'kind':args.kind,**({'sign':args.sign} if args.kind=='ternary' else {'family':args.family})}
    proof=discover(domain,args.seconds,args.radius);args.output.write_text(json.dumps(proof,separators=(',',':'))+'\n')
    print(json.dumps(verify(proof),indent=2))
