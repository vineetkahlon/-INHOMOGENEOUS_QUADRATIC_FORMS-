"""Farkas-tree prover/checker over rational polytopes.  Rows: (dict var->Fraction, rhs Fraction) meaning sum <= rhs.
Witness value f(theta)=g.theta+h must satisfy 0 < f < T(theta), T affine. Failure branches: f<=0 ; f>=T."""
from fractions import Fraction as F
import numpy as np, math, json
from scipy.optimize import linprog

class Model:
    def __init__(s,vars_): s.V=list(vars_); s.ix={v:i for i,v in enumerate(s.V)}; s.rows=[]
    def le(s,coef,rhs):   # sum coef*var <= rhs
        s.rows.append(({k:F(v) for k,v in coef.items() if v!=0},F(rhs)))
    def ge(s,coef,rhs): s.le({k:-F(v) for k,v in coef.items()},-F(rhs))
    def box(s,v,lo,hi): s.le({v:1},hi); s.ge({v:1},lo)

def dense(rows,V,ix):
    B=np.zeros((len(rows),len(V))); c=np.zeros(len(rows))
    for i,(co,r) in enumerate(rows):
        for k,v in co.items(): B[i,ix[k]]=float(v)
        c[i]=float(r)
    return B,c

def center(rows,V,ix):
    B,c=dense(rows,V,ix); nrm=np.linalg.norm(B,axis=1)+1e-12
    A=np.hstack([B,nrm[:,None]]); obj=np.zeros(len(V)+1); obj[-1]=-1
    res=linprog(obj,A_ub=A,b_ub=c,bounds=[(None,None)]*len(V)+[(None,1)],method="highs")
    if res.status!=0: return None,None
    return res.x[:-1],res.x[-1]

def exact_solve(Mcols,rhs):
    """solve sum_k w_k Mcols[k] = rhs exactly (Fractions), Mcols: list of vectors. return list or None"""
    m=len(rhs); n=len(Mcols)
    A=[[Mcols[k][i] for k in range(n)]+[rhs[i]] for i in range(m)]
    piv=[]; r=0
    for col in range(n):
        p=next((i for i in range(r,m) if A[i][col]!=0),None)
        if p is None: continue
        A[r],A[p]=A[p],A[r]; pv=A[r][col]; A[r]=[x/pv for x in A[r]]
        for i in range(m):
            if i!=r and A[i][col]!=0:
                f=A[i][col]; A[i]=[x-f*y for x,y in zip(A[i],A[r])]
        piv.append(col); r+=1
        if r==m: break
    for i in range(r,m):
        if A[i][n]!=0: return None
    w=[F(0)]*n
    for i,col in enumerate(piv): w[col]=A[i][n]
    return w

def farkas(rows,V,ix):
    B,c=dense(rows,V,ix); m=len(rows)
    # w>=0, B^T w = 0, c^T w = -1
    Aeq=np.vstack([B.T,c[None,:]]); beq=np.zeros(len(V)+1); beq[-1]=-1
    res=linprog(np.ones(m),A_eq=Aeq,b_eq=beq,bounds=[(0,None)]*m,method="highs-ds")
    if res.status!=0: return None
    S=[i for i in range(m) if res.x[i]>1e-11]
    cols=[[rows[i][0].get(v,F(0)) for v in V]+[rows[i][1]] for i in S]
    w=exact_solve(cols,[F(0)]*len(V)+[F(-1)])
    if w is None or any(x<0 for x in w):
        # fallback: rationalise and accept only if exact identity holds
        w=[F(res.x[i]).limit_denominator(10**9) for i in S]
    cert={S[k]:w[k] for k in range(len(S)) if w[k]!=0}
    return cert if check_leaf(rows,V,cert) else None

def check_leaf(rows,V,cert):
    if any(x<0 for x in cert.values()): return False
    tot={v:F(0) for v in V}; rhs=F(0)
    for i,wt in cert.items():
        co,r=rows[i]
        for k,v in co.items(): tot[k]+=wt*v
        rhs+=wt*r
    return all(x==0 for x in tot.values()) and rhs<0

def prove(root,V,witnesses,target,maxnodes=200000,log=None):
    """witnesses(theta_float_dict) -> list of (label,g:dict,h) candidates ranked; target: (dict,const)"""
    ix={v:i for i,v in enumerate(V)}; stats={"nodes":0,"leaves":0}
    tco,tc=target
    def rec(rows):
        stats["nodes"]+=1
        if stats["nodes"]>maxnodes: raise RuntimeError("node limit")
        cert=farkas(rows,V,ix)
        if cert is not None:
            stats["leaves"]+=1; return {"leaf":{str(k):str(v) for k,v in cert.items()}}
        x,eps=center(rows,V,ix)
        if x is None: raise RuntimeError("LP failure")
        th={v:x[ix[v]] for v in V}
        best=None
        for lab,g,h in witnesses(th):
            fv=sum(float(g.get(k,0))*th[k] for k in g)+float(h)
            Tv=sum(float(tco.get(k,0))*th[k] for k in tco)+float(tc)
            mg=min(fv,Tv-fv)
            if best is None or mg>best[0]: best=(mg,lab,g,h)
        if best is None or best[0]<=0:
            raise RuntimeError("no witness at point %s"%{k:round(v,6) for k,v in th.items()})
        _,lab,g,h=best
        r1=(dict(g),-h)                                   # f<=0
        r2=({k:tco.get(k,F(0))-g.get(k,F(0)) for k in set(g)|set(tco)},h-tc)   # T-f<=0
        return {"w":lab,"c":[rec(rows+[r1]),rec(rows+[r2])]}
    tree=rec(list(root)); return tree,stats

def check_tree(root,V,tree,witness_affine,target):
    """independent replay: witness_affine(label)->(g,h); exact."""
    tco,tc=target; n=[0,0]
    def rec(rows,node):
        n[0]+=1
        if "leaf" in node:
            n[1]+=1; cert={int(k):F(v) for k,v in node["leaf"].items()}
            assert all(0<=k<len(rows) for k in cert); 
            if not check_leaf(rows,V,cert): raise ValueError("bad leaf")
            return
        g,h=witness_affine(node["w"]); g={k:F(v) for k,v in g.items()}; h=F(h)
        r1=(g,-h); r2=({k:tco.get(k,F(0))-g.get(k,F(0)) for k in set(g)|set(tco)},h-tc)
        if len(node["c"])!=2: raise ValueError("missing child")
        rec(rows+[r1],node["c"][0]); rec(rows+[r2],node["c"][1])
    rec(list(root),tree); return n
