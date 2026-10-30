"""Generic exact dyadic cover builder and independent checker for fixed families (see families.py)."""
from fractions import Fraction as F
import math

def _bounds(n,j,const,box):
    a0,a1,k0,k1,m0,m1=box
    klo,khi=(k0,k1) if n>0 else (k1,k0)
    return n*n*a0+n*klo+const+n*j, n*n*a1+n*khi+const+n*j      # min,max of f+mu

def leaf_ok(n,j,const,box,modes):
    lo,hi=_bounds(n,j,const,box)
    if lo-box[5]<=0: return False
    for c1,c0 in modes:
        m=box[4] if 1+c1>=0 else box[5]
        if hi-(1+c1)*m<c0: return True
    return False

def inside(box,nb):
    return any(box[0]>=a and box[1]<=b and box[2]>=k and box[3]<=l and box[4]>=mu for a,b,k,l,mu in nb)

def in_pinch_ball(box,ball):
    """box (alpha,kappa ranges) inside the closed l1-ball of radius r about (a*,k*) modulo Z^2"""
    a0,a1,k0,k1=box[:4]; pa,pk,r=ball
    for i in (-1,0,1):
        for j in (-1,0,1):
            if max(abs(x-pa-i)+abs(y-pk-j) for x in (a0,a1) for y in (k0,k1))<=r: return True
    return False

class Family:
    def __init__(s,K,nmax,mu0,mu1,modes,tau_modes=(),nbhd=(),pinch=()):
        s.K=K; s.nmax=nmax; s.mu0=F(mu0); s.mu1=F(mu1); s.modes=list(modes); s.tau_modes=list(tau_modes); s.nbhd=list(nbhd)
        s.pinch=list(pinch)   # list of (alpha*, kappa*, r): pinch balls proved by Lemma (ladder)
        s.cache={n:K(n) for n in range(-nmax,nmax+1) if n}
    def modes_for(s,box): return s.modes+(s.tau_modes if inside(box,s.nbhd) else [])
    def const(s,n,w):
        d=dict(s.cache[n]); return d[tuple(w)]

def find(fam,box):
    md=fam.modes_for(box)
    for n in sorted(fam.cache,key=abs):
        for w,cst in fam.cache[n]:
            lo,_=_bounds(n,0,cst,box); lo-=box[5]
            j0=math.floor(-lo/n) if n>0 else math.ceil(-lo/n)
            for j in (j0-1,j0,j0+1):
                if leaf_ok(n,j,cst,box,md): return [n,list(w),j]
    return None

def build(fam,maxdepth=40):
    cnt={"leaves":0,"nodes":0}
    def rec(box,depth):
        cnt["nodes"]+=1
        for k,b in enumerate(fam.pinch):
            if in_pinch_ball(box,b): cnt["pinch"]=cnt.get("pinch",0)+1; return {"P":k}
        w=find(fam,box)
        if w: cnt["leaves"]+=1; return {"x":w}
        if depth>=maxdepth: raise RuntimeError("stuck at %s"%[str(x) for x in box])
        a0,a1,k0,k1,m0,m1=box; am=(a0+a1)/2; km=(k0+k1)/2; mm=(m0+m1)/2
        quad=[(a0,am,k0,km,m0,m1),(am,a1,k0,km,m0,m1),(a0,am,km,k1,m0,m1),(am,a1,km,k1,m0,m1)]
        mus=[(a0,a1,k0,k1,m0,mm),(a0,a1,k0,k1,mm,m1)]
        if m1>m0:
            sq=sum(1 for b in quad if find(fam,b))/4; sm=sum(1 for b in mus if find(fam,b))/2
            if sm>sq or (sm==sq and (a1-a0)<F(1,512)):
                return {"m":[rec(b,depth+1) for b in mus]}
        return {"q":[rec(b,depth+1) for b in quad]}
    return rec((F(0),F(1),F(0),F(1),fam.mu0,fam.mu1),0),cnt

def check(fam,tree):
    """Independent replay: reconstruct every box from the path, recompute the witness constant from
    the family definition, check strict inequalities exactly, require complete children."""
    st={"nodes":0,"leaves":0}
    def rec(node,box):
        st["nodes"]+=1
        if "P" in node:
            k=int(node["P"])
            if not (0<=k<len(fam.pinch) and in_pinch_ball(box,fam.pinch[k])): raise ValueError("pinch leaf outside ball")
            st["pinch_leaves"]=st.get("pinch_leaves",0)+1; return
        if "x" in node:
            n,w,j=node["x"]; n=int(n); j=int(j)
            if n==0 or abs(n)>fam.nmax: raise ValueError("bad n")
            if not leaf_ok(n,j,fam.const(n,w),box,fam.modes_for(box)): raise ValueError("leaf fails at %s"%[str(x) for x in box])
            st["leaves"]+=1; return
        a0,a1,k0,k1,m0,m1=box; am=(a0+a1)/2; km=(k0+k1)/2; mm=(m0+m1)/2
        if "q" in node:
            ch=node["q"]
            if len(ch)!=4: raise ValueError("quad node needs 4 children")
            for c,b in zip(ch,[(a0,am,k0,km,m0,m1),(am,a1,k0,km,m0,m1),(a0,am,km,k1,m0,m1),(am,a1,km,k1,m0,m1)]): rec(c,b)
        elif "m" in node:
            ch=node["m"]
            if len(ch)!=2: raise ValueError("mu node needs 2 children")
            for c,b in zip(ch,[(a0,a1,k0,k1,m0,mm),(a0,a1,k0,k1,mm,m1)]): rec(c,b)
        else: raise ValueError("unknown node")
    rec(tree,(F(0),F(1),F(0),F(1),fam.mu0,fam.mu1)); return st
