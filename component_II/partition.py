"""Exact partition of a t-range: each interval closed by (i) Macbeath (H,N) with width < delta_low,
or (ii) uniform max-gap (window |z|<=M) < delta_low. Exact rationals throughout."""
from gap import *
import sys, json

def cbrt_down(x, den=10**5):
    k=int(float(x)**(1/3)*den)+2
    while F(k,den)**3 > x: k-=1
    return F(k,den)

def macbeath(a,b,dl,Nmax=6):
    """find (H,N) with max_{t in[a,b]} |H-N^2 t|+1/2 < dl; prefer no resonance inside."""
    best=None
    for N in range(1,Nmax+1):
        for H2 in range(1,2*N*N*1+2):
            H=F(H2,2); w=max(abs(H-N*N*a),abs(H-N*N*b))+F(1,2)
            if w<dl:
                res=H/(N*N); inside = a<=res<=b
                cand=(inside,N,H,w)
                if best is None or cand<best: best=cand
    return best

def partition(tlo,thi,dlow,minw=F(1,10**5),Ms=(4,6)):
    out=[]; stack=[(tlo,thi)]
    while stack:
        a,b=stack.pop()
        dl=dlow(a,b)
        m=macbeath(a,b,dl)
        if m and not m[0]:
            out.append(dict(a=a,b=b,how="macbeath",H=m[2],N=m[1],width=m[3],dlow=dl)); continue
        done=False
        for M in Ms:
            g,_=exact_max(a,b,M)
            if g<dl:
                out.append(dict(a=a,b=b,how="gap",M=M,gap=g,dlow=dl)); done=True; break
        if done: continue
        if m:   # Macbeath works except exactly at resonance H/N^2 inside
            out.append(dict(a=a,b=b,how="macbeath_res",H=m[2],N=m[1],width=m[3],dlow=dl,res=m[2]/(m[1]**2))); continue
        if b-a<minw:
            out.append(dict(a=a,b=b,how="FAIL",dlow=dl)); continue
        c=(a+b)/2; stack+= [(c,b),(a,c)]
    out.sort(key=lambda r:r["a"])
    return out

def show(out):
    for r in out:
        s={k:(str(v) if isinstance(v,F) else v) for k,v in r.items()}
        print(s, flush=True)
