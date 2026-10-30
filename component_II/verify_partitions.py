"""Exact check of the t-partitions (G-level, x2=+-1): each interval closed by a uniform gap bound or by Macbeath."""
from fractions import Fraction as F
import json
from gap import exact_max
from partition import cbrt_down
def sqrt_down(x,den=10**4):
    import math
    k=int(math.sqrt(float(x))*den)+2
    while F(k,den)**2>x: k-=1
    return F(k,den)
CASES={  # (L, ell, delta_low(t), t_lo, t_hi)
 "c22_L2":(2,F(3803,10000)), "c22_L3":(3,F(3294,10000))}
def setup(key):
    L,ell=CASES[key]
    if key.startswith("c22"):
        assert ell*ell*288*(L+1)<125                          # ell <= sqrt(125/(288(L+1)))
        dlow=lambda t: max(L*t,F(25,16)*ell,F(25,16)*cbrt_down(F(5,18)*t))   # delta>Lt, delta>=25d/16, t<=18d^3/5
        return dlow,F(25,16)*ell/(L+1),F(4,L+3)/L
    assert ell*ell*128<25                                      # ell <= 5 sqrt2/16
    dlow=lambda t: max(L*t,F(5,4)*ell,F(5,4)*cbrt_down(F(5,12)*t)) # delta>Lt, delta>=5d/4, t<=12d^3/5
    return dlow,max(F(15,32)*ell,F(5,4)*ell/(L+1)),min(F(12,5)*F(8,15)**3,F(4,L*L+4*L-1))
def main():
    P=json.load(open("certs/partitions.json")); res={}
    for key,rows in P.items():
        dlow,tlo,thi=setup(key); cur=tlo; resonances=[]
        for r in rows:
            a,b=F(r["a"]),F(r["b"])
            if a!=cur: raise ValueError("%s: gap/overlap at %s"%(key,a))
            dl=dlow(a)
            if F(r["dlow"])>dl: raise ValueError("dlow overstated")
            if r["how"]=="gap":
                g,_=exact_max(a,b,int(r["M"]))
                if not g<dl or g!=F(r["gap"]): raise ValueError("gap claim fails %s"%r)
            elif r["how"] in ("macbeath","macbeath_res"):
                H,N=F(r["H"]),int(r["N"]); w=max(abs(H-N*N*a),abs(H-N*N*b))+F(1,2)   # convex in t
                if not w<dl: raise ValueError("Macbeath width fails %s"%r)
                res_t=H/(N*N); inside=a<=res_t<=b
                if inside!=(r["how"]=="macbeath_res"): raise ValueError("resonance flag wrong")
                if inside: resonances.append((str(res_t),str(H),N))
            else: raise ValueError("unclosed interval %s"%r)
            cur=b
        if cur!=thi: raise ValueError("%s does not reach t_hi"%key)
        res[key]=dict(intervals=len(rows),range=[str(tlo),str(thi)],resonances=resonances)
    return res
if __name__=="__main__":
    for k,v in main().items(): print(k,v)
