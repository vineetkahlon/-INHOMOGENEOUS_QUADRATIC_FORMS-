from fractions import Fraction as F
from itertools import combinations
import math, sys

def residues(beta,q,M):
    return sorted(((beta*z - q*z*z) % 1) for z in range(-M,M+1))

def maxgap(beta,q,M):
    r=residues(beta,q,M)
    g=max(b-a for a,b in zip(r,r[1:]))
    return max(g, r[0]+1-r[-1])

def lines(qlo,qhi,M):
    L=[]  # each line: A*beta + B*q = C
    for i in range(-M,M+1):
        for j in range(i+1,M+1):
            A=F(i-j); B=F(-(i*i-j*j))
            vals=[A*b+B*q for b in (F(0),F(1,2)) for q in (qlo,qhi)]
            for k in range(math.floor(min(vals))-1, math.ceil(max(vals))+2):
                L.append((A,B,F(k)))
    L+= [(F(1),F(0),F(0)),(F(1),F(0),F(1,2)),(F(0),F(1),qlo),(F(0),F(1),qhi)]
    return list(set(L))

def exact_max(qlo,qhi,M):
    Ls=lines(qlo,qhi,M); best=F(0); arg=None
    for (a1,b1,c1),(a2,b2,c2) in combinations(Ls,2):
        det=a1*b2-a2*b1
        if det==0: continue
        beta=(c1*b2-c2*b1)/det; q=(a1*c2-a2*c1)/det
        if 0<=beta<=F(1,2) and qlo<=q<=qhi:
            g=maxgap(beta,q,M)
            if g>best: best,arg=g,(beta,q)
    return best,arg
