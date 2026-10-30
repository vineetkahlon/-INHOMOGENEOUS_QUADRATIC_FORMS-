"""Exact proof of the two-piece gap bound (Lemma 'two-piece'):
   g_4(t,beta) <= 1/2 - t on [1/16,1/10],   g_4(t,beta) <= 4t on [1/10,1/6].
On each cell of the collision arrangement every consecutive gap is affine, hence g_4 - bound is convex there;
so it suffices to check all pairwise intersections of the arrangement lines in [0,1/2] x [t_-,t_+]."""
from fractions import Fraction as F
from itertools import combinations
from gap import lines, maxgap
PIECES=[(F(1,16),F(1,10),(F(1,2),F(-1))),(F(1,10),F(1,6),(F(0),F(4)))]   # bound = c0 + c1 t
def check_piece(lo,hi,bound,M=4):
    c0,c1=bound; Ls=lines(lo,hi,M); verts=set()
    for (a1,b1,k1),(a2,b2,k2) in combinations(Ls,2):
        det=a1*b2-a2*b1
        if det==0: continue
        beta=(k1*b2-k2*b1)/det; t=(a1*k2-a2*k1)/det
        if 0<=beta<=F(1,2) and lo<=t<=hi: verts.add((beta,t))
    worst=max(maxgap(b,t,M)-(c0+c1*t) for b,t in verts)
    if worst>0: raise ValueError("two-piece bound fails on [%s,%s]"%(lo,hi))
    return len(verts),worst
def main():
    return [dict(interval=[str(lo),str(hi)],vertices=n,max_excess=str(w)) for lo,hi,b in PIECES for n,w in [check_piece(lo,hi,b)]]
if __name__=="__main__":
    print(main())
    from gap import exact_max
    print("sup_beta g4(31/250) =",exact_max(F(31,250),F(31,250),4)[0])
