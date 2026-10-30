"""Exact checks for the six-variable supplement (Proposition 1, Lemma 2, Theorem 3)."""
from fractions import Fraction as F
from itertools import product
from math import prod
def lhs(ms,n): return prod(F((m+3)**(n-i),m+1) for i,m in enumerate(ms,1))
def rhs(n,C,k): return F(2**(n-k),C)*4**(k*n-k*(k+1)//2)
def main():
    c=[]
    c.append(rhs(5,8,1)==512 and rhs(5,8,2)==4**7)                         # paper's thresholds are instances
    c.append(rhs(6,16,1)==2048 and rhs(6,16,2)==4**9)
    c.append(lhs([4],6)>rhs(6,16,1) and not lhs([3],6)>rhs(6,16,1))
    c.append(lhs([3,2],6)==405000 and lhs([2,3],6)==337500 and lhs([3,2],6)>4**9 and lhs([2,3],6)>4**9)
    c.append(not lhs([3,1],6)>4**9 and not lhs([2,2],6)>4**9)
    c.append(all((m+3)**3>=32*(m+1) for m in (1,2,3)))
    c.append(F(3,4)**3*F(128,27)==2 and F(9,4)*F(27,128)==F(243,512) and F(243,512)>F(4,9))
    c.append(F(8,3)*F(2,3)**5==F(256,729) and F(16,27)**2==F(256,729) and F(16,27)<F(2,3))
    c.append(4*F(2,3)**5==F(128,243) and F(128,243)<F(2,3) and F(31,48)**5<F(1,8))
    c.append(F(3,8)+F(1,18)==F(31,72) and 3*F(31,72)-1==F(7,24) and 3*F(3,8)-1==F(1,8) and 4*F(1,32)==F(1,8) and F(7,24)<F(31,32))
    h=F(1,2); res={}
    for hs in range(4):
        pairs=[(F(0),h)]*hs+[(h,F(0))]*(3-hs)
        r=sorted({sum(2*e*(y+x)-(y+x)**2/2 for (e,x),y in zip(pairs,ys))%2 for ys in product((0,1),repeat=3)})
        res[hs]=r
        if hs<=2: c.append(r==sorted(((-F(hs,8)+k*h)%2) for k in range(4)))
    c.append(res[3]==[F(5,8),F(13,8)])                                     # h=3: both parities, difference 1
    for (m,_) in [(4,0)]: pass
    # monotonicity of (m+3)^p/(m+1) for p>=3, m>=1: cross-difference positive on a range plus polynomial argument
    c.append(all(F((m+4)**p,m+2)>F((m+3)**p,m+1) for p in range(3,8) for m in range(1,60)))
    if not all(c): raise SystemExit("six-variable check failed: %s"%[i for i,x in enumerate(c) if not x])
    return len(c)
if __name__=="__main__": print(main(),"six-variable checks passed")
