"""(3,1): exact interval tests (Table 1) and residue-gap classes (Table 2)."""
from fractions import Fraction as F
from math import gcd
ROWS=[(F(5946,10000),F(5956,10000),8,F(19,2)),(F(5956,10000),F(5964,10000),10,F(15)),(F(5964,10000),F(6051,10000),4,F(5,2)),
      (F(6051,10000),F(6166,10000),3,F(3,2)),(F(6166,10000),F(6209,10000),5,F(9,2)),(F(6209,10000),F(6240,10000),4,F(3))]
GAPS=[(F(3,20),10,F(15),F(1,4)),(F(5,32),4,F(5,2),F(1,4)),(F(9,50),5,F(9,2),F(4,25)),(F(3,16),4,F(3),F(1,2)),(F(1,6),3,F(3,2),F(2,3))]
def rgcd(x,y): return F(gcd(x.numerator*y.denominator,y.numerator*x.denominator),x.denominator*y.denominator)
def main():
    assert F(5946,10000)**4*8<1 and F(624,1000)+4*F(624,1000)**5>1          # 8^{-1/4} > 0.5946 ; d<0.624
    for (l,u,N,H),(l2,_,_,_) in zip(ROWS,ROWS[1:]+[(F(6240,10000),0,0,0)]):
        assert max(abs(H-2*N*N*l**5),abs(H-2*N*N*u**5))+F(1,2)<l             # width < d on [l,u]
        assert u==l2
    out=[]
    for a,N,H,claim in GAPS:
        g=rgcd(F(1,N),2*a); P=2*int(1/a)*N*int(1/g)                           # a period multiple for x
        worst=0
        for k in range(int(1/g)):
            beta=H/N+k*g
            r=sorted({(beta*x-a*x*x)%1 for x in range(0,4*P)})
            worst=max(worst,max([b-c for c,b in zip(r,r[1:])]+[r[0]+1-r[-1]]))
        assert worst==claim, (a,worst)
        solved = worst**5 < a/2                                                 # d^5=a/2
        out.append((str(a),str(worst),solved))
    assert [s for *_,s in out]==[True,True,True,True,False]
    # a=1/6: with n=2 the value is 2j+4alpha+2kappa+K(w); K(w) mod 2 has periods 12,4,4 in x,y,z
    h=F(1,2); res=[]
    for b,c in [(F(1,6),F(0)),(F(0),h)]:
        for e,cy in [(F(0),h),(h,F(0))]:
            for f,cz in [(F(0),h),(h,F(0))]:
                K=lambda x,y,z: 2*b*(x+c)-F(1,6)*(x+c)**2+2*e*(y+cy)-h*(y+cy)**2+2*f*(z+cz)-h*(z+cz)**2
                for x in range(12):
                    for y in range(4):
                        for z in range(4):
                            v=K(x,y,z)
                            assert (K(x+12,y,z)-v)%2==0 and (K(x,y+4,z)-v)%2==0 and (K(x,y,z+4)-v)%2==0
                r=sorted({K(x,y,z)%2 for x in range(12) for y in range(4) for z in range(4)})
                g=max([q-p for p,q in zip(r,r[1:])]+[r[0]+2-r[-1]])
                allhalf=(b,c,e,cy,f,cz)==(0,h,0,h,0,h)
                assert (g==F(2,3)) if allhalf else (g==F(1,3))
                res.append(((str(b),str(e),str(f)),str(g)))
    return out,res
if __name__=="__main__": print(main())
