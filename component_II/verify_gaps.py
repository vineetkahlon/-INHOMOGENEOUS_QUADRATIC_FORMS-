"""Exact uniform residue-gap bounds used in the paper (max over beta and over the q-interval, |z|<=M)."""
from fractions import Fraction as F
from gap import exact_max
ROWS=[ # name, qlo, qhi, M, claimed exact max
 ("U1 (2,3) L=1, F, A low", F(7,50),F(157,1000),4,F(57,100)),
 ("U2 (2,3) L=1, F, A high",F(22,125),F(1,5),4,F(3,5)),
 ("U3 (2,3) L=3, t low",    F(73,500),F(63,400),4,F(23,40)),
 ("U4 (2,3) L=3, t high",   F(22,125),F(27,125),4,F(3,5)),
 ("U5 (2,2) L=3 point t=1/5",F(1,5),F(1,5),4,F(3,5)),
]
def main():
    for nm,lo,hi,M,cl in ROWS:
        g,_=exact_max(lo,hi,M)
        if g!=cl: raise ValueError("%s: got %s"%(nm,g))
    import verify_twopiece
    return len(ROWS), verify_twopiece.main()
if __name__=="__main__":
    n,t=main(); print(n,"gap rows ok; two-piece bound:",t)
