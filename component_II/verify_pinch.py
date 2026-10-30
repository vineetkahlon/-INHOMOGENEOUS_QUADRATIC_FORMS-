"""Exact verification of the pinch-isolation data (Lemma 'ladder' and Theorem 'failure sets', Section 6).
For a pinch point p* of a resonant family f = n j + n^2 alpha + n kappa + n e (z+c) - q (z+c)^2, a rung (v,n) means
v is a value of f at p* for both n and -n.  Ascending/descending ladders start at n=1; consecutive ratios rho
must satisfy the handoff inequalities for all admissible (mu, Delta), Delta >= d0 + 3 mu, mu in [mu0, mu1].
Infinite tails are explicit polynomial families; identities are checked as polynomial identities via degree bounds."""
from fractions import Fraction as F
from partition import cbrt_down
h=F(1,2)

def val(q,e,c,al,ka,n,z,j): return n*j+n*n*al+n*ka+n*e*(z+c)-q*(z+c)**2

def check_rung(fam,p,v,n,wit_pos,wit_neg):
    q,e,c=fam; al,ka=p
    (z1,j1),(z2,j2)=wit_pos,wit_neg
    if val(q,e,c,al,ka,n,z1,j1)!=v or val(q,e,c,al,ka,-n,z2,j2)!=v: raise ValueError("rung fails %s"%((v,n),))

def poly_identity(fn,deg,lo=0):
    """fn(w) is a polynomial of degree <= deg in w; it vanishes identically iff it vanishes at deg+1 points"""
    return all(fn(F(w))==0 for w in range(lo,lo+deg+1))

PARAMS={  # family: d0, mu0, mu1 (Table 7)
 "L2_t1/4": (max(F(3803,10000),cbrt_down(F(5,72)))),
 "L3_t1/6": (max(F(3294,10000),cbrt_down(F(5,108)))),
}
def mu_range(L,t0,d0):
    dmax=min((L+1)*t0,1-F(L*L-1,4)*t0); return F(3,16)*d0,(dmax-d0)/3

def handoff(prefix,tail_ratio,d0,mu0,mu1,kind):
    """prefix: list of (v,n); after the prefix the ladder continues with value v_tail and consecutive ratio <= tail_ratio.
    Check (A) rho^2 (mu - v_k) < mu + Delta - v_{k+1}  or (D) rho^2 (v_k - mu - Delta) < v_{k+1} - mu,
    with Delta = d0+3mu (worst case; both sides monotone in Delta), at mu0 and mu1 (affine in mu)."""
    pairs=[(prefix[k][0],prefix[k+1][0],F(prefix[k+1][1],prefix[k][1])) for k in range(len(prefix)-1)]
    pairs.append((prefix[-1][0],prefix[-1][0],tail_ratio))
    for v,v2,rho in pairs:
        for mu in (mu0,mu1):
            D=d0+3*mu
            ok = rho**2*(mu-v) < mu+D-v2 if kind=="A" else rho**2*(v-mu-D) < v2-mu
            if not ok: raise ValueError("handoff %s fails: %s"%(kind,(v,v2,rho,mu)))
    return len(pairs)

def main():
    out={}
    # ---------- P1, P2: family (q,e,c)=(1/4,0,0), L=2, t0=1/4 ----------
    fam=(F(1,4),F(0),F(0)); d0=PARAMS["L2_t1/4"]; mu0,mu1=mu_range(2,F(1,4),d0)
    # descending ladder: v=3/4 at n=(z^2+3)/4, z odd, j=1 (P1); j = 1-(n+1)/2 shift for P2 (n odd)
    for z in (1,3,5,7,9):
        n=(z*z+3)//4; assert (z*z+3)%4==0 and n%2==1
        check_rung(fam,(0,0),F(3,4),n,(z,1),(z,-1))
        check_rung(fam,(h,h),F(3,4),n,(z,1-(n+1)//2),(z,-1+(n-1)//2))
    assert poly_identity(lambda z: F((z*z+3),4)*1-z*z/4-F(3,4),2)            # n j - z^2/4 = 3/4 identically
    # ratio ((z+2)^2+3)/(z^2+3) is decreasing in z>=1, max at z=1 equals 3
    crossP=lambda z: ((z+2)**2+3)**2-((z+4)**2+3)*(z*z+3)
    assert poly_identity(lambda z: crossP(z)-(8*z*z+32*z-8),4) and F(12,4)==3      # >0 for z>=1: ratio decreasing, max 3
    nD=handoff([(F(3,4),1)],F(3),d0,mu0,mu1,"D")
    # ascending ladder P1: v=0 at every n (z=0,j=0): ratio <= 2;  P2: v=0 at n odd (j=-(n+1)/2) and n=2
    for n in range(1,12):
        check_rung(fam,(0,0),F(0),n,(0,0),(0,0))
        if n%2==1: check_rung(fam,(h,h),F(0),n,(0,-(n+1)//2),(0,(n-1)//2))
    check_rung(fam,(h,h),F(0),2,(2,-1),(2,0))
    nA=handoff([(F(0),1)],F(2),d0,mu0,mu1,"A")                                # ratios 2,3/2,5/3,... <= 2
    r12=min(F(3,4)-mu1, mu0+d0+3*mu0-0)
    assert r12>h
    out["P1,P2"]=dict(d0=str(d0),mu=[str(mu0),str(mu1)],radius_lower=str(r12),handoffs=nA+nD)
    # values of G+mu at P1,P2 in (0,1]: exactly {3/4,1}  (n(n+1)/2 integer, z^2 = 0,1 mod 4)
    assert {(z*z)%4 for z in range(4)}=={0,1}
    # ---------- P3: family (1/6,1/6,0), L=3, t0=1/6, p*=(7/12,1/2) ----------
    fam=(F(1,6),F(1,6),F(0)); p=(F(7,12),h); d0=PARAMS["L3_t1/6"]; mu0,mu1=mu_range(3,F(1,6),d0)
    R=lambda n,z,j: 12*n*j+7*n*n+6*n+2*n*z-2*z*z
    assert all(12*val(*fam,*p,n,z,j)==R(n,z,j) for n in range(-6,7) if n for z in range(-6,7) for j in range(-3,4))
    # symmetry V_{-n}=V_n:  R(-n,-z,-j-1) = R(n,z,j)
    # R(-n,-z,-j-1)-R(n,z,j) has degree <= 2 in each of n,z,j, so vanishing on a 3x3x3 grid proves the identity
    assert all(R(-n,-z,-j-1)==R(n,z,j) for n in range(3) for z in range(3) for j in range(3))
    def wit(n,target):
        for z in range(-12*abs(n),12*abs(n)+1):
            r=7*n*n+6*n+2*n*z-2*z*z
            if (target-r)%(12*n)==0: return (z,(target-r)//(12*n))
        raise ValueError
    def rung3(v,n):
        z,j=wit(n,int(12*v)); check_rung(fam,p,v,n,(z,j),(-z,-j-1))
    # descending prefix: 3/4@1, 5/6@2, 3/4@3, 5/6@6, 3/4@9,11,17,19,27 ; tail 3/4 @ n=6w^2+3 (z=3w, j=-(7w^2+w+4)/2), w>=2
    pre=[(F(3,4),1),(F(5,6),2),(F(3,4),3),(F(5,6),6),(F(3,4),9),(F(3,4),11),(F(3,4),17),(F(3,4),19),(F(3,4),27),(F(3,4),57)]
    for v,n in pre: rung3(v,n)
    assert poly_identity(lambda w: R(6*w*w+3,3*w,-(7*w*w+w+4)/2)-9,4)
    assert all((7*w*w+w+4)%2==0 for w in range(2))                            # j integral (period 2 in w)
    # tail ratio (6(w+1)^2+3)/(6w^2+3) decreasing (difference of cross products = 4w^2+8w >= 0), at w=3: 99/57
    cross=lambda w: (6*(w+1)**2+3)**2-(6*(w+2)**2+3)*(6*w*w+3)      # >=0 iff ratio decreasing at w
    assert poly_identity(lambda w: cross(w)-(72*w*w+144*w),4)            # = 72w(w+2) >= 0
    nD=handoff(pre,F(99,57),d0,mu0,mu1,"D")
    # ascending prefix: 1/12@1, 1/12@3 ; tail 1/12 @ n=2z^2+1, z=0,2 mod 3 (j=-1-(7z^2+z)/6): 9,19,51,73,...
    preA=[(F(1,12),1),(F(1,12),3),(F(1,12),9),(F(1,12),19)]
    for v,n in preA: rung3(v,n)
    assert poly_identity(lambda z: R(2*z*z+1,z,-1-(7*z*z+z)/6)-1,4)
    assert all((7*z*z+z)%6==0 for z in range(6) if z%3 in (0,2))               # j integral (period 6 in z)
    tail=[2*z*z+1 for z in range(3,400) if z%3 in (0,2)]
    tr=max(F(b,a) for a,b in zip(tail,tail[1:])); assert tr==F(51,19)       # max consecutive ratio of the tail
    # beyond z=400: consecutive tail elements have z'-z in {1,2}, so the ratio is <= (2(z+2)^2+1)/(2z^2+1),
    # which decreases in z (cross-difference 8z^2+16z-4... checked below) and is < 51/19 for z>=6
    assert poly_identity(lambda z: (2*(z+2)**2+1)*(2*(z+1)**2+1)-(2*(z+3)**2+1)*(2*z*z+1)-(16*z*z+48*z+8),4)  # >0
    assert F(2*402**2+1,2*400**2+1) < F(51,19)
    nA=handoff(preA,F(51,19),d0,mu0,mu1,"A")
    r3=min(F(3,4)-mu1, mu0+d0+3*mu0-F(1,12)); assert r3>h
    # values of G+mu at P3 in (0,1]: exactly {1/12,3/4,5/6}  (congruence argument; also in verify_scalars)
    out["P3"]=dict(d0=str(d0),mu=[str(mu0),str(mu1)],radius_lower=str(r3),handoffs=nA+nD)
    return out
if __name__=="__main__": print(main())
