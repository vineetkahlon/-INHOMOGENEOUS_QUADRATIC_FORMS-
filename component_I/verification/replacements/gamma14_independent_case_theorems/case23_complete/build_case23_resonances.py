#!/usr/bin/env python3
"""Certificate discovery for the three remaining resonant L=1,2,3 cases.

Discovery uses floating point *only* to locate candidate rectangles.
The separate exact checker, verify_case23_complete.py, validates all inequalities
and all dyadic partition trees using Python fractions.
"""
from fractions import Fraction as F
import json,sys,time

PARAMS = [
    (1, F(1,6), F(1,6), [(F(1,6),F(0)),(F(0),F(1,2))], F(3,5),F(1,25),F(53,1000)),
    (2, F(1,4), F(1,4), [(F(0),F(0)),(F(1,4),F(1,2))], F(7,10),F(31,192)*F(28,75),F(23,96)*F(32,75)),
    (3, F(1,6), F(1,6), [(F(1,6),F(0)),(F(0),F(1,2))], F(61,100),F(31,192)*F(163,500),F(7,100)),
]

def discover(t,D,mu_lo,mu_hi,e,c,maxdepth=15):
    # A common affine model for F (L=1) and G (L=2,3).
    D=float(D); mlo=float(mu_lo); mhi=float(mu_hi);e=float(e);c=float(c);t=float(t)
    spread=mhi-mlo
    witness_pool=[(n,z,n*e*(z+c)-t*(z+c)**2)
        for nn in range(1,6) for n in (nn,-nn) for z in range(-8,9)]
    out=[]
    def solve(path,a0,a1,k0,k1):
        ah=(a1-a0)/2;kh=(k1-k0)/2;am=(a0+a1)/2;km=(k0+k1)/2
        for n,z,offset in witness_pool:
            span=n*n*ah+abs(n)*kh+spread/2
            if span>=D/2-1e-10:continue
            center=n*n*am+n*km+offset-(mlo+mhi)/2
            j=round((D/2-center)/n)
            mid=center+n*j
            if mid-span>1e-10 and mid+span<D-1e-10:
                out.append({'path':path,'witness':[n,z,j]});return
        if len(path)>=maxdepth:raise RuntimeError('No witness at '+path)
        am=(a0+a1)/2;km=(k0+k1)/2
        for k in range(4):
            aa0,aa1=(am,a1) if k&1 else (a0,am)
            kk0,kk1=(km,k1) if k&2 else (k0,km)
            solve(path+str(k),aa0,aa1,kk0,kk1)
    solve('',0.,1.,0.,1.)
    return out

def main(output):
    result={'description':'Independent resonance covers for (m,K)=(2,3), branches L=1,2,3', 'cases':[]}
    for L,t,unused,classes,D,mlo,mhi in PARAMS:
        for e,c in classes:
            st=time.monotonic()
            leaves=discover(t,D,mlo,mhi,e,c)
            print(f'L={L} e={e} c={c}: {len(leaves)} leaves, {time.monotonic()-st:.2f}s',flush=True)
            result['cases'].append({'L':L,'t':str(t),'e':str(e),'c':str(c),'D':str(D),'mu_lo':str(mlo),'mu_hi':str(mhi),'leaves':leaves})
    with open(output,'w') as f:json.dump(result,f,separators=(',',':'),indent=1)

if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else 'case23_new_certificates.json')
