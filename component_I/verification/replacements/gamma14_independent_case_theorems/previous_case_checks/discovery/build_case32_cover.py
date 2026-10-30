#!/usr/bin/env python3
"""Discover rational quadtree covers for reduced (3,2) resonance families.
The separate verify_case32_cover.py rebuilds and verifies every inequality exactly.
This builder may use floating point to *find* data; floats are never trusted by checker.
"""
import math, json, sys, time
from fractions import Fraction as Q
from itertools import product
cases=[
  {'L':1,'t':Q(1,2),'D':Q(9,10),'mu_lo':Q(363,4000),'mu_hi':Q(5,48), 'classes':[(Q(0),Q(1,2)),(Q(1,2),Q(0))]},
  {'L':2,'t':Q(1,4),'D':Q(37,50),'mu_lo':Q(117,1600),'mu_hi':(Q(32,75)/3 + 8*Q(32,75)**4)/4, 'classes':[(Q(0),Q(0)),(Q(1,4),Q(1,2))]},
  {'L':3,'t':Q(1,6),'D':Q(16,25),'mu_lo':Q(102,1600),'mu_hi':(Q(16,45)/3 + 12*Q(16,45)**4)/4, 'classes':[(Q(1,6),Q(0)),(Q(0),Q(1,2))]}
]

# Floating discovery.
def discover(case,e,c,maxdepth=12):
    D=float(case['D']); M=float(case['mu_hi']-case['mu_lo']);lo_mu=float(case['mu_lo']); t=float(case['t']);e=float(e);c=float(c)
    nvals=list(range(1,6))+list(range(-1,-6,-1))
    W=[(n,z,n*e*(z+c)-t*(z+c)**2) for n in nvals for z in range(-8,9)]
    leaves=[]; branch=0; depthmax=0
    def go(path, al,au, kl,ku):
        nonlocal branch, depthmax
        amid=(al+au)/2; kmid=(kl+ku)/2; halfalpha=(au-al)/2; halfkappa=(ku-kl)/2
        for n,z,offset in W:
            halfspan=n*n*halfalpha+abs(n)*halfkappa+M/2
            if halfspan >= D/2-1e-11:continue
            center=n*n*amid+n*kmid+offset-lo_mu-M/2
            j=round((D/2-center)/n)
            mid=center+n*j
            if mid-halfspan>1e-10 and mid+halfspan<D-1e-10:
                leaves.append({'path':path,'witness':[n,z,j]})
                depthmax=max(depthmax,len(path))
                return
        if len(path)>=maxdepth:raise ValueError(f'bad depth={maxdepth} near {(al,au,kl,ku)} e,c={e,c}')
        branch+=1
        # square quadtree layout: bit 0 alpha, bit1 kappa
        for k in range(4):
            if k&1:a,b=amid,au
            else:a,b=al,amid
            if k&2:x,y=kmid,ku
            else:x,y=kl,kmid
            go(path+str(k),a,b,x,y)
    st=time.time()
    go('',0,1,0,1)
    print(' L',case['L'],'e,c',e,c,'leaves',len(leaves),'internal',branch,'maxdepth',depthmax,'duration',round(time.time()-st,2),flush=True)
    return leaves

if __name__=='__main__':
    data={'description':'Rational quadtree certificates for (m,K)=(3,2), using completed ternary form G and resonant t.',
          'cases':[]}
    for case in cases:
        for e,c in case['classes']:
            if case['L']==1 and e==Q(1,2):continue
            cert={'L':case['L'],'t':str(case['t']),'D':str(case['D']),'mu_hi':str(case['mu_hi']),'mu_lo':str(case['mu_lo']),'e':str(e),'c':str(c),'leaves':discover(case,e,c)}
            data['cases'].append(cert)
    with open(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/case32_certificate.json','w') as f:json.dump(data,f,indent=0,separators=(',',':'))
