#!/usr/bin/env python3
from build_case32_cover import discover
from fractions import Fraction as Q
import json,sys
cases=[
 {'L':5,'t':Q(0),'D':Q(19,20),'mu_lo':Q(1173,19200),'mu_hi':Q(3899,45000)},
 {'L':6,'t':Q(0),'D':Q(99,100),'mu_lo':Q(18113,336000),'mu_hi':Q(298411,4100625)}
]
out={'description':'Universal third-completed H covers for (m,K)=(2,3), L=5,6', 'cases':[]}
for c in cases:
 r=discover(c,Q(0),Q(0),maxdepth=13)
 for row in r:
  n,z,j=row['witness'];row['witness']=[n,j]
 out['cases'].append({'L':c['L'],'D':str(c['D']),'mu_lo':str(c['mu_lo']),'mu_hi':str(c['mu_hi']),'leaves':r})
with open(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/case23_L56_certificate.json','w') as f:json.dump(out,f,indent=1)
