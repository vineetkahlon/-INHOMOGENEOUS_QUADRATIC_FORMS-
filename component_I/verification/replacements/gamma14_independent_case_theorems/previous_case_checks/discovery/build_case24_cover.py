#!/usr/bin/env python3
"""Find exact rational rectangle-cover candidates for (2,4), L=1 and L=3 resonance."""
import json,sys
from fractions import Fraction as Q
from build_case32_cover import discover
cases=[
 {'L':1,'t':Q(1,2),'D':Q(19,20),'mu_lo':Q(301,4800),'mu_hi':Q(27,280),
  'classes':[(Q(0),Q(1,2)),(Q(1,2),Q(0))]},
 {'L':3,'t':Q(1,6),'D':Q(33,50),'mu_lo':Q(7,160),'mu_hi':Q(9,140),
  'classes':[(Q(1,6),Q(0)),(Q(0),Q(1,2))]},
]
if __name__=='__main__':
 out={'description':'(m,K)=(2,4), 4 exact resonance trees','cases':[]}
 for c in cases:
  for e,shift in c['classes']:
   rows=discover(c,e,shift,maxdepth=14)
   out['cases'].append({**{k:str(c[k]) for k in ('t','D','mu_lo','mu_hi')},'L':c['L'],'e':str(e),'c':str(shift),'leaves':rows})
 with open(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/case24_certificate.json','w') as f:json.dump(out,f,indent=1)
