#!/usr/bin/env python3
"""Construct 312-leaf exact rectangle certificates for the m=3,K=1 resonance."""
from verify_diag_one_sixth import F, cover_one
import json
from pathlib import Path

def make_path(box,depth):
    l,u,v,z=box
    k=l.numerator*(2**depth)//l.denominator
    j=v.numerator*(2**depth)//v.denominator
    assert l==F(k,2**depth) and u==F(k+1,2**depth)
    assert v==F(j,2**depth) and z==F(j+1,2**depth)
    return ''.join(str(((j>>i)&1)*2+((k>>i)&1)) for i in range(depth-1,-1,-1))

OUT={"method":"Dyadic rational square covering, strict 0 < Q < 3/5", "cases":[]}
for b in (F(0),F(1,3),F(2,3)):
    c=F(1,2);e=f=F(0);cy=cz=F(1,2)
    count,depth,leaves=cover_one(b,c,e,cy,f,cz)
    rows=[]
    for box,dep,witness in leaves:
        rows.append({"path":make_path(box,dep),"witness":[int(z) for z in witness]})
    rows.sort(key=lambda v:v['path'])
    assert len(rows)==count
    OUT['cases'].append({"b":str(b),"c":str(c),"e":str(e),"cy":str(cy),"f":str(f),"cz":str(cz),"leaves":rows})
path=Path(__file__).with_name('case31_resonance_certificate.json')
path.write_text(json.dumps(OUT,separators=(',',':'),sort_keys=True))
print('wrote',path,'leaves',sum(len(c['leaves']) for c in OUT['cases']))
