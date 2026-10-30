#!/usr/bin/env python3
"""Standalone exact rational verifier for (m,K)=(2,4) resonance covers."""
import sys,json
from fractions import Fraction as Q
from verify_case32_cover import rect, verify_tree, require
PARAMS={
 (1,Q(0),Q(1,2)):(Q(1,2),Q(19,20),Q(301,4800),Q(27,280)),
 (1,Q(1,2),Q(0)):(Q(1,2),Q(19,20),Q(301,4800),Q(27,280)),
 (3,Q(1,6),Q(0)):(Q(1,6),Q(33,50),Q(7,160),Q(9,140)),
 (3,Q(0),Q(1,2)):(Q(1,6),Q(33,50),Q(7,160),Q(9,140)),
}

def check_scalars():
 # For K=4, Delta >= 35d/16, t<=6d^3.
 # L>=4 would require d^2>=35/[96(L+1)] and d<=64/[35(L+3)].
 # The ratio (L+3)^2/(L+1) increases for L>=1, so checking L=4 suffices.
 require(Q(35,96*5)>(Q(64,35*7))**2,'L>=4 not excluded')
 # L1: d^2 >=35/192 and d<=16/35; t>=Delta/2, t<=6d^3
 require(Q(35,192)>Q(42,100)**2,'L1 d lower')
 require(Q(35,16)*Q(42,100)>Q(2,3),'L1 t<1/2 side')
 require(6*Q(16,35)**3<Q(35,16)*Q(42,100),'L1 t>1/2 side')
 # L2: d^2 >=35/288 and d<=64/175; t>=Delta/3>1/4
 require(Q(35,288) > Q(12,35)**2,'L2 t not bounded away from 1/4')
 require(Q(35,288) > Q(34,100)**2, 'L2 lower bound')
 require(24*Q(64,175)**3-Q(1,2)<Q(35,16)*Q(34,100), 'L2 Macbeath width')
 # L3: d^2>=35/384 and d<=32/105; t>=Delta/4, t<=6d^3
 require(Q(35,384)>Q(3,10)**2, 'L3 lower')
 require(Q(35,16)*Q(3,10)>Q(8,13), 'L3 Macbeath low side')
 require(54*Q(32,105)**3-1<Q(35,16)*Q(3,10),'L3 Macbeath high side')
 # Resonant t=1/2 in L1 gives d^5=aA>=(d/3)(d/4)=d^2/12.
 # Resonant t=1/6 in L3 gives 3d^5=aA>=d^2/12 -> d^3>=1/36.
 require(Q(1,12)>Q(43,100)**3,'L1 resonance d lower')
 require(Q(1,36)>Q(3,10)**3,'L3 resonance d lower')
 for L,t,D,mu_lo,mu_hi,res_dlo_cub,res_dhi in (
       (1,Q(1,2),Q(19,20),Q(301,4800),Q(27,280),Q(1,12),Q(16,35)),
       (3,Q(1,6),Q(33,50),Q(7,160),Q(9,140),Q(1,36),Q(32,105))):
    require((D*Q(16,35))**3<res_dlo_cub, f'D too high L={L}')
    require(mu_hi == Q(27,128)*res_dhi, 'mu upper mismatch')
    dlo=Q(43,100) if L==1 else Q(3,10)
    require(mu_lo==Q(7,48)*dlo,'mu lower mismatch')
    keys=[(k,e,c) for (k,e,c) in PARAMS if k==L]
    require(len(keys)==2,'must have two families')
    for key in keys:
        require(PARAMS[key]==(t,D,mu_lo,mu_hi),'model mismatch')

def check_reps():
 for L,e,c in PARAMS:
    require(L in (1,3), 'bad L')
    require(c in (0,Q(1,2)), 'bad c')
    for sign in (1,-1):
        require((L*(sign*e-c/L-Q(1,2))).denominator==1,'invalid +/- resonance')
 for L in (1,3):
    for c in (Q(0),Q(1,2)):
        e=(Q(1,2)+c/L)%Q(1,L)
        require((L,e,c) in PARAMS,'missing normal form')

def main(path):
 check_scalars();check_reps()
 with open(path) as f:data=json.load(f)
 require(set(data)=={'description','cases'} and len(data['cases'])==4,'bad inventory')
 seen=set();rows=[]
 for item in data['cases']:
    require(set(item)=={'L','t','D','mu_lo','mu_hi','e','c','leaves'},'extra/missing fields')
    L=item['L'];require(type(L) is int,'noninteger L')
    e,c=Q(item['e']),Q(item['c']);key=L,e,c
    require(key in PARAMS and key not in seen,'missing/duplicate/foreign case')
    seen.add(key)
    params=PARAMS[key]
    require(tuple(Q(item[k]) for k in ('t','D','mu_lo','mu_hi'))==params,'altered parameters')
    t,D,ml,mh=params
    rows.append(verify_tree(item['leaves'],L,e,c,t,D,ml,mh))
 require(seen==set(PARAMS),'some cases omitted')
 expected={(1,Q(0),Q(1,2)):28,(1,Q(1,2),Q(0)):64,(3,Q(1,6),Q(0)):70,(3,Q(0),Q(1,2)):52}
 for s in rows:require(s['leaves']==expected[(s['L'],Q(s['e']),Q(s['c']))],'unexpected leaf count')
 return {'status':'PASS','families':rows,'total_leaves':sum(row['leaves'] for row in rows)}

if __name__=='__main__':
 print(json.dumps(main(sys.argv[1] if len(sys.argv)>1 else 'case24_certificate.json'),indent=2))
