#!/usr/bin/env python3
"""Exact, independent certificate checker for BKR case (m,K)=(3,1).

Only Python standard library is used. Verifies the finite ingredients; the
mathematical Macbeath reductions and their domain inclusions are supplied
in bkr_cases_41_31_note.pdf.

Run: python3 -S -B verify_bkr_case31.py case31_resonance_certificate.json
"""
import sys, json
from fractions import Fraction as F
from math import gcd, lcm
from itertools import product

DELTA=F(3,5)
M3=F(1,6)
D_LO=F(5946,10000)
D_HI=F(624,1000)
ROWS=((D_LO,F(5956,10000),8,F(19,2)),
      (F(5956,10000),F(5964,10000),10,F(15)),
      (F(5964,10000),F(6051,10000),4,F(5,2)),
      (F(6051,10000),F(6166,10000),3,F(3,2)),
      (F(6166,10000),F(6209,10000),5,F(9,2)),
      (F(6209,10000),D_HI,4,F(3)))
EF=((F(0),F(1,2)),(F(1,2),F(0)))
BC=tuple((F(i,6),F(0)) for i in (1,3,5))+tuple((F(i,3),F(1,2)) for i in range(3))

def need(c,msg):
    if not c:raise RuntimeError(msg)


def stepgcd(a,b):
    den=lcm(a.denominator,b.denominator)
    return F(gcd(int(a*den),int(b*den)),den)


def gap(vs,modulus):
    need(len(vs)>0,'empty residues')
    s=sorted(set(vs))
    return max((s[(i+1)%len(s)]-s[i])%modulus for i in range(len(s)))


def circle_resonance(n,H):
    a=H/(n*n)
    g=stepgcd(F(1,n),2*a)
    need(1/g==int(1/g),'subgroup period nonintegral')
    betas={((H/n+j*g)%1) for j in range(int(1/g))}
    worst=F(0)
    for beta in betas:
        p=2*lcm(a.denominator,beta.denominator)
        vals=[(beta*x-a*x*x)%1 for x in range(p)]
        worst=max(worst,gap(vals,F(1)))
    return a,len(betas),worst


def check_intervals():
    need(D_LO**4 < F(1,8),'lower d-bound must be rounded outward')
    need(D_HI+4*D_HI**5>1,'upper d-bound must be rounded outward')
    current=D_LO; resonance=[]
    for l,u,n,H in ROWS:
        need(l==current and u>l,'gap/overlap in Macbeath interval rows')
        need(max(abs(H-2*n*n*l**5),abs(H-2*n*n*u**5))+F(1,2)<l,
             f'Macbeath width not certified on [{l},{u}] N={n} H={H}')
        a=H/(n*n)
        if 2*l**5<=a<=2*u**5:resonance.append((n,H))
        current=u
    need(current==D_HI,'upper endpoint not covered')
    need(resonance==[(10,F(15)),(4,F(5,2)),(3,F(3,2)),(5,F(9,2)),(4,F(3))],
         f'unexpected resonance inventory: {resonance}')
    summary=[]
    for n,H in resonance:
        a,bcount,worst=circle_resonance(n,H)
        if a !=M3:
            need(worst**5<a/2,f'uncertified one-dimensional resonance a={a}, max-gap={worst}')
        else:
            need(worst==F(2,3),'incorrect exceptional resonance gap')
        summary.append({'a':str(a),'beta_classes':bcount,'gap':str(worst),'resolved':a!=M3})
    return summary


def offset(n,b,c,e,cy,f,cz,x,y,z):
    return (n*b*(x+c)-M3*(x+c)**2
            +n*e*(y+cy)-F(1,2)*(y+cy)**2
            +n*f*(z+cz)-F(1,2)*(z+cz)**2)


def n2_modular_gap(b,c,e,cy,f,cz):
    vs=[offset(2,b,c,e,cy,f,cz,x,y,z)%2
        for x,y,z in product(range(12),range(4),range(4))]
    return gap(vs,F(2))


def rect(path):
    al=F(0);au=F(1);kl=F(0);ku=F(1)
    need(isinstance(path,str) and len(path)<=24,'invalid path')
    for ch in path:
        need(ch in '0123','invalid quadtree digit')
        ia=int(ch)&1;ik=int(ch)>>1
        amid=(al+au)/2;kmed=(kl+ku)/2
        if ia:al=amid
        else:au=amid
        if ik:kl=kmed
        else:ku=kmed
    return al,au,kl,ku


def check_cover(leaves,b,c,e,cy,f,cz):
    trie={}
    for row in leaves:
        need(set(row)=={'path','witness'},'unexpected leaf structure')
        path=row['path'];ws=row['witness']
        need(isinstance(ws,list) and len(ws)==5 and all(type(z) is int for z in ws),'invalid witness')
        n,x,y,z,j=ws
        need(n in (1,-1,2,-2) and 0<=x<12 and 0<=y<4 and 0<=z<4,'witness outside authorized set')
        lo_a,hi_a,lo_k,hi_k=rect(path)
        val=offset(n,b,c,e,cy,f,cz,x,y,z)+n*j
        lo=n*n*lo_a+n*(lo_k if n>0 else hi_k)+val
        hi=n*n*hi_a+n*(hi_k if n>0 else lo_k)+val
        need(0<lo and hi<DELTA,f'witness does not strictly cover leaf {path}: [{lo},{hi}]')
        node=trie
        for ch in path:
            need('leaf' not in node,'descendant below leaf')
            node=node.setdefault(ch,{})
        need(not node,'duplicate or ancestor below existing descendant')
        node['leaf']=True
    def completeness(node):
        if 'leaf' in node:
            need(set(node)=={'leaf'},'incorrect terminal')
            return 1
        need(set(node)==set('0123'),f'incomplete quadrants at a node: {node.keys()}')
        return sum(completeness(node[k]) for k in '0123')
    count=completeness(trie)
    need(count==len(leaves),'trie count mismatch')
    return count


def check_exception(certificate_path):
    with open(certificate_path,'r') as fh:obj=json.load(fh)
    need(set(obj)=={'method','cases'},'unexpected top-level certificate structure')
    need(obj['method']=='Dyadic rational square covering, strict 0 < Q < 3/5', 'wrong proof target')
    provided={}
    for case in obj['cases']:
        need(set(case)=={'b','c','e','cy','f','cz','leaves'},'unexpected case fields')
        key=tuple(F(case[name]) for name in ('b','c','e','cy','f','cz'))
        need(key not in provided,'duplicate case')
        provided[key]=case['leaves']
    count_easy=0;count_hard=0;leaves=0
    for b,c in BC:
        # Exact derivation: exceptional beta_+=b-c/3 and beta_-=-b-c/3 belong to 1/2+(1/3)Z.
        need((3*(b-c/3)-F(3,2)).denominator==1,'bad beta_+ class')
        need((3*(-b-c/3)-F(3,2)).denominator==1,'bad beta_- class')
        for (e,cy),(f,cz) in product(EF,repeat=2):
            # Both signs of x_2 force +/- e-cy ==1/2 (mod 1), idem f,cz.
            for ee,cc in ((e,cy),(f,cz)):
                need((ee-cc-F(1,2)).denominator==1 and (-ee-cc-F(1,2)).denominator==1,
                     'bad forced dyad class')
            key=(b,c,e,cy,f,cz)
            G=n2_modular_gap(*key)
            if G<DELTA:
                need(key not in provided,'redundant supplied hard-case rectangle proof')
                need(G==F(1,3),'unexpected easy-case phase gap')
                count_easy+=1
            else:
                need(G==F(2,3),'unrecognized phase-gap structure')
                need(key in provided,'missing hard case')
                count_hard+=1
                leaves+=check_cover(provided[key],*key)
    need(len(provided)==count_hard==3,'extraneous resonance case certificates')
    need(count_easy==21,'unexpected case inventory')
    need(leaves==312,'unexpected number of exact certified leaves')
    need(DELTA**5<F(1,12),'3/5 must be strictly below d at a=1/6')
    return {'resonant_shift_classes':24,'solved_by_n2_modular_gap':count_easy,
            'solved_by_rectangle_cover':count_hard,'exact_rectangle_leaves':leaves}


def main():
    path=sys.argv[1] if len(sys.argv)>1 else 'case31_resonance_certificate.json'
    summary={'interval_pieces':len(ROWS),'resonances':check_intervals(), 'a_one_sixth':check_exception(path)}
    print(json.dumps(summary,indent=2))
    print('PASS: all finite ingredients for independent (3,1) proof verified exactly.')

if __name__=='__main__':main()
