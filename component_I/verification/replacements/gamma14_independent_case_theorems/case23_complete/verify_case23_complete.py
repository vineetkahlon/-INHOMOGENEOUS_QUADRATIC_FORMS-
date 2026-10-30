#!/usr/bin/env python3
"""Exact, standard-library verification for (m,K)=(2,3), all L=1,...,6.

There are two certificate formats:
  * Seven 2D dyadic quadtree affine F/G covers at the rational resonances,
  * Two 2D H covers for L=5,6 from the 2026-10-26 supplement.

Five uniform 1D quadratic gap lemmas are verified via an exact arrangement
of rational lines in the (beta,t) plane. These are analytic finite proof
steps, not numerical sampling. The checker also verifies stated rational
inequalities reducing all parameter domains to these finite models.

IMPORTANT: This program only verifies the explicitly reconstructed finite
claims; the mathematical arguments connecting them to the canonical form,
Dumir's lemma and Macbeath's lemma must also be reviewed.
"""
from fractions import Fraction as F
from itertools import combinations
import json, sys, pathlib

BASE=pathlib.Path(__file__).resolve().parent

def require(condition,msg):
    if not condition:raise AssertionError(msg)

def rect(path):
    require(type(path) is str and len(path)<=25,'invalid quadtree path')
    a0,a1,k0,k1=F(0),F(1),F(0),F(1)
    for ch in path:
        require(ch in '0123','invalid digit in quadtree path')
        digit=int(ch); am=(a0+a1)/2;km=(k0+k1)/2
        a0,a1=(am,a1) if digit&1 else (a0,am)
        k0,k1=(km,k1) if digit&2 else (k0,km)
    return a0,a1,k0,k1

def exact_arrangement_max_gap(tlo,thi,zmax=4):
    """Find global maximum circular gap of {beta*z-t*z*z mod 1}.

    beta in [0,1/2], t in [tlo,thi]. Symmetry z -> -z handles beta<0.
    Collision lines are (i-j)*beta-(i*i-j*j)*t=k, k integer.
    Between such lines the order on the unit circle is fixed; every
    circular gap is affine, hence attains its maximum at a cell vertex.
    The vertices are exactly intersections of the collision and four
    boundary lines. Everything here is computed as exact fractions.
    """
    require(F(0)<tlo<thi<F(1),'bad t interval')
    bl,bh=F(0),F(1,2)
    zrange=tuple(range(-zmax,zmax+1))
    lines=set()
    for i,j in combinations(zrange,2):
        a=i-j; b=j*j-i*i
        values=(a*bl+b*tlo,a*bl+b*thi,a*bh+b*tlo,a*bh+b*thi)
        lo=min(values);hi=max(values)
        for k in range(lo.numerator//lo.denominator-1,
                       hi.numerator//hi.denominator+2):
            if lo<=k<=hi:lines.add((a,b,F(k)))
    lines=sorted(lines)+[(1,0,bl),(1,0,bh),(0,1,tlo),(0,1,thi)]
    points=set()
    for (a,b,k),(c,d,l) in combinations(lines,2):
        determinant=a*d-b*c
        if determinant==0:continue
        beta=(k*d-b*l)/F(determinant)
        t=(a*l-k*c)/F(determinant)
        if bl<=beta<=bh and tlo<=t<=thi:points.add((beta,t))
    require(len(points)>=4,'not enough arrangement vertices')
    mx=F(0);attainers=[]
    for beta,t in points:
        residues=sorted({(beta*z-t*z*z)%1 for z in zrange})
        gaps=[residues[q+1]-residues[q] for q in range(len(residues)-1)]
        gaps.append(F(1)+residues[0]-residues[-1])
        gap=max(gaps)
        if gap>mx:mx=gap;attainers=[(beta,t)]
        elif gap==mx:attainers.append((beta,t))
    return mx,len(lines),len(points),attainers

GAP_BOXES={
    'L1_A_low':  (F(7,50),F(157,1000),F(57,100),F(57,100)),
    'L1_A_high': (F(22,125),F(1,5),F(3,5),F(3,5)),
    'L3_t_low':  (F(73,500),F(63,400),F(23,40),F(23,40)),
    'L3_t_high': (F(22,125),F(27,125),F(3,5),F(3,5)),
    'L4_t_full': (F(13,125),F(7,50),F(1,2),F(1,2))
}

def check_gap_boxes():
    results=[]
    for name,(lo,hi,expected,allowed) in GAP_BOXES.items():
        gap,nlines,nvertices,where=exact_arrangement_max_gap(lo,hi)
        require(gap==expected and gap<=allowed,'wrong exact maximum at '+name)
        results.append({'box':name,'max_circular_gap':str(gap),
                        'collision_and_boundary_lines':nlines,
                        'exact_arrangement_vertices':nvertices,
                        'maximizers':[list(map(str,v)) for v in where[:5]]})
    return results

# New resonance models, two families each at L=1,2,3.
NEW={
 (1,F(1,6),F(0)):(F(1,6),F(3,5),F(1,25),F(53,1000),106),
 (1,F(0),F(1,2)):(F(1,6),F(3,5),F(1,25),F(53,1000),88),
 (2,F(0),F(0)):(F(1,4),F(7,10),F(31,192)*F(28,75),F(23,96)*F(32,75),76),
 (2,F(1,4),F(1,2)):(F(1,4),F(7,10),F(31,192)*F(28,75),F(23,96)*F(32,75),88),
 (3,F(1,6),F(0)):(F(1,6),F(61,100),F(31,192)*F(163,500),F(7,100),124),
 (3,F(0),F(1,2)):(F(1,6),F(61,100),F(31,192)*F(163,500),F(7,100),64)
}

def check_tree(leaves,known,L,e,c,t,D,mlo,mhi,kind):
    require(type(leaves) is list and leaves,'empty leaves')
    trie={};depth=0;nmax=zmax=jmax=0
    for item in leaves:
        require(type(item) is dict and set(item)=={'path','witness'},'bad certificate leaf')
        path=item['path'];w=item['witness']
        require(type(w) is list and all(type(q) is int for q in w),'noninteger witness')
        if kind=='H':
            require(len(w)==2,'H witness dimension')
            n,j=w;z=None
        else:
            require(len(w)==3,'FG witness dimension')
            n,z,j=w
        require(1<=abs(n)<=5 and abs(j)<=40 and (z is None or abs(z)<=8),'witness out of range')
        al,au,kl,ku=rect(path)
        off=n*j if z is None else n*j+n*e*(z+c)-t*(z+c)**2
        lo=off+n*n*al+n*(kl if n>0 else ku)-mhi
        hi=off+n*n*au+n*(ku if n>0 else kl)-mlo
        require(F(0)<lo and hi<D, f'{kind} invalid cover L{L} e{e} c{c} p{path} w{w} lo{lo} hi{hi}')
        depth=max(depth,len(path));nmax=max(nmax,abs(n));jmax=max(jmax,abs(j))
        if z is not None:zmax=max(zmax,abs(z))
        node=trie
        for ch in path:
            require('leaf' not in node,'descendant beneath existing leaf')
            node=node.setdefault(ch,{})
        require(not node,'duplicate or ancestor leaf')
        node['leaf']=1
    def visit(node):
        if 'leaf' in node:
            require(set(node)=={'leaf'},'malformed terminal node');return 0,1
        require(set(node)==set('0123'),'incomplete quadtree partition')
        sub=[visit(node[x]) for x in '0123']
        return 1+sum(t[0] for t in sub),sum(t[1] for t in sub)
    intern,leaves_count=visit(trie)
    require(leaves_count==len(leaves)==known,'wrong leaf inventory')
    return {'L':L,'type':kind,'e':str(e) if e is not None else None,
            'c':str(c) if c is not None else None,'leaves':leaves_count,
            'internal':intern,'max_depth':depth,'max_abs_n':nmax,
            'max_abs_z':zmax,'max_abs_j':jmax}

def check_new(path):
    with open(path) as f: data=json.load(f)
    require(set(data)=={'description','cases'} and len(data['cases'])==6,'new inventory fields')
    seen=set();out=[]
    for case in data['cases']:
        require(set(case)=={'L','t','D','mu_lo','mu_hi','e','c','leaves'},'new case fields')
        L=case['L'];require(type(L) is int,'bad index')
        e,c=F(case['e']),F(case['c']); key=(L,e,c)
        require(key in NEW and key not in seen,'duplicate, missing or extraneous case')
        seen.add(key)
        t,D,mlo,mhi,nleaves=NEW[key]
        require(tuple(F(case[x]) for x in ('t','D','mu_lo','mu_hi'))==(t,D,mlo,mhi),'parameter tampering')
        kind='F' if L==1 else 'G'
        out.append(check_tree(case['leaves'],nleaves,L,e,c,t,D,mlo,mhi,kind))
    require(seen==set(NEW),'missing resonance case')
    return out

OLD={
    5:(F(19,20),F(1173,19200),F(3899,45000),F(51,200),F(4,15),82),
    6:(F(99,100),F(18113,336000),F(298411,4100625),F(59,250),F(32,135),58)
}
def check_first_L1_resonance(path):
    with open(path) as f: data=json.load(f)
    require(set(data)=={'description','case'},'L1-t fields')
    case=data['case']
    require(set(case)=={'L','t','D','mu_lo','mu_hi','e','c','leaves'},'L1-t case fields')
    require(case['L']==1 and F(case['e'])==F(0) and F(case['c'])==F(1,2),'wrong L1-t family')
    D=F(22,25);mlo=F(31,192)*F(47,100);mhi=F(23,96)*F(8,15)
    require(tuple(F(case[x]) for x in ('t','D','mu_lo','mu_hi'))==(F(1,2),D,mlo,mhi),'L1-t bounds altered')
    return check_tree(case['leaves'],124,1,F(0),F(1,2),F(1,2),D,mlo,mhi,'G')

def check_old(path):
    with open(path) as f:data=json.load(f)
    require(set(data)=={'description','cases'} and len(data['cases'])==2,'old inventory fields')
    seen=set();out=[]
    for case in data['cases']:
        require(set(case)=={'L','D','mu_lo','mu_hi','leaves'},'old case fields')
        L=case['L'];require(type(L) is int and L in OLD and L not in seen,'bad old index')
        seen.add(L)
        D,mlo,mhi,dlo,dhi,nleaves=OLD[L]
        require(tuple(F(case[x]) for x in ('D','mu_lo','mu_hi'))==(D,mlo,mhi),'old model tampering')
        out.append(check_tree(case['leaves'],nleaves,L,None,None,None,D,mlo,mhi,'H'))
    require(seen==set(OLD),'missing L5/6')
    return out

def check_scalars():
    # Original pair: d/3<=a<d/2; delta=d+3a/4;
    # delta/4<=A<delta/3; Delta=delta+2A;
    # Delta>=15d/8, t<=24d^3/5, Delta_L>=15(L+3)d/32.
    # Cases with Delta>=1, or Delta_L>=1, are already settled by rounding.
    require(F(47,100)**3<F(5,48),'L1 t=1/2 implies d>47/100')
    require(F(15,8)*F(47,100)>F(22,25),'L1-t G target below Delta')
    require(F(5,16)*F(47,100)>F(7,50),'L1 A lower outside gap band')
    delta_lo=F(5,4)*F(47,100)
    require(delta_lo>F(57,100),'L1 low A circle gap target')
    require(F(1,2)+9*(F(1,6)-F(157,1000))<delta_lo,'L1 middle left Macbeath')
    require(F(1,2)+9*(F(22,125)-F(1,6))<delta_lo,'L1 middle right Macbeath')
    require((F(12,25))**4<F(22,375),'L1 high A yields d>12/25, delta>3/5')
    require(F(97,200)**4<F(1,18),'L1 A=1/6 yields d>97/200')
    df=F(511,1000)
    require(df+F(9,2)*df**5>F(2,3),'L1 A=1/6 yields d<511/1000')
    require(F(1,25)<F(3,2)*F(97,200)**5,'L1 F mu minimum')
    require(F(3,2)*df**5<F(53,1000),'L1 F mu maximum')
    require(F(3,5)<F(97,200)+F(9,2)*F(97,200)**5,'L1 F target below delta')
    # L2: d^2>=25/192; Delta >=15d/8 >9/14.
    require((F(8,15)*F(9,14))**2 <F(25,192),'L2 Delta>9/14')
    require(F(28,75)**3 < F(5,96),'L2 t>=1/4 yields d>28/75 and Delta>7/10')
    require(F(15,8)*F(28,75)==F(7,10),'L2 target via strict d>28/75')
    require(F(31,192)*F(28,75) < F(23,96)*F(32,75),'L2 mu interval')
    # L3 d^2 >=25/256 so d>=5/16; t in [73/500,27/125].
    require(F(15,32)*F(5,16)>F(73,500),'L3 t below lower band')
    require(F(24,5)*F(16,45)**3<F(27,125),'L3 t above upper band')
    delta3_lo=F(15,8)*F(5,16)
    require(delta3_lo>F(23,40),'L3 low t circle gap')
    require(F(1,2)+abs(F(3,2)-9*F(63,400))<delta3_lo,'L3 middle left')
    require(F(1,2)+abs(F(3,2)-9*F(22,125))<delta3_lo,'L3 middle right')
    require(F(11,300)>F(8,25)**3,'L3 high t implies d>8/25,Delta>3/5')
    require(F(15,8)*F(8,25)==F(3,5),'L3 high t target')
    require(F(163,500)**3<F(5,144),'L3 t=1/6 d lower')
    require(F(15,8)*F(163,500)>F(61,100),'L3 resonance G target')
    d3=F(16,45)
    muupper=max((d3/3+9*d3**4)/4,(d3/2+6*d3**4)/4)
    require(muupper<F(7,100),'L3 convex mu upper')
    # L4 d^2>=25/320=5/64, d<32/105.
    require((F(13,125)*F(8,3))**2<F(5,64),'L4 t lower')
    require(F(24,5)*F(32,105)**3<F(7,50),'L4 t upper')
    require((F(8,15)*F(1,2))**2<F(5,64),'L4 Delta>1/2')
    # L5,L6 old package bounds including exact model mu values.
    for L,(D,mlo,mhi,dlo,dhi,expected) in OLD.items():
        require(dlo**2 < F(25,64*(L+1)),'L5/6 d lower rounded outward')
        require(F(15*(L+3),32)*dlo>D,'L5/6 H target')
        coeff=(F(1,3)+F(5,16)+F(15,8*(L+1)))/4
        require(mlo==coeff*dlo,'L5/6 mu lower')
        require(mhi==(F(23,24)*dhi+F(24,5)*dhi**3)/4,'L5/6 mu upper')
        require(dhi==F(32,15*(L+3)),'L5/6 d upper')
    return {'status':'PASS','scope':'rational bounds used in all L=1..6 reductions'}

def main():
    scalars=check_scalars()
    gaps=check_gap_boxes()
    new=check_new(BASE/'case23_new_certificates.json')
    extra=check_first_L1_resonance(BASE/'case23_L1_t_certificate.json')
    new=[extra]+new
    old=check_old(BASE/'case23_L56_certificate.json')
    return {'status':'PASS','scope':'(m,K)=(2,3), L=1,...,6',
            'scalar_comparisons':scalars,
            'uniform_gap_lemmas':gaps,
            'resonance_certificates':new,
            'third_completion_certificates':old,
            'totals':{'exact_gap_boxes':len(gaps),'quadtree_certificates':len(new)+len(old),
                      'terminal_rectangles':sum(x['leaves'] for x in new+old),
                      'internal_nodes':sum(x['internal'] for x in new+old)}}

if __name__=='__main__':
    print(json.dumps(main(),indent=2))
