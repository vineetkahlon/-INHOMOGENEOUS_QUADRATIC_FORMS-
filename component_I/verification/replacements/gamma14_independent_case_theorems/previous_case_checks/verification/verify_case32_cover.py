#!/usr/bin/env python3
"""Independent exact rational verifier of the (3,2) resonance covers.

Only the standard library is used; no floating point, LP, NumPy or scipy.
Checks four complete dyadic trees against reconstructed affine G witnesses,
and audits all rational bounds used to reduce the (3,2) case to these trees.
Run: python3 -S -B verify_case32_cover.py case32_certificate.json
"""
import json
import sys
from fractions import Fraction as Q

PARAMS = {
    (2,Q(0),Q(0)):         (Q(1,4),Q(37,50),Q(117,1600),Q(3222152,31640625)),
    (2,Q(1,4),Q(1,2)):     (Q(1,4),Q(37,50),Q(117,1600),Q(3222152,31640625)),
    (3,Q(1,6),Q(0)):       (Q(1,6),Q(16,25),Q(102,1600),Q(106036,1366875)),
    (3,Q(0),Q(1,2)):       (Q(1,6),Q(16,25),Q(102,1600),Q(106036,1366875)),
}

def require(c, message):
    if not c:raise AssertionError(message)

def check_scalars():
    """Algebraic comparisons only; radicals compared by squaring positive rationals."""
    require(Q(3,32) > (Q(32,105))**2, "L >= 4 not excluded")
    # Non-resonant Macbeath coverage L=1,2,3
    require(Q(15,64) > Q(12,25)**2,"L1 d lower")
    require((Q(15,8)*Q(12,25)) > Q(2,3), "L1 width")
    # L2 d in [sqrt(5/32), 32/75]; t>=Delta/3, t<=4d^3
    require(Q(5,32)>Q(39,100)**2, "L2 d lower")
    require(Q(15,8)*Q(39,100) > Q(9,14), "L2 small t width")
    hi=Q(32,75)
    require((16*hi**2-Q(15,8))*hi<Q(1,2), "L2 large t width")
    # L3 d in [sqrt(15/128),16/45]; t>=Delta/4, t<=4d^3
    require(Q(15,128)>Q(34,100)**2, "L3 d lower")
    require(Q(15,8)*Q(34,100)>Q(8,13), "L3 small t width")
    hi=Q(16,45)
    require(36*hi**3-1 < Q(15,8)*Q(34,100), "L3 large t width")
    bounds=[]
    for L,lo_d,hi_d,target,mu_lo in (
            (2,Q(5,32),Q(32,75),Q(37,50),Q(117,1600)),
            (3,Q(15,128),Q(16,45),Q(16,25),Q(102,1600))):
        require((target*Q(8,15))**2 < lo_d,"target not below Delta")
        min_d=Q(39,100) if L==2 else Q(34,100)
        require(min_d**2<lo_d and mu_lo==Q(3,16)*min_d,"mu lower invalid")
        mu_hi=(hi_d/3+4*L*hi_d**4)/4
        m=PARAMS[(L, Q(0),Q(0))] if L==2 else PARAMS[(L,Q(1,6),Q(0))]
        require(m==(Q(1,2*L),target,mu_lo,mu_hi),"hardcoded model bounds mismatch")
        bounds.append({'L':L,'d_upper':str(hi_d),'t':str(Q(1,2*L)),
                       'certified_width':str(target),'mu_interval':[str(mu_lo),str(mu_hi)]})
    return bounds

def rect(path):
    require(type(path) is str and len(path)<=20, "invalid path")
    a0,a1,k0,k1=Q(0),Q(1),Q(0),Q(1)
    for char in path:
        require(char in '0123',"non-quadtree branch")
        digit=int(char)
        am=(a0+a1)/2;km=(k0+k1)/2
        if digit&1:a0=am
        else:a1=am
        if digit&2:k0=km
        else:k1=km
    return a0,a1,k0,k1

def verify_tree(leaves,L,e,c,t,D,mu_lo,mu_hi):
    require(type(leaves) is list and len(leaves)>0,"empty certificate")
    trie={};maxdepth=0; max_n=max_z=max_j=0
    for item in leaves:
        require(type(item) is dict and set(item)=={'path','witness'},"bad leaf fields")
        p=item['path'];w=item['witness']
        require(type(w) is list and len(w)==3 and all(type(x) is int for x in w),"bad witness")
        n,z,j=w
        require(1<=abs(n)<=5 and abs(z)<=8 and abs(j)<=40,"witness outside agreed search set")
        al,au,kl,ku=rect(p)
        offset=n*j+n*e*(z+c)-t*(z+c)**2
        # Since mu occurs with coefficient -1, lower/upper are attained at mu_hi/mu_lo.
        lo=n*n*al+n*(kl if n>0 else ku)+offset-mu_hi
        hi=n*n*au+n*(ku if n>0 else kl)+offset-mu_lo
        require(Q(0)<lo and hi<D, f"FAIL: L={L}, e={e}, c={c}, path={p}, witness={w}, value interval=[{lo},{hi}]")
        maxdepth=max(maxdepth,len(p));max_n=max(max_n,abs(n));max_z=max(max_z,abs(z));max_j=max(max_j,abs(j))
        current=trie
        for char in p:
            require('leaf' not in current,"descendant beneath leaf")
            current=current.setdefault(char,{})
        require(not current,"duplicate leaf or ancestor above descendants")
        current['leaf']=True
    def visit(node):
        if 'leaf' in node:
            require(set(node)=={'leaf'},"leaf contains descendants")
            return (0,1)
        require(set(node)==set('0123'),f"incomplete quadtree: {set(node)}")
        rows=[visit(node[str(k)]) for k in range(4)]
        return (1+sum(r[0] for r in rows), sum(r[1] for r in rows))
    intern,leaf_ct=visit(trie)
    require(leaf_ct==len(leaves),"unreachable certificate row")
    return {'L':L,'e':str(e),'c':str(c),'internal':intern,'leaves':leaf_ct,'maxdepth':maxdepth,
            'max_abs_n':max_n,'max_abs_z':max_z,'max_abs_j':max_j}

def check_reps():
    # The resonance condition for n=+/-1 has +/-e - c/L == 1/2 (mod 1/L).
    # An integer shear can change e by an integral multiple of 1/L.
    for L,e,c in PARAMS:
        for sign in (1,-1):
            require((L*(sign*e-c/L-Q(1,2))).denominator==1, "invalid forced Macbeath congruence")
        require(2*c in (0,1),"bad c shift class")
    for L in (2,3):
        pairs={ (e,c) for k,e,c in PARAMS if k==L }
        # For each c=0 or 1/2, e is prescribed modulo 1/L by the two signs.
        for c in (Q(0),Q(1,2)):
            e0=(Q(1,2)+c/L) % Q(1,L)
            require((e0,c) in pairs,"missing shear-reduced resonance class")
        require(len(pairs)==2,"extra resonance class")

def main(path):
    bounds=check_scalars();check_reps()
    with open(path,encoding='utf-8') as f: data=json.load(f)
    require(type(data) is dict and set(data)=={'description','cases'},"bad top level")
    require(type(data['cases']) is list and len(data['cases'])==4,"wrong number of cases")
    seen=set();summary=[]
    for obj in data['cases']:
        require(set(obj)=={'L','t','D','mu_lo','mu_hi','e','c','leaves'},"wrong model fields")
        L=obj['L'];require(type(L) is int,"noninteger L")
        e,c=Q(obj['e']),Q(obj['c']);key=(L,e,c)
        require(key in PARAMS and key not in seen,"missing/duplicate/foreign family")
        seen.add(key)
        t,D,ml,mh=PARAMS[key]
        require(tuple(Q(obj[k]) for k in ('t','D','mu_lo','mu_hi'))==(t,D,ml,mh),"model parameters were modified")
        summary.append(verify_tree(obj['leaves'],L,e,c,t,D,ml,mh))
    require(seen==set(PARAMS),"missing required resonance cases")
    expected={(2,Q(0),Q(0)):40,(2,Q(1,4),Q(1,2)):40,(3,Q(1,6),Q(0)):70,(3,Q(0),Q(1,2)):52}
    for s in summary:
        require(s['leaves']==expected[(s['L'],Q(s['e']),Q(s['c']))],"unexpected certificate inventory")
    return {'status':'PASS','rational_bounds':bounds,'certificate_families':summary,
            'total_terminal_rectangles':sum(x['leaves'] for x in summary),
            'total_internal_nodes':sum(x['internal'] for x in summary)}

if __name__=='__main__':
    if len(sys.argv)>2:raise SystemExit('usage: verify_case32_cover.py [certificate.json]')
    path=sys.argv[1] if len(sys.argv)>1 else 'case32_certificate.json'
    answer=main(path)
    print(json.dumps(answer,indent=2))
