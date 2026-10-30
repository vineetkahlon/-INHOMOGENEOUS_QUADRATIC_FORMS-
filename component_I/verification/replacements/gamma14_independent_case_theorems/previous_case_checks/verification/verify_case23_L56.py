#!/usr/bin/env python3
"""Exact finite proof of BKR (m,K)=(2,3) branches L=5 and L=6.
The other L values are NOT covered. Python standard library only.
"""
import json,sys
from fractions import Fraction as Q
from verify_case32_cover import rect,require
DATA={5:(Q(19,20),Q(1173,19200),Q(3899,45000),Q(51,200),Q(4,15),82),
      6:(Q(99,100),Q(18113,336000),Q(298411,4100625),Q(59,250),Q(32,135),58)}

def check_scalars():
    for L,(D,mu_lo,mu_hi,dlo,dhi,nleaf) in DATA.items():
        require(dlo**2<Q(25,64*(L+1)), f'L{L} lower d enclosure')
        require(D<Q(15*(L+3),32)*dlo, f'L{L} witness target not below Delta_L')
        mu_min_coeff=(Q(1,3)+Q(5,16)+Q(15,8*(L+1)))/4
        require(mu_lo==mu_min_coeff*dlo,f'L{L} mu lower')
        require(mu_hi==(Q(23,24)*dhi+Q(24,5)*dhi**3)/4,f'L{L} mu upper')
        require(dhi==Q(32,15*(L+3)), f'L{L} d upper')

def check_tree(c,L):
    D,mu_lo,mu_hi,_,_,expected=DATA[L]
    require(set(c)=={'L','D','mu_lo','mu_hi','leaves'},'wrong case fields')
    require((Q(c['D']),Q(c['mu_lo']),Q(c['mu_hi']))==(D,mu_lo,mu_hi),'altered model bounds')
    trie={};mx_depth=0;mx_n=0
    for item in c['leaves']:
        require(set(item)=={'path','witness'},'wrong leaf fields')
        p=item['path']; v=item['witness']
        require(type(v) is list and len(v)==2 and all(type(k) is int for k in v),'wrong witness')
        n,j=v
        require(1<=abs(n)<=5 and abs(j)<=10,'unexpected witness bound')
        al,au,kl,ku=rect(p)
        low=n*j+n*n*al+n*(kl if n>0 else ku)-mu_hi
        high=n*j+n*n*au+n*(ku if n>0 else kl)-mu_lo
        require(Q(0)<low and high<D,f'bad box L{L} path {p} w{v}: [{low},{high}]')
        node=trie
        for ch in p:
            require('leaf' not in node,'descendant under a leaf')
            node=node.setdefault(ch,{})
        require(not node,'duplicate/ancestor of existing leaf')
        node['leaf']=1
        mx_n=max(mx_n,abs(n));mx_depth=max(mx_depth,len(p))
    def recur(node):
        if 'leaf' in node:
            require(set(node)=={'leaf'},'leaf has descendants')
            return 1
        require(set(node)==set('0123'),'incomplete dyadic partition')
        return sum(recur(node[str(k)]) for k in range(4))
    count=recur(trie)
    require(count==len(c['leaves'])==expected,'wrong leaf inventory')
    return {'L':L,'leaves':count,'max_depth':mx_depth,'max_abs_n':mx_n}

def main(path):
    check_scalars()
    with open(path) as f:z=json.load(f)
    require(set(z)=={'description','cases'} and len(z['cases'])==2,'wrong number of cases')
    seen=set();summary=[]
    for c in z['cases']:
        L=c.get('L');require(type(L) is int and L in DATA and L not in seen,'bad case')
        seen.add(L);summary.append(check_tree(c,L))
    require(seen==set(DATA),'missing branch')
    return {'status':'PASS','branch_status':'(2,3) L=5,6 only','cases':summary}
if __name__=='__main__':print(json.dumps(main(sys.argv[1] if len(sys.argv)>1 else 'case23_L56_certificate.json'),indent=2))
