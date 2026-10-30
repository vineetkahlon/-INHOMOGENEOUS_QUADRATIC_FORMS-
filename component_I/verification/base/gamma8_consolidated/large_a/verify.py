"""Exact verification of the Section 5 family-cover lemma."""
import json
import hashlib
from itertools import product
from pathlib import Path
from fractions import Fraction as F
from model import DIM,FAMILIES,DELTAS,root_rows,witness_rows
from exact_linear import verify_leaf


def require(x,msg):
    if not x:raise ValueError(msg)


def verify(proof):
    index=proof['family']
    require(type(index) is int and 0<=index<len(FAMILIES),'Bad family')
    require(proof.get('complete') is True and not proof.get('pending'),'Incomplete proof')
    require(proof['delta']==str(DELTAS[FAMILIES[index]['kind']]),'Bad target')
    nodes=proof['nodes'];stack=[(0,root_rows())];seen=set();leaves=0;ns=set();maxw=0
    while stack:
        idx,rows=stack.pop()
        require(type(idx) is int and 0<=idx<len(nodes) and idx not in seen,'Bad tree')
        seen.add(idx);node=nodes[idx]
        if set(node)=={'certificate'}:
            verify_leaf(rows,node['certificate'],DIM);leaves+=1
        elif set(node)=={'n','w','j','children'}:
            n,w,j,children=(node[k] for k in ('n','w','j','children'))
            require(type(n) is int and n!=0 and type(j) is int,'Bad integers')
            require(len(w)==3 and all(type(x) is int for x in w),'Bad ternary witness')
            require(len(children)==2,'Missing branch')
            ns.add(n);maxw=max(maxw,*map(abs,w))
            cuts=witness_rows(index,n,w,j)
            for child,cut in zip(children,cuts):stack.append((child,rows+[cut]))
        else:raise ValueError('Unknown node')
    require(len(seen)==len(nodes),'Unreachable nodes')
    return {'family':index,'kind':FAMILIES[index]['kind'],'nodes':len(nodes),
            'leaves':leaves,'x2_values':sorted(ns),'max_abs_w':maxw,'verified_exactly':True}


def run():
    root=Path(__file__).resolve().parent
    require(len(FAMILIES)==72,'Wrong family count')
    require(len({(f['kind'],f['b'],f['c']) for f in FAMILIES})==72,'Repeated family')
    # Independently enumerate the rational solutions of the failure congruences.
    expected=set()
    first=((F(0),F(1,2)),(F(1,2),F(0)))
    for r in (1,2,3):
        for (b0,c0),(b1,c1),b2,j in product(first,first,(F(0),F(1,2)),range(2*r)):
            c2=F(j,2*r)
            if all((sign*b2-r*c2-F(r,2)).denominator==1 for sign in (-1,1)):
                expected.add((r,(b0,b1,b2),(c0,c1,c2)))
    for (b0,c0),b1,b2,j,k in product(first,(F(0),F(1,2)),(F(0),F(1,2)),range(6),range(6)):
        c1,c2=F(j,6),F(k,6)
        if (2*c1+c2-b1).denominator==1 and (c1+2*c2-b2).denominator==1:
            expected.add(('hex',(b0,b1,b2),(c0,c1,c2)))
    require(expected=={(f['kind'],f['b'],f['c']) for f in FAMILIES},'Incomplete congruence cover')
    for r in (1,2,3):
        require(DELTAS[r]**5<F(r,4),'Target is not below d')
    require(DELTAS['hex']**5<F(3,4),'Hex target is not below d')
    # Exact comparisons in the condensed Section 4 proof.
    for m,k_lower in ((4,3),(3,2),(2,2)):
        require(F(3*(m+3)**4,128*(m+1))>k_lower**2,'Incorrect K lower bound')
    for m,K in ((4,3),(3,3),(2,5)):
        require(F((m+3)**4*(K+3)**3,(m+1)*(K+1))>4**7,'Jackson threshold failed')
    require(F(675,2048)>F(1,4),'(3,2) width comparison failed')
    require(F(78125,18432)>4,'L lower bound failed')
    require(F(175,192)**2/3>F(1,4),'Last width comparison failed')
    summaries=[]
    for i in range(72):
        path=root/'certificates'/f'family_{i:02d}.json'
        raw=path.read_bytes();proof=json.loads(raw)
        require(proof['family']==i,'Missing required family')
        summaries.append({'file':str(path.relative_to(root)),
                          'sha256':hashlib.sha256(raw).hexdigest(),**verify(proof)})
    return {'result':'All 72 forced families and Section 4 scalar comparisons verified exactly.',
            'nodes':sum(s['nodes'] for s in summaries),
            'leaves':sum(s['leaves'] for s in summaries),'families':summaries}


if __name__=='__main__':print(json.dumps(run(),indent=2))
