"""Exact checker for 0<f<=d. Failure branches retain their strict endpoints."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from binary_resonance_model import DIM,root_rows,witness_rows

def verify_leaf(rows,weights):
    coefficients=[F(0)]*DIM;rhs=F(0);strict_mass=F(0)
    for i,raw in weights:
        if type(i) is not int or not 0<=i<len(rows):raise ValueError('Invalid row index')
        weight=F(raw)
        if weight<0:raise ValueError('Negative multiplier')
        a,b,strict=rows[i]
        coefficients=[c+weight*x for c,x in zip(coefficients,a)]
        rhs+=weight*b
        if strict:strict_mass+=weight
    if any(coefficients):raise ValueError('Nonzero coefficient sum')
    if not (rhs<0 or rhs==0 and strict_mass>0):raise ValueError('Invalid strict/non-strict contradiction')
    return rhs==0

def verify(proof):
    if proof.get('complete') is not True or proof.get('target')!='0<f<=d':raise ValueError('Incomplete or wrong target')
    if proof.get('model_version')!=2:raise ValueError('Wrong model version')
    family=proof['family']
    if family not in (0,1):raise ValueError('Bad family')
    root=[(a,b,False) for a,b in root_rows(proof.get('shape_box'))]
    nodes=proof['nodes'];stack=[(0,root)];seen=set();leaves=zeros=0;zs=set();maxj=0
    while stack:
        idx,rows=stack.pop()
        if type(idx) is not int or not 0<=idx<len(nodes) or idx in seen:raise ValueError('Invalid tree')
        seen.add(idx);node=nodes[idx]
        if set(node)=={'certificate'}:
            zeros+=verify_leaf(rows,node['certificate']);leaves+=1
        elif set(node)=={'n','z','j','children'}:
            n,z,j,children=node['n'],node['z'],node['j'],node['children']
            if type(n) is not int or not n:raise ValueError('Bad n')
            if len(z)!=2 or any(type(x) is not int for x in z) or type(j) is not int or len(children)!=2:raise ValueError('Invalid witness')
            zs.add((n,*z));maxj=max(maxj,abs(j));cuts=witness_rows(family,n,z,j)
            for k in (0,1):stack.append((children[k],rows+[(*cuts[k],bool(k))]))
        else:raise ValueError('Unrecognized node')
    if len(seen)!=len(nodes):raise ValueError('Unreachable nodes')
    return {'family':family,'nodes':len(nodes),'leaves':leaves,'strict_zero_leaves':zeros,'vectors':len(zs),'x2_values':sorted({z[0] for z in zs}),'max_abs_j':maxj,'verified_exactly':True}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path);args=ap.parse_args()
    print(json.dumps(verify(json.loads(args.certificate.read_text())),indent=2))
