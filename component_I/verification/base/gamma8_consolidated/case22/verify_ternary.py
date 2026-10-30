"""Independent rational checker for the ternary outer-polytope lemma."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from model import DIM,root_rows,witness_rows

def verify_leaf(rows,certificate):
    coeff=[F(0)]*DIM;constant=F(0)
    for index,raw in certificate:
        if type(index) is not int or not 0<=index<len(rows):raise ValueError('Bad row index')
        w=F(raw)
        if w<0:raise ValueError('Negative multiplier')
        a,b=rows[index]
        coeff=[c+w*x for c,x in zip(coeff,a)];constant+=w*b
    if any(coeff) or constant>=0:raise ValueError('Invalid Farkas identity')

def verify(proof):
    if proof.get('complete') is not True:raise ValueError('Incomplete search is not a proof')
    nodes=proof['nodes'];seen=set();stack=[(0,root_rows(proof['L'],proof.get('model_version',1),proof.get('d_interval')))]
    leaves=0;zs=set();max_j=0
    while stack:
        idx,rows=stack.pop()
        if type(idx) is not int or not 0<=idx<len(nodes) or idx in seen:raise ValueError('Bad tree')
        seen.add(idx);node=nodes[idx]
        if set(node)=={'certificate'}:
            verify_leaf(rows,node['certificate']);leaves+=1
        elif set(node)=={'z','j','children'}:
            z,j,children=node['z'],node['j'],node['children']
            if len(z)!=3 or any(type(x) is not int for x in z) or type(j) is not int or len(children)!=2:
                raise ValueError('Bad integer witness')
            zs.add(tuple(z));max_j=max(max_j,abs(j))
            cuts=witness_rows(z,j)
            for child,cut in zip(children,cuts):stack.append((child,rows+[cut]))
        else:raise ValueError('Unrecognized proof node')
    if len(seen)!=len(nodes):raise ValueError('Unreachable node')
    return {'L':proof['L'],'nodes':len(nodes),'leaves':leaves,'vectors':len(zs),'max_abs_j':max_j,'verified_exactly':True}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path);args=ap.parse_args()
    print(json.dumps(verify(json.loads(args.certificate.read_text())),indent=2))
