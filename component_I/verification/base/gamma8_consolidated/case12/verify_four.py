"""Standalone exact verifier for the four residual-family lemma."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from four_families_model import DIM,DELTA,root_rows,witness_rows
from exact_linear import verify_leaf

def verify(proof):
    if proof.get('complete') is not True or proof['delta']!=str(DELTA):raise ValueError('Incomplete or wrong theorem')
    family=proof['family'];nodes=proof['nodes'];stack=[(0,root_rows())];seen=set();leaves=0;pairs=set();maxj=0
    while stack:
        idx,rows=stack.pop()
        if type(idx) is not int or not 0<=idx<len(nodes) or idx in seen:raise ValueError('Invalid tree')
        seen.add(idx);node=nodes[idx]
        if set(node)=={'certificate'}:verify_leaf(rows,node['certificate'],DIM);leaves+=1
        elif set(node)=={'n','z','j','children'}:
            n,z,j,children=(node[k] for k in ('n','z','j','children'))
            if type(n) is not int or n==0 or type(j) is not int or len(children)!=2:raise ValueError('Invalid witness')
            if len(z)!=2 or any(type(x) is not int for x in z):raise ValueError('Invalid vector')
            pairs.add((n,*z));maxj=max(maxj,abs(j))
            for child,cut in zip(children,witness_rows(family,n,z,j)):stack.append((child,rows+[cut]))
        else:raise ValueError('Unrecognized proof node')
    if len(seen)!=len(nodes):raise ValueError('Unreachable node')
    return {'family':family,'nodes':len(nodes),'leaves':leaves,'witness_pairs':len(pairs),
            'x2_values':sorted({x[0] for x in pairs}),'max_abs_j':maxj,'verified_exactly':True}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('certificate',type=Path);args=ap.parse_args()
    print(json.dumps(verify(json.loads(args.certificate.read_text())),indent=2))
