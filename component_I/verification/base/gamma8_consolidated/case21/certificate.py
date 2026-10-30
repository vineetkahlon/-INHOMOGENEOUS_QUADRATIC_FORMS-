"""Proof interface and exact, closed-branch checker. Standard library only."""
from fractions import Fraction as F
import ternary_model as T
import diagonal_model as D
from exact_linear import verify_leaf

EXPECTED={
    'low_plus.json':{'kind':'ternary','sign':1},
    'low_minus.json':{'kind':'ternary','sign':-1},
    **{f'diagonal_{i}.json':{'kind':'diagonal','family':i} for i in range(3)}
}

def model(domain):
    if domain=={'kind':'ternary','sign':1} or domain=={'kind':'ternary','sign':-1}:
        if type(domain['sign']) is not int:raise ValueError('Invalid sign')
        return T,T.root_rows(domain['sign'])
    if set(domain)=={'kind','family'} and domain['kind']=='diagonal':
        D.family_parameters(domain['family'])
        return D,D.root_rows()
    raise ValueError('Unknown domain')

def affine(domain,witness):
    if domain['kind']=='ternary':
        z,j=witness['z'],witness['j']
        return T.witness_affine(z,j)
    return D.witness_affine(domain['family'],witness['n'],witness['z'],witness['j'])

def cuts(domain,witness):
    row,c=affine(domain,witness)
    upper=[-x for x in row];upper[0]+=1
    # 0 < f < d fails exactly when f <= 0 or f >= d.
    return [(row,-c),(upper,c)]

def verify(proof):
    if proof.get('complete') is not True or proof.get('target')!='0<f<d':raise ValueError('Incomplete or wrong target')
    if proof.get('model_version')!=1:raise ValueError('Wrong model version')
    domain=proof['domain'];module,root=model(domain)
    nodes=proof['nodes'];stack=[(0,root)];seen=set();leaves=0;zs=set();maxj=0;ns=set();maxcoord=0
    while stack:
        idx,rows=stack.pop()
        if type(idx) is not int or not 0<=idx<len(nodes) or idx in seen:raise ValueError('Invalid tree')
        seen.add(idx);node=nodes[idx]
        if set(node)=={'certificate'}:
            verify_leaf(rows,node['certificate'],module.DIM);leaves+=1
        else:
            required={'z','j','children'}|({'n'} if domain['kind']=='diagonal' else set())
            if set(node)!=required:raise ValueError('Unrecognized node')
            z,j,children=node['z'],node['j'],node['children']
            if len(z)!=3 or any(type(x) is not int for x in z) or type(j) is not int or len(children)!=2:raise ValueError('Invalid witness')
            n=node.get('n',1)
            if type(n) is not int or not n:raise ValueError('Invalid x2')
            ns.add(n);zs.add((n,*z));maxj=max(maxj,abs(j));maxcoord=max(maxcoord,*map(abs,z))
            for child,cut in zip(children,cuts(domain,node)):stack.append((child,rows+[cut]))
    if len(seen)!=len(nodes):raise ValueError('Unreachable nodes')
    return {'domain':domain,'nodes':len(nodes),'leaves':leaves,'witness_vectors':len(zs),'x2_values':sorted(ns),
            'max_abs_j':maxj,'max_abs_integer_negative_coordinate':maxcoord,'strict_upper_bound':True,'verified_exactly':True}
