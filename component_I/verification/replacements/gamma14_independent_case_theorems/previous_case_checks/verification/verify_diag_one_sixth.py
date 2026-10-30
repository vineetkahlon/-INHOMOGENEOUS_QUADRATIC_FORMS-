#!/usr/bin/env python3
"""Independent exact proof for the a=1/6 resonant diagonal family of Gamma(1,4).

Uses only Python stdlib and rational arithmetic. This builds and exhaustively verifies
rectangle covers; successful termination certifies a strong 0<Q<3/5 statement.
No other theorem of Bhardwaj--Kathuria--Raka is used.
"""
from fractions import Fraction as F
from itertools import product
from functools import lru_cache
import json

TARGET=F(3,5)
A=F(1,6)
EF=((F(0),F(1,2)), (F(1,2),F(0)))
BC=tuple((F(i,6),F(0)) for i in (1,3,5)) + tuple((F(i,3),F(1,2)) for i in range(3))

@lru_cache(None)
def phases(n, b, c, e, cy, f, cz):
    """Distinct offsets modulo |n| for w_i from finite, complete periods."""
    modulus=abs(n)
    # Periodicity: for each n in {1,-1,2,-2}, moving x by 12 changes terms by a multiple of |n|;
    # moving y,z by 4 changes their terms by a multiple of |n|.
    out={}
    for x, y, z in product(range(12),range(4),range(4)):
        C=n*b*(x+c)-A*(x+c)**2 + n*e*(y+cy)-F(1,2)*(y+cy)**2 + n*f*(z+cz)-F(1,2)*(z+cz)**2
        v=C % modulus
        if v not in out:out[v]=(x,y,z,C)
    return tuple((v,)+out[v] for v in sorted(out))


def maxgap(residues,modulus):
    vs=sorted(set(residues))
    return max((vs[(i+1)%len(vs)]-vs[i])%modulus for i in range(len(vs)))


def witness_for_box(box,witnesses):
    al,au,kl,ku=box
    for n,x,y,z,C in witnesses:
        lo=n*n*al + n*(kl if n>0 else ku) + C
        hi=n*n*au + n*(ku if n>0 else kl) + C
        if hi-lo>=TARGET:continue
        m=abs(n)
        q=(-lo // m)+1
        if hi+m*q <TARGET:
            j=q if n>0 else -q
            assert lo+n*j>0 and hi+n*j < TARGET
            return n,x,y,z,j
    return None


def cover_one(b,c,e,cy,f,cz, max_depth=18):
    states=[]
    for n in (1,-1,2,-2):
        for v,x,y,z,C in phases(n,b,c,e,cy,f,cz):
            states.append((n,x,y,z,C))
    # Try covering the square [0,1]^2 in alpha,kappa. We traverse a dyadic quadtree;
    # on each certified leaf, the min/max of a rational affine witness are checked exactly.
    stack=[((F(0),F(1),F(0),F(1)),0)]
    leaves=0;max_seen=0
    output=[]
    while stack:
        box,depth=stack.pop()
        w=witness_for_box(box,states)
        if w is not None:
            leaves+=1
            output.append((box,depth,w))
            continue
        max_seen=max(max_seen,depth)
        if depth>=max_depth:
            raise RuntimeError(f'Uncovered at depth {max_depth}: {box} \n one class={b,c,e,cy,f,cz}')
        l,u,v,z=box
        p=(l+u)/2; q=(v+z)/2
        stack.extend([((l,p,v,q),depth+1),((p,u,v,q),depth+1),
                      ((l,p,q,z),depth+1),((p,u,q,z),depth+1)])
    return leaves,max_seen,output


def main():
    report=[]
    total_leaves=0
    for b,c in BC:
        for (e,cy),(f,cz) in product(EF,repeat=2):
            offsets=phases(2,b,c,e,cy,f,cz)
            gap=maxgap((z[0] for z in offsets),2)
            entry={'b':str(b),'c':str(c),'e':str(e),'cy':str(cy),'f':str(f),'cz':str(cz),'n2_phase_gap':str(gap)}
            if gap < TARGET:
                entry['proof']='n=2 phase gap'
                entry['leaves']=0
            else:
                assert gap==F(2,3)
                count,depth,cover=cover_one(b,c,e,cy,f,cz)
                entry['proof']='exact dyadic rectangle cover for n=1,-1,2,-2'
                entry['leaves']=count;entry['max_depth']=depth
                total_leaves+=count
            report.append(entry)
            print(entry,flush=True)
    assert len(report)==24
    print(json.dumps({'num_classes':24, 'classes_by_n2_gap':{g:sum(e['n2_phase_gap']==g for e in report) for g in set(e['n2_phase_gap'] for e in report)}, 'total_box_leaves':total_leaves},indent=2))

if __name__=='__main__':main()
