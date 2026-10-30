"""The four Section 5 geometries and all forced shift/linear families."""
from fractions import Fraction as F
from itertools import product

DIM=2
PAIRS=((F(0),F(1,2)),(F(1,2),F(0)))
DELTAS={1:F(151,200),2:F(87,100),3:F(15,16),'hex':F(15,16)}


def families():
    out=[]
    for r in (1,2,3):
        last=[]
        for b in (F(0),F(1,2)):
            for j in range(r):
                c=((b-F(r,2)+j)/r)%1
                last.append((b,c))
        for pairs in product(PAIRS,PAIRS,last):
            out.append({'kind':r,'b':tuple(p[0] for p in pairs),
                        'c':tuple(p[1] for p in pairs)})
    for b0,c0 in PAIRS:
        for b1,b2 in product((F(0),F(1,2)),repeat=2):
            for j in range(3):
                c1=((2*b1-b2+2*j)/3)%1
                c2=((-b1+2*b2-j)/3)%1
                out.append({'kind':'hex','b':(b0,b1,b2),'c':(c0,c1,c2)})
    return out


FAMILIES=families()


def phi(kind,w):
    x,y,z=w
    if kind=='hex':return x*x/2+y*y+y*z+z*z
    return (x*x+y*y+kind*z*z)/2


def root_rows():
    return [([F(1),F(0)],F(1,2)),([F(-1),F(0)],F(1,2)),
            ([F(0),F(1)],F(1,2)),([F(0),F(-1)],F(1,2))]


def witness_affine(index,n,w,j=0):
    family=FAMILIES[index]
    shift=tuple(F(x)+c for x,c in zip(w,family['c']))
    offset=n*j+n*sum(b*x for b,x in zip(family['b'],w))-phi(family['kind'],shift)
    return [F(n*n),F(n)],offset


def witness_rows(index,n,w,j):
    row,offset=witness_affine(index,n,w,j)
    delta=DELTAS[FAMILIES[index]['kind']]
    return [(row,-offset),([-x for x in row],offset-delta)]
