"""Fixed rational model for the four residual (1,2,1) shift families."""
from fractions import Fraction as F
DIM=3
DELTA=F(7,12)

def root_rows():
    rows=[]
    for i,(lo,hi) in enumerate([(F(29,100),F(1,3)),(F(-1,2),F(1,2)),(F(-1,2),F(1,2))]):
        row=[F(0)]*DIM;row[i]=1;rows.append((row,hi))
        row=[F(0)]*DIM;row[i]=-1;rows.append((row,-lo))
    return rows

def witness_affine(family,n,z,j=0):
    if family not in (0,1,2,3):raise ValueError('Bad family')
    i,k=divmod(family,2)
    c4=F(i,2);ell4=F(i,4)
    c5=F(k,2);ell5=F(1-k,2)
    y,w=z
    tau=c4*c4/4+c5*c5/2
    constant=n*j-tau+n*ell4*y+n*ell5*w-y*y*F(1,4)-w*w*F(1,2)-c4*y/2-c5*w
    return [F(-1,4),F(n*n),F(n)],constant

def witness_rows(family,n,z,j):
    row,offset=witness_affine(family,n,z,j)
    return [(row,-offset),([-x for x in row],offset-DELTA)]
