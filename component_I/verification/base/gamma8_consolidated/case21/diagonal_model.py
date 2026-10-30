"""Uniform diagonal square-block model after Macbeath has forced h=k=0."""
from fractions import Fraction as F
NAMES=('d','a','u','v','b','alpha','kappa')
DIM=len(NAMES)
def vector(**kw):
    row=[F(0)]*DIM
    for k,v in kw.items():row[NAMES.index(k)]=F(v)
    return row

def root_rows():
    out=[]
    def add(rhs=0,**kw):out.append((vector(**kw),F(rhs)))
    def bound(n,ll,hh):add(hh,**{n:1});add(-ll,**{n:-1})
    lo,hi=F(6389,10000),F(7,10)
    bound('d',lo,hi)
    add(d=F(1,3),a=-1);add(a=1,d=F(-1,2));add(1,d=1,a=1)
    add(F(5,16),a=1)
    add(u=-1);add(u=1,a=F(-1,2));add(v=1,u=F(-1,2))
    for xi in [F(j,8) for j in range(5)]:add(a=-xi*xi,u=2*xi,v=-1)
    for name in ['b','alpha','kappa']:bound(name,F(-1,2),F(1,2))
    slope=(hi**5-lo**5)/(hi-lo);intercept=lo**5-slope*lo
    add(a=1,d=-2*slope,rhs=2*intercept)
    for j in range(9):
        xi=lo+F(j,8)*(hi-lo)
        add(a=-1,d=10*xi**4,rhs=8*xi**5)
    return out

def family_parameters(family):
    if type(family) is not int or not 0<=family<3:raise ValueError('Bad family')
    # Symmetry identifies the two mixed cases.
    bits=[(0,0),(0,1),(1,1)][family]
    return tuple((F(bit,2),F(1-bit,2)) for bit in bits)

def witness_affine(family,n,z,j=0):
    x,y,w=z
    (by,cy),(bz,cz)=family_parameters(family)
    row=vector(a=-x*x,u=-2*x,v=-1,b=n*x,alpha=n*n,kappa=n)
    offset=n*j+n*(by*y+bz*w)-((F(y)+cy)**2+(F(w)+cz)**2)/2
    return row,offset

def witness_rows(family,n,z,j):
    row,offset=witness_affine(family,n,z,j);other=[-x for x in row];other[0]+=1
    return [(row,-offset),(other,offset)]
