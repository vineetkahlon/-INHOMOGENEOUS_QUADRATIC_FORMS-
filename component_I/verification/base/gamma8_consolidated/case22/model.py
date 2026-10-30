"""Rational outer models for BKR's (m,K)=(2,2), at constant eight.

Every physical point maps into the model; the converse is not asserted.
The verifier uses only this file and Python's standard library.
"""
from fractions import Fraction as F
from itertools import product,combinations
from math import isqrt

NAMES = ('d','a','A','t','p','q','r','s','u','b','v','beta1','beta2','beta3','nu')
DIM = len(NAMES)

def vector(**kw):
    out = [F(0)]*DIM
    for key, value in kw.items(): out[NAMES.index(key)] = F(value)
    return out

def d_bounds(L):
    # Rational lower bound on sqrt(125/(288*(L+1))).
    lo = F(isqrt(125*10**8//(288*(L+1))),10000)
    return lo,F(64,25*(L+3))

def upper_cuberoot(q,den=1000000):
    low=0;high=den
    while F(high,den)**3<q:high*=2
    while high-low>1:
        mid=(low+high)//2
        if F(mid,den)**3<q:low=mid
        else:high=mid
    return F(high,den)

def affine_fit(points):
    matrix=[[F(x) for x in [*p,1,p[0]*p[1]*p[2]]] for p in points]
    for j in range(4):
        pivot=next((i for i in range(j,4) if matrix[i][j]),None)
        if pivot is None:return None
        matrix[j],matrix[pivot]=matrix[pivot],matrix[j]
        value=matrix[j][j];matrix[j]=[x/value for x in matrix[j]]
        for i in range(4):
            if i==j:continue
            value=matrix[i][j];matrix[i]=[x-value*y for x,y in zip(matrix[i],matrix[j])]
    return tuple(row[-1] for row in matrix)

def product_facets(bounds):
    """All affine lower/upper bounds supported by the eight box corners.

    Multiaffinity makes the corner check sufficient on the entire box.
    """
    points=list(product(*bounds));lower=set();upper=set()
    for chosen in combinations(points,4):
        row=affine_fit(chosen)
        if row is None:continue
        errors=[sum(c*x for c,x in zip(row,[*p,1]))-p[0]*p[1]*p[2] for p in points]
        if max(errors)<=0:lower.add(row)
        if min(errors)>=0:upper.add(row)
    return sorted(lower),sorted(upper)

def coefficient_bounds(L,box=None):
    lo,hi=box or d_bounds(L)
    aa=(lo/3,hi/2)
    bb=(5*lo/12,min(11*hi/16,F(99,70)*hi*hi))
    tt=(25*lo/(16*(L+1)),min(121*hi/(64*L),F(18,5)*hi**3))
    if L>1:tt=(tt[0],min(tt[1],4*(1-25*lo/16)/(L*L-1)))
    return aa,bb,tt

def root_rows(L,version=3,box=None):
    if type(L) is not int or L not in (1,2,3,4): raise ValueError('Bad L')
    if version not in (1,2,3):raise ValueError('Bad model version')
    lo,hi=d_bounds(L)
    if box is not None:
        bl,bh=map(F,box)
        if not lo<=bl<bh<=hi:raise ValueError('Bad d interval')
        lo,hi=bl,bh
    out=[]
    def add(rhs=0,**kw): out.append((vector(**kw),F(rhs)))
    def bound(name,lo,hi): add(hi,**{name:1});add(-lo,**{name:-1})
    bound('d',lo,hi)
    add(d=F(1,3),a=-1);add(a=1,d=F(-1,2))
    add(d=F(1,3),a=F(1,4),A=-1)
    add(A=1,d=F(-1,2),a=F(-3,8))
    add(A=-1,a=F(3,4))
    add(A=F(3,4),t=-1)
    add(A=1,t=-1,v=-1)
    add(d=1,a=F(3,4),A=F(3,4),t=-(L+1))
    add(d=-1,a=F(-3,4),A=F(-3,4),t=L)
    add(1,d=1,a=F(3,4),A=F(3,4),t=F(L*L-1,4))
    add(1,d=1,a=1)
    # p=ah, q=ak, r=ah^2, s=ahk, u=ak^2.
    for lin,square in [('p','r'),('q','u')]:
        add(a=F(-1,2),**{lin:1});add(a=F(-1,2),**{lin:-1})
        add(a=F(-1,4),**{square:1})
        for xi in map(F,('-1/2','-1/4','0','1/4','1/2')):
            add(a=-xi*xi,**{lin:2*xi,square:-1})
    for i,j in [(1,1),(-1,-1),(1,-1),(-1,1)]:
        add(a=F(-1,4),p=F(i,2),q=F(j,2),s=-i*j)
    # b=A*lambda, v=A*lambda^2, with 0<=lambda<=1/2.
    add(b=-1);add(b=1,A=F(-1,2));add(v=1,A=F(-1,4))
    for xi in map(F,('0','1/8','1/4','3/8','1/2')):
        add(v=-1,b=2*xi,A=-xi*xi)
    for name in ('beta1','beta2','beta3'):bound(name,F(-1,2),F(1,2))
    bound('nu',F(0),F(1))
    # Finite necessary minimum constraints; phi(e1)=a already.
    for z in product(range(-2,3),repeat=3):
        if not 0<sum(x*x for x in z)<=5:continue
        if next(x for x in z if x)!=abs(next(x for x in z if x)):continue
        row=phi(z);row=[-x for x in row];row[1]+=1
        out.append((row,F(0)))
    if version>=2:
        # t <= (18/5)d^3; convex secant is an upper bound on d^3.
        add(t=1,d=-F(18,5)*(lo*lo+lo*hi+hi*hi),rhs=-F(18,5)*lo*hi*(lo+hi))
        # A <= sqrt(2)d^2 < (99/70)d^2; another convex secant.
        add(A=1,d=-F(99,70)*(lo+hi),rhs=-F(99,70)*lo*hi)
        # a <= d^(5/3); exact rational upper endpoints on its secant.
        a0,a1=upper_cuberoot(lo**5),upper_cuberoot(hi**5)
        slope=(a1-a0)/(hi-lo)
        add(a=1,d=-slope,rhs=a0-slope*lo)
    if version>=3:
        add(v=1,b=F(-1,2))
        lower,upper=product_facets(coefficient_bounds(L,(lo,hi)))
        slope=(hi**5-lo**5)/(hi-lo);intercept=lo**5-slope*lo
        for ca,cA,ct,c0 in lower:
            add(a=2*ca,A=2*cA,t=2*ct,d=-slope,rhs=intercept-2*c0)
        for xi in [lo,(3*lo+hi)/4,(lo+hi)/2,(lo+3*hi)/4,hi]:
            for ca,cA,ct,c0 in upper:
                add(a=-2*ca,A=-2*cA,t=-2*ct,d=5*xi**4,rhs=4*xi**5+2*c0)
    return out

def phi(z):
    x,y,w=z
    return vector(a=x*x,p=2*x*y,q=2*x*w,A=y*y,r=y*y,s=2*y*w,b=2*y*w,
                  t=w*w,v=w*w,u=w*w)

def witness_affine(z,j=0):
    c=[-v for v in phi(z)]
    for i,x in zip((11,12,13),z):c[i]=F(x)
    c[14]=F(1)
    return c,F(j)

def witness_rows(z,j):
    c,offset=witness_affine(z,j)
    other=[-x for x in c];other[0]+=1
    return [(c,-offset),(other,offset)]
