"""Exact rational outer models for (m,K)=(1,2), at constant eight."""
from fractions import Fraction as F
from itertools import product,combinations
NAMES=('d','a','A','t','p','q','r','s','u','b','v','beta1','beta2','beta3','nu')
DIM=len(NAMES)

def vector(**kw):
    out=[F(0)]*DIM
    for key,value in kw.items():out[NAMES.index(key)]=F(value)
    return out

def lower_sqrt(q,den=10000):
    from math import isqrt
    q=F(q)
    return F(isqrt(q.numerator*den*den//q.denominator),den)

def upper_budget_root(L,den=1000000):
    lo=0;hi=den
    def poly(x):return 27*(L*L+3)**3*x**5-8192*(1-x)**3
    while hi-lo>1:
        mid=(lo+hi)//2
        if poly(F(mid,den))<0:lo=mid
        else:hi=mid
    return F(hi,den)

def d_bounds(L):
    if type(L) is not int or L not in (1,2,3,4):raise ValueError('Bad L')
    lo=lower_sqrt(max(F(1,8),F(123,256*(L+1))))
    hi=F(2,3) if L==1 else min(F(2,3),F(128,41*(L+3)),F(128,155+9*L*L),upper_budget_root(L))
    return lo,hi

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

SQRT43_UPPER=F(577351,500000)

def coefficient_bounds(L,box=None):
    lo,hi=box or d_bounds(L)
    aa=(lo/2,min(2*hi/3,1-lo,F(413,1000),upper_cuberoot(hi**5)))
    bb=(3*lo/8,min(hi/2,SQRT43_UPPER*hi*hi))
    tt=(max(9*lo/32,41*lo/(32*(L+1)),F(3,2)*lo**3),
        min(F(8,3)*hi**3,11*hi/(8*L),4*(1-41*lo/32)/(L*L-1) if L>1 else 11*hi/8))
    return aa,bb,tt

def root_rows(L,version=2,box=None,shape_box=None,region=None):
    if type(L) is not int or L not in (1,2,3,4) or version!=2:raise ValueError('Bad model')
    lo,hi=d_bounds(L)
    if box is not None:
        bl,bh=map(F,box)
        if not lo<=bl<bh<=hi:raise ValueError('Bad d interval')
        lo,hi=bl,bh
    out=[]
    def add(rhs=0,**kw):out.append((vector(**kw),F(rhs)))
    def bound(name,low,high):add(high,**{name:1});add(-low,**{name:-1})
    bound('d',lo,hi)
    add(d=F(1,2),a=-1);add(a=1,d=F(-2,3));add(F(413,1000),a=1)
    add(d=F(3,8),A=-1);add(A=1,d=F(-1,2))
    add(A=-1,a=F(3,4));add(A=F(3,4),t=-1);add(A=1,t=-1,v=-1)
    add(d=1,A=F(3,4),t=-(L+1));add(d=-1,A=F(-3,4),t=L)
    add(1,d=1,a=1)
    add(1,d=1,A=F(3,4),t=F(L*L-1,4))
    if L>=2:add(A=3*(L-1),d=-4);add(a=9*(L-1),d=-16)
    for lin,square in [('p','r'),('q','u')]:
        add(a=F(-1,2),**{lin:1});add(a=F(-1,2),**{lin:-1})
        add(a=F(-1,4),**{square:1})
        for xi in map(F,('-1/2','-1/4','0','1/4','1/2')):
            add(a=-xi*xi,**{lin:2*xi,square:-1})
    for i,j in [(1,1),(-1,-1),(1,-1),(-1,1)]:
        add(a=F(-1,4),p=F(i,2),q=F(j,2),s=-i*j)
    add(b=-1);add(b=1,A=F(-1,2));add(v=1,b=F(-1,2))
    for xi in map(F,('0','1/8','1/4','3/8','1/2')):add(v=-1,b=2*xi,A=-xi*xi)
    for name in ('beta1','beta2','beta3'):bound(name,F(-1,2),F(1,2))
    bound('nu',F(0),F(1))
    for z in product(range(-2,3),repeat=3):
        if not 0<sum(x*x for x in z)<=5:continue
        if next(x for x in z if x)<0:continue
        row=[-x for x in phi(z)];row[1]+=1;out.append((row,F(0)))
    add(t=1,d=-F(8,3)*(lo*lo+lo*hi+hi*hi),rhs=-F(8,3)*lo*hi*(lo+hi))
    for xi in [lo,(lo+hi)/2,hi]:add(t=-1,d=F(9,2)*xi*xi,rhs=3*xi**3)
    add(A=1,d=-SQRT43_UPPER*(lo+hi),rhs=-SQRT43_UPPER*lo*hi)
    a0,a1=upper_cuberoot(lo**5),upper_cuberoot(hi**5)
    slope=(a1-a0)/(hi-lo);add(a=1,d=-slope,rhs=a0-slope*lo)
    lower,upper=product_facets(coefficient_bounds(L,(lo,hi)))
    slope=(hi**5-lo**5)/(hi-lo);intercept=lo**5-slope*lo
    for ca,cA,ct,c0 in lower:add(a=2*ca,A=2*cA,t=2*ct,d=-slope,rhs=intercept-2*c0)
    for xi in [lo,(3*lo+hi)/4,(lo+hi)/2,(lo+3*hi)/4,hi]:
        for ca,cA,ct,c0 in upper:add(a=-2*ca,A=-2*cA,t=-2*ct,d=5*xi**4,rhs=4*xi**5+2*c0)
    if shape_box is not None:
        hb,kb,lb=[tuple(map(F,pair)) for pair in shape_box]
        if not (F(-1,2)<=hb[0]<hb[1]<=F(1,2) and F(-1,2)<=kb[0]<kb[1]<=F(1,2) and 0<=lb[0]<lb[1]<=F(1,2)):
            raise ValueError('Bad shape box')
        for lin,square,(ll,uu) in [('p','r',hb),('q','u',kb)]:
            add(a=ll,**{lin:-1});add(a=-uu,**{lin:1})
            add(a=ll*uu,**{square:1,lin:-(ll+uu)})
            for xi in [ll,(ll+uu)/2,uu]:add(a=-xi*xi,**{lin:2*xi,square:-1})
        h0,h1=hb;k0,k1=kb
        add(a=-h0*k0,p=k0,q=h0,s=-1);add(a=-h1*k1,p=k1,q=h1,s=-1)
        add(a=h1*k0,p=-k0,q=-h1,s=1);add(a=h0*k1,p=-k1,q=-h0,s=1)
        l0,l1=lb
        add(A=l0,b=-1);add(A=-l1,b=1);add(A=l0*l1,b=-(l0+l1),v=1)
        for xi in [l0,(l0+l1)/2,l1]:add(A=-xi*xi,b=2*xi,v=-1)
    if region is not None:
        if L==1 and region=='L1_outer':add(1,d=1,A=F(3,4),t=1)
        elif L==3 and region=='L3_low':add(2,d=1,A=F(3,4),t=9)
        elif L==3 and region=='L3_high':add(-1,d=1,A=F(3,4),t=-9)
        elif L==1 and region=='L1_resonance':bound('t',F(1,2),F(1,2))
        elif L==1 and region=='L1_A_low':
            bound('t',F(1,2),F(1,2));add(F(3,2),d=1,A=4)
        elif L==1 and region=='L1_A_high':
            bound('t',F(1,2),F(1,2));add(F(-1,2),d=1,A=-4)
        elif L==3 and region=='L3_resonance':bound('t',F(1,6),F(1,6))
        else:raise ValueError('Bad region')
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
