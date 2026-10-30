"""Exact rational outer models for Section 8, (m,K)=(2,1), at constant eight."""
from fractions import Fraction as F
from itertools import product,combinations
NAMES=('d','a','A','t','p','q','r','s','u','b','v','beta1','beta2','beta3','nu')
DIM=len(NAMES)

def vector(**kw):
    out=[F(0)]*DIM
    for key,value in kw.items():out[NAMES.index(key)]=F(value)
    return out

def d_bounds():
    return F(4419,10000), F(8,15)

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

SQRT2_UPPER=F(1414214,1000000)
def coefficient_bounds():
    lo,hi=d_bounds()
    aa=(lo/3,min(hi/2,1-lo,upper_cuberoot(hi**5)))
    bb=(5*lo/8,min(11*hi/8,SQRT2_UPPER*hi*hi,1-5*lo/4))
    tt=(15*lo/32,F(12,5)*hi**3)
    return aa,bb,tt

def root_rows(sign):
    if type(sign) is not int or sign not in (-1,1):raise ValueError('Bad shear sign')
    lo,hi=d_bounds()
    shape_box=[(F(0),F(1,2)),(F(0),F(1,2)) if sign==1 else (F(-1,2),F(0)),(F(0),F(1,2))]
    out=[]
    def add(rhs=0,**kw):out.append((vector(**kw),F(rhs)))
    def bound(name,low,high):add(high,**{name:1});add(-low,**{name:-1})
    bound('d',lo,hi)
    add(d=F(1,3),a=-1);add(a=1,d=F(-1,2))
    add(d=F(1,2),a=F(3,8),A=-1);add(A=1,d=-1,a=F(-3,4))
    add(A=-1,a=F(3,4));add(A=F(3,4),t=-1);add(A=1,t=-1,v=-1)
    add(1,d=1,a=1)
    add(1,A=1,d=1,a=F(3,4))
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
    a0,a1=upper_cuberoot(lo**5),upper_cuberoot(hi**5)
    slope=(a1-a0)/(hi-lo);add(a=1,d=-slope,rhs=a0-slope*lo)
    slope=(hi**5-lo**5)/(hi-lo);intercept=lo**5-slope*lo
    add(t=1,d=-F(12,5)*(lo*lo+lo*hi+hi*hi),rhs=-F(12,5)*lo*hi*(lo+hi))
    add(A=1,d=-SQRT2_UPPER*(lo+hi),rhs=-SQRT2_UPPER*lo*hi)
    lower,upper=product_facets(coefficient_bounds())
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
