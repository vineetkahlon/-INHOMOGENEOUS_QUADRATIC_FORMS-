"""Coupled binary F model for the t=1/6 Macbeath exceptions."""
from fractions import Fraction as F
NAMES=('d','a','A','p','v','r','s','u','w','mu','alpha','kappa')
DIM=len(NAMES)

def vector(**kw):
    row=[F(0)]*DIM
    for name,value in kw.items():row[NAMES.index(name)]=F(value)
    return row

def upper_root(q,power,den=1000000):
    lo=0;hi=den
    while hi-lo>1:
        mid=(hi+lo)//2
        if F(mid,den)**power<q:lo=mid
        else:hi=mid
    return F(hi,den)

def d_bounds():return F(3968,10000),upper_root(F(2,27)**3,10)

def root_rows(shape_box=None):
    lo,hi=d_bounds();out=[]
    def add(rhs=0,**kw):out.append((vector(**kw),F(rhs)))
    def bound(name,ll,hh):add(hh,**{name:1});add(-ll,**{name:-1})
    bound('d',lo,hi)
    add(d=F(1,2),a=-1);add(a=1,d=F(-2,3));add(1,a=1,d=1)
    add(d=F(3,8),A=-1);add(A=1,d=F(-1,2));add(a=F(3,4),A=-1)
    add(F(2,9),A=1);add(upper_root(F(2,27),2),a=1)
    add(F(-1,2),d=-1,A=F(-3,4));add(F(2,3),d=1,A=F(3,4))
    bound('alpha',F(-1,2),F(1,2));bound('kappa',F(-1,2),F(1,2))
    lbox,qbox,cbox=[tuple(map(F,pair)) for pair in (shape_box or [(0,F(1,2)),(F(-1,2),F(1,2)),(F(-1,2),F(1,2))])]
    if not (0<=lbox[0]<lbox[1]<=F(1,2) and F(-1,2)<=qbox[0]<qbox[1]<=F(1,2) and F(-1,2)<=cbox[0]<cbox[1]<=F(1,2)):
        raise ValueError('Bad shape box')
    for name,(ll,hh) in [('p',lbox),('r',qbox),('u',cbox)]:
        add(A=ll,**{name:-1});add(A=-hh,**{name:1})
    for lin,square,(ll,hh) in [('p','v',lbox),('u','mu',cbox)]:
        add(A=ll*hh,**{square:1,lin:-(ll+hh)})
        for j in range(5):
            xi=ll+F(j,4)*(hh-ll)
            add(A=-xi*xi,**{lin:2*xi,square:-1})
    add(F(1,6),A=1,v=-1)
    def product_rows(other,prod,box):
        l0,l1=lbox;q0,q1=box
        add(A=-l0*q0,p=q0,**{other:l0,prod:-1})
        add(A=-l1*q1,p=q1,**{other:l1,prod:-1})
        add(A=l1*q0,p=-q0,**{other:-l1,prod:1})
        add(A=l0*q1,p=-q1,**{other:-l0,prod:1})
    product_rows('r','s',qbox);product_rows('u','w',cbox)
    for xi in map(F,('-2','-1','-1/2','0','1/2','1','2')):add(mu=-1,w=2*xi,v=-xi*xi)
    # a <= d^(5/3), bounded by an outward-rounded secant.
    aa0,aa1=upper_root(lo**5,3),upper_root(hi**5,3)
    slope=(aa1-aa0)/(hi-lo);add(a=1,d=-slope,rhs=aa0-slope*lo)
    # a A = 3 d^5; use exact bilinear box bounds and power secants/tangents.
    al,ah=lo/2,min(2*hi/3,upper_root(F(2,27),2),upper_root(hi**5,3))
    Al,Ah=3*lo/8,min(hi/2,F(2,9))
    secant=(hi**5-lo**5)/(hi-lo);intercept=lo**5-secant*lo
    for ca,cA,c0 in [(Al,al,-al*Al),(Ah,ah,-ah*Ah)]:
        add(a=ca,A=cA,d=-3*secant,rhs=3*intercept-c0)
    for xi in [lo,(3*lo+hi)/4,(lo+hi)/2,(lo+3*hi)/4,hi]:
        for ca,cA,c0 in [(Al,ah,-ah*Al),(Ah,al,-al*Ah)]:
            add(a=-ca,A=-cA,d=15*xi**4,rhs=12*xi**5+c0)
    return out

def witness_affine(family,n,z,j=0):
    if family not in (0,1):raise ValueError('Bad family')
    y,w=z;c=F(family,2);e=F(1-family,6)
    row=vector(a=F(-1,4),A=-y*y,p=-2*y*w,r=2*n*y,u=-2*y,
               v=-w*w,s=2*n*w,w=-2*w,mu=-1,alpha=n*n,kappa=n)
    offset=n*j+e*n*w-(F(w)+c)**2/6
    return row,offset

def witness_rows(family,n,z,j):
    row,offset=witness_affine(family,n,z,j);other=[-x for x in row];other[0]+=1
    return [(row,-offset),(other,offset)]
