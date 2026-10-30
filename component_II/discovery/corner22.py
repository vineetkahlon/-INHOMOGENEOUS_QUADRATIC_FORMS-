"""Coupled F-level model for the (2,2) resonant corners.  F with x2=n, x4=y+c4, x5=z+c:
 f = n j + alpha n^2 + kappa n + e n z - A y^2 - 2 p y z + 2 r n y - 2 u y - v z^2 + 2 s n z - 2 w z - m - t0 (z+c)^2 - a/4
 target delta2 = d + 3a/4.  Region: (2,2) at fixed t=t0 (resonant), d+a+A <= tau, optional (a+A)/4 >= mustar."""
from farkas import *
from partition import cbrt_down
import sys, json, math
V=["d","a","A","P","p","v","r","s","u","w","m","Q2","X","aG","kG"]

def root_model(L,t0,tau,mustar=None,pinch=(F(0),F(0)),eps=F(1,128)):
    M=Model(V)
    dl=cbrt_down(F(5,18)*t0)     # aA = d^5/(2 t0) >= (d/3)(5d/12)
    dh=F(4,7)*tau                # 7d/4 <= d+a+A <= tau
    if mustar is not None: pass
    M.box("d",dl,dh)
    M.ge({"a":1,"d":F(-1,3)},0); M.le({"a":1,"d":F(-1,2)},0)
    # delta2/3 <= A <= delta2/2 ; A>=3a/4 ; 3A/4<=t0
    M.ge({"A":1,"d":F(-1,3),"a":F(-1,4)},0); M.le({"A":1,"d":F(-1,2),"a":F(-3,8)},0)
    M.ge({"A":1,"a":F(-3,4)},0); M.le({"A":F(3,4)},t0)
    # L t0 <= delta <= (L+1) t0 ; budget
    dco={"d":1,"a":F(3,4),"A":F(3,4)}
    M.ge(dco,L*t0); M.le(dco,min((L+1)*t0, 1-F(L*L-1,4)*t0))
    M.le({"d":1,"a":1,"A":1},tau)
    if mustar is not None: M.ge({"a":F(1,4),"A":F(1,4)},mustar)
    # P=aA: McCormick on box
    aL,aU=dl/3,dh/2; AL,AU=F(5,12)*dl,min(F(4,3)*t0,tau-F(4,3)*dl)
    M.box("a",aL,aU); M.box("A",AL,AU)
    M.ge({"P":1,"A":-aL,"a":-AL},-aL*AL); M.ge({"P":1,"A":-aU,"a":-AU},-aU*AU)
    M.le({"P":1,"A":-aU,"a":-AL},-aU*AL); M.le({"P":1,"A":-aL,"a":-AU},-aL*AU)
    # P = d^5/(2t0): secant above, tangents below
    k=1/(2*t0); sl=(dh**5-dl**5)/(dh-dl)
    M.le({"P":1,"d":-k*sl},k*(dl**5-sl*dl))
    for i in range(7):
        x0=dl+(dh-dl)*i/6; M.ge({"P":1,"d":-k*5*x0**4},k*(x0**5-5*x0**5))
    # moments: lambda in [0,1/2], q,C in [-1/2,1/2]
    M.ge({"p":1},0); M.le({"p":1,"A":F(-1,2)},0); M.le({"v":1,"p":F(-1,2)},0)
    for xi in [F(0),F(1,8),F(1,4),F(3,8),F(1,2)]: M.ge({"v":1,"p":-2*xi,"A":xi*xi},0)
    M.ge({"v":1,"A":-1},-t0)                       # t + A lambda^2 >= A
    for X in ("r","u"): M.le({X:1,"A":F(-1,2)},0); M.ge({X:1,"A":F(1,2)},0)
    for S,T in (("s","r"),("w","u")):
        M.ge({S:1,"p":F(1,2)},0); M.ge({S:1,"p":F(-1,2),T:F(-1,2),"A":F(1,4)},0)
        M.le({S:1,"p":F(1,2),T:F(-1,2),"A":F(-1,4)},0); M.le({S:1,"p":F(-1,2)},0)
    M.le({"m":1,"A":F(-1,4)},0)
    for xi in [F(-1,2),F(-1,4),F(0),F(1,4),F(1,2)]: M.ge({"m":1,"u":-2*xi,"A":xi*xi},0)
    for xi in [F(-2),F(-1),F(-1,2),F(0),F(1,2),F(1),F(2)]: M.ge({"m":1,"w":-2*xi,"v":xi*xi},0)
    # Q2=A q^2, X=A q C
    M.ge({"Q2":1},0); M.le({"Q2":1,"A":F(-1,4)},0)
    for xi in [F(-1,2),F(-1,4),F(0),F(1,4),F(1,2)]: M.ge({"Q2":1,"r":-2*xi,"A":xi*xi},0)
    for xi in [F(-2),F(-1),F(-1,2),F(0),F(1,2),F(1),F(2)]:
        M.ge({"v":1,"s":-2*xi,"Q2":xi*xi},0); M.ge({"m":1,"X":-2*xi,"Q2":xi*xi},0)
    M.ge({"X":1,"r":F(1,2),"u":F(1,2),"A":F(1,4)},0); M.ge({"X":1,"r":F(-1,2),"u":F(-1,2),"A":F(1,4)},0)
    M.le({"X":1,"r":F(1,2),"u":F(-1,2),"A":F(-1,4)},0); M.le({"X":1,"r":F(-1,2),"u":F(1,2),"A":F(-1,4)},0)
    M.box("aG",pinch[0]-eps,pinch[0]+eps); M.box("kG",pinch[1]-eps,pinch[1]+eps)
    return M

def witness_affine_factory(t0,e,c):
    def wa(lab):
        n,y,z,j=lab
        g={"aG":n*n,"Q2":-n*n,"kG":n,"X":2*n,"A":-y*y,"p":-2*y*z,"r":2*n*y,"u":-2*y,"v":-z*z,"s":2*n*z,"w":-2*z,"m":-1,"a":F(-1,4)}
        h=n*j+e*n*z+e*n*c-t0*(z+c)**2
        return {k:F(v) for k,v in g.items() if v!=0},F(h)
    return wa

def witnesses_factory(t0,e,c,target):
    wa=witness_affine_factory(t0,e,c); tco,tc=target
    def W(th):
        T=sum(float(tco[k])*th[k] for k in tco)+float(tc); out=[]
        for n in [1,-1,2,-2,3,-3,4,-4]:
            for y in range(-3,4):
                for z in range(-3,4):
                    g,h=wa((n,y,z,0))
                    rest=sum(float(g[k])*th[k] for k in g)+float(h)
                    j0=math.floor((T/2-rest)/n)
                    for j in (j0,j0+1):
                        g2,h2=wa((n,y,z,j)); out.append(((n,y,z,j),g2,h2))
        return out
    return W

if __name__=="__main__":
    name=sys.argv[1]
    cfg={"L2_t14_00":(2,F(1,4),F(0),F(0),F(153,200),None,(F(0),F(0))),
         "L2_t14_00_half":(2,F(1,4),F(0),F(0),F(153,200),None,(F(1,2),F(1,2))),
         "L3_t16_160":(3,F(1,6),F(1,6),F(0),F(153,200),F(3,40),(F(7,12),F(1,2)))}[name]
    L,t0,e,c,tau,mus,pinch=cfg
    M=root_model(L,t0,tau,mus,pinch)
    target=({"d":F(1),"a":F(3,4)},F(0))
    tree,st=prove(M.rows,V,witnesses_factory(t0,e,c,target),target)
    print(name,"root rows",len(M.rows),st,flush=True)
    n=check_tree(M.rows,V,tree,witness_affine_factory(t0,e,c),target); print("exact replay nodes,leaves",n)
    json.dump({"cfg":[str(x) for x in (L,t0,e,c,tau,mus,pinch[0],pinch[1])],"tree":tree},open("corner_%s.json"%name,"w"))
