"""All fixed families of the verification package (each justified in the paper; section noted)."""
from fractions import Fraction as F
from families import single, diag_block, gram_block
from cover import Family
from partition import cbrt_down
h=F(1,2)
def large_a_classes():
    """Sec. 4 (large ternary minimum): 72 forced classes.  Returns list of (name, Gram, b, c, Delta)."""
    out=[]
    pairs=[(F(0),h),(h,F(0))]
    for r,Dl in [(1,F(151,200)),(2,F(87,100)),(3,F(15,16))]:
        M=[[h,0,0],[0,h,0],[0,0,F(r,2)]]
        for (b1,c1) in pairs:
            for (b2,c2) in pairs:
                for b3 in (F(0),h):
                    for j in range(r):
                        c3=((b3-F(r,2)+j)/r)%1
                        out.append(("phi%d"%r,M,[b1,b2,b3],[c1,c2,c3],Dl))
    M=[[h,0,0],[0,F(1),h],[0,h,F(1)]]
    for (b1,c1) in pairs:
        for b2 in (F(0),h):
            for b3 in (F(0),h):
                sols=[(c2,c3) for c2 in [F(k,6) for k in range(6)] for c3 in [F(k,6) for k in range(6)]
                      if (2*c2+c3-b2)%1==0 and (c2+2*c3-b3)%1==0]
                assert len(sols)==3, sols
                for c2,c3 in sols: out.append(("phiH",M,[b1,b2,b3],[c1,c2,c3],F(15,16)))
    assert len(out)==72
    return out

def all_families():
    fams={}
    for i,(nm,M,b,c,Dl) in enumerate(large_a_classes()):
        fams["large_a/%02d_%s"%(i,nm)]=Family(gram_block(M,b,c,[-1,0,1]),4,0,0,[(F(0),Dl)])
    # Sec. (2,2) L=1 and (1,2) L=1: t=1/2, A=1/4, lambda=0; mu=a/4
    for (e4,c4) in [(F(0),F(0)),(F(1,4),h)]:
        for (e5,c5) in [(h,F(0)),(F(0),h)]:
            K=diag_block([(F(1,4),e4,c4),(h,e5,c5)],range(-3,4))
            fams["c22_L1/e4_%s_e5_%s"%(e4,e5)]=Family(K,4,F(1,24),F(3,44),[(F(0),F(2,3))])
            fams["c12_L1/e4_%s_e5_%s"%(e4,e5)]=Family(K,4,F(29,400),F(1,12),[(F(0),F(7,12))])
    # Sec. (3,1): a=1/6 resonance.  Seven of the eight classes are settled by n=2 residues (verify_31.py);
    # only the all-half class b=e=f=0, c=c_y=c_z=1/2 needs a cover; target 3/5
    K=diag_block([(F(1,6),F(0),h),(h,F(0),h),(h,F(0),h)],range(-3,4))
    fams["c31/b_0_e_0_f_0"]=Family(K,2,0,0,[(F(0),F(3,5))])
    # Shared one-variable resonant rows for (3,2),(2,4),(2,3)
    rows=[("R0a",h,(F(0),h),F(3,32),F(3,32),F(15,16)),("R0b",h,(h,F(0)),F(3,32),F(3,32),F(15,16)),("R1",h,(F(0),h),F(301,4800),F(27,280),F(19,20)),("R2",h,(h,F(0)),F(301,4800),F(27,280),F(19,20)),
          ("R3",h,(F(0),h),F(1457,19200),F(23,180),F(22,25)),
          ("R4",F(1,6),(F(1,6),F(0)),F(1,25),F(53,1000),F(3,5)),("R5",F(1,6),(F(0),h),F(1,25),F(53,1000),F(3,5)),
          ("R6",F(1,4),(F(0),F(0)),F(217,3600),F(23,225),F(7,10)),("R7",F(1,4),(F(1,4),h),F(217,3600),F(23,225),F(7,10)),
          ("R8",F(1,6),(F(1,6),F(0)),F(5053,96000),F(351,5000),F(61,100)),("R9",F(1,6),(F(0),h),F(5053,96000),F(351,5000),F(61,100))]
    for nm,q,(e,c),m0,m1,W in rows:
        fams["shared/"+nm]=Family(single(q,e,c,10),6,m0,m1,[(F(0),W)])
    # (2,2) L=2,3 resonances (G-level): modes delta>=d_lo+3mu, delta>121mu/19, delta>L t
    # family (2,1/4,0,0) has no cover: Theorem 6.x (pinch balls of radius 1/2 about (0,0),(1/2,1/2) cover the torus)
    def p22(L,t0,dlo):
        dmax=min((L+1)*t0,1-F(L*L-1,4)*t0); return F(3,16)*dlo,(dmax-dlo)/3,[(F(3),dlo),(F(121,19),F(0)),(F(0),L*t0)]
    for L,t0,dlo,forms in [(2,F(1,5),max(F(3803,10000),cbrt_down(F(1,18))),[(F(0),F(0)),(F(1,10),F(1,4)),(F(1,5),h)]),
                           (2,F(1,4),max(F(3803,10000),cbrt_down(F(5,72))),[(F(1,4),h)]),          # (0,0): Theorem, no cover
                           (3,F(1,6),max(F(3294,10000),cbrt_down(F(5,108))),[(F(1,6),F(0)),(F(0),h)])]:
        m0,m1,md=p22(L,t0,dlo)
        for e,c in forms:
            pin=[(F(7,12),h,h)] if (L,t0,e)==(3,F(1,6),F(1,6)) else []
            fams["c22_res/L%d_t%s_e%s_c%s"%(L,t0,e,c)]=Family(single(t0,e,c,12),12,m0,m1,md,pinch=pin)
    return fams
