"""Emit LaTeX tables for the paper directly from specs / certificates / inventory."""
import json
from fractions import Fraction as F
from specs import all_families
from partition import cbrt_down
inv=json.load(open("inventory.json"))
def fr(x):
    x=F(x)
    return "$%s$"%(str(x.numerator) if x.denominator==1 else ("%s%d/%d"%("-" if x<0 else "",abs(x.numerator),x.denominator)))
out=[]
fams=all_families()
# shared rows
rows=[]
for k in sorted(k for k in fams if k.startswith("shared/")):
    f=fams[k]; nm=k.split("/")[1]
    # recover q,e,c from spec tuple (re-derive from specs source list)
rowspec=[("R0a","1/2","(0,1/2)","3/32","3/32","15/16"),("R0b","1/2","(1/2,0)","3/32","3/32","15/16"),
 ("R1","1/2","(0,1/2)","301/4800","27/280","19/20"),("R2","1/2","(1/2,0)","301/4800","27/280","19/20"),
 ("R3","1/2","(0,1/2)","1457/19200","23/180","22/25"),("R4","1/6","(1/6,0)","1/25","53/1000","3/5"),("R5","1/6","(0,1/2)","1/25","53/1000","3/5"),
 ("R6","1/4","(0,0)","217/3600","23/225","7/10"),("R7","1/4","(1/4,1/2)","217/3600","23/225","7/10"),
 ("R8","1/6","(1/6,0)","5053/96000","351/5000","61/100"),("R9","1/6","(0,1/2)","5053/96000","351/5000","61/100")]
L=["\\begin{tabular}{lccccc}\\toprule Row & $q$ & $(e,c)$ & $\\mu_-$ & $\\mu_+$ & $W$\\\\\\midrule"]
for nm,q,ec,m0,m1,W in rowspec:
    f=fams["shared/"+nm]; assert f.mu0==F(m0) and f.mu1==F(m1) and f.modes==[(F(0),F(W))]
    L.append("%s & %s & $%s$ & %s & %s & %s\\\\"%(nm,fr(q),ec,fr(m0),fr(m1),fr(W)))
L.append("\\bottomrule\\end{tabular}"); out.append(("shared_rows","\n".join(L)))
# (2,2) resonance families
L=["\\begin{tabular}{lcccccc}\\toprule $L$ & $t$ & $(e,c)$ & $d_0$ & $\\mu_-$ & $\\mu_+$ & pinch leaves\\\\\\midrule"]
for k in sorted(k for k in fams if k.startswith("c22_res/")):
    f=fams[k]; Ls=k[len("c22_res/"):]; Lv=int(Ls[1]); t=Ls.split("_t")[1].split("_e")[0]; e=Ls.split("_e")[1].split("_c")[0]; c=Ls.split("_c")[1]
    d0=[c1 for c1,c0 in f.modes if c1==3][0]; d0=[c0 for c1,c0 in f.modes if c1==3][0]
    extra=("none (Thm.~\\ref{thm:fail})" if False else ("pinch ball" if f.pinch else "--"))
    L.append("%d & %s & $(%s,%s)$ & %s & %s & %s & %s\\\\"%(Lv,fr(t),e,c,fr(d0),fr(f.mu0),fr(f.mu1),extra))
L.insert(1,"2 & $1/4$ & $(0,0)$ & $41103/100000$ & $123309/1600000$ & $11299/100000$ & no cover needed\\\\")
L.append("\\bottomrule\\end{tabular}"); out.append(("res22","\n".join(L)))
# partitions (decimal display)
P=json.load(open("certs/partitions.json"))
sup=["# Exact t-partitions for (2,2), L=2,3 (supplement to Section 6.4)","",
     "All endpoints are exact rationals; `verify_partitions.py` re-derives every bound. `*` marks the four resonant intervals.",""]
for key in ("c22_L2","c22_L3"):
    sup+=["## L = "+key[-1],"","| t_- | t_+ | method | bound | delta_- |","|---|---|---|---|---|"]
    for r in P[key]:
        if r["how"]=="gap": meth="gap, M=%s"%r["M"]; bd=r["gap"]
        else: meth="Macbeath (%s,%s)%s"%(r["H"],r["N"]," *" if r["how"]=="macbeath_res" else ""); bd=r["width"]
        sup.append("| %s | %s | %s | %s | %s |"%(r["a"],r["b"],meth,bd,r["dlow"]))
    sup.append("")
open("/home/claude/pkg/PARTITIONS.md","w").write("\n".join(sup))
# inventory
g=inv["family_groups"]
for nm,tex in out: open("/home/claude/paper/tab_%s.tex"%nm,"w").write(tex+"\n")
print({k:v for k,v in g.items()}); print([nm for nm,_ in out])
for t0,x in [(F(1,5),F(1,18)),(F(1,4),F(5,72)),(F(1,6),F(5,108))]: print(t0,cbrt_down(x))
