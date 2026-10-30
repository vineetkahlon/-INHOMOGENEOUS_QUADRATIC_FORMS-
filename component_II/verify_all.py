"""Master verifier for the new certificate component (exact rational arithmetic only; numpy/scipy not imported).
Run:  python3 verify_all.py > inventory.json"""
import json, os, sys, hashlib, time
from fractions import Fraction as F
def sha(p): return hashlib.sha256(open(p,"rb").read()).hexdigest()[:16]
def main():
    t0=time.time(); inv={}
    from specs import all_families
    from cover import check
    fams=all_families()
    files={f[:-5] for f in os.listdir("certs") if f.endswith(".json") and "__" in f}
    want={k.replace("/","__") for k in fams}
    if files!=want: raise SystemExit("inventory mismatch: missing %s extra %s"%(want-files,files-want))
    groups={}
    for k,fam in sorted(fams.items()):
        p="certs/"+k.replace("/","__")+".json"
        st=check(fam,json.load(open(p)))
        g=k.split("/")[0]; G=groups.setdefault(g,{"families":0,"nodes":0,"leaves":0})
        G["families"]+=1; G["nodes"]+=st["nodes"]; G["leaves"]+=st["leaves"]; G["pinch_leaves"]=G.get("pinch_leaves",0)+st.get("pinch_leaves",0)
        inv[k]=dict(st,sha=sha(p))
    import verify_partitions, verify_gaps, verify_31, verify_corners_exact
    part=verify_partitions.main(); ngap,twopiece=verify_gaps.main(); t31,t31res=verify_31.main(); cor=verify_corners_exact.main()
    import verify_scalars, verify_pinch, verify_sixvar; nsc=verify_scalars.main(); pin=verify_pinch.main(); six=verify_sixvar.main()
    out=dict(family_groups=groups,families=inv,partitions=part,gap_rows=ngap,two_piece_bound=twopiece,
             table_31=t31,table_31_n2_residues=t31res,corners=cor,pinch_ladders=pin,sixvar_supplement_checks=six,scalar_checks=nsc,seconds=round(time.time()-t0,1),status="ALL CHECKS PASSED")
    return out
if __name__=="__main__":
    print(json.dumps(main(),indent=1))
