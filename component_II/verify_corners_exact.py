"""Replay the three corner Farkas trees at the exact pinch points (Section 6, Appendix A.4) exactly;
root rows are regenerated from corner_model.root_model with the G-parameters fixed (eps = 0)."""
import json
from fractions import Fraction as F
from farkas_check import check_tree
import corner_model as C
CFG={"P1":(2,F(1,4),F(0),F(0),None,(F(0),F(0))),
     "P2":(2,F(1,4),F(0),F(0),None,(F(1,2),F(1,2))),
     "P3":(3,F(1,6),F(1,6),F(0),F(1,12),(F(7,12),F(1,2)))}
def main():
    out={}
    for name,(L,t0,e,c,mus,pinch) in CFG.items():
        M=C.root_model(L,t0,F(3,4),mus,pinch,eps=F(0))
        tree=json.load(open("certs/corner_%s.json"%name))["tree"]
        n=check_tree(M.rows,C.V,tree,C.witness_affine_factory(t0,e,c),({"d":F(1),"a":F(3,4)},F(0)))
        out[name]=dict(nodes=n[0],leaves=n[1],root_rows=len(M.rows))
    return out
if __name__=="__main__": print(main())
