import sys, json
from fractions import Fraction as F
import corner22 as C
from farkas import prove, check_tree
name=sys.argv[1]
cfg={"P1":(2,F(1,4),F(0),F(0),None,(F(0),F(0))),"P2":(2,F(1,4),F(0),F(0),None,(F(1,2),F(1,2))),
     "P3":(3,F(1,6),F(1,6),F(0),F(1,12),(F(7,12),F(1,2)))}[name]
L,t0,e,c,mus,pinch=cfg
M=C.root_model(L,t0,F(3,4),mus,pinch,eps=F(0))
target=({"d":F(1),"a":F(3,4)},F(0))
tree,st=prove(M.rows,C.V,C.witnesses_factory(t0,e,c,target),target)
n=check_tree(M.rows,C.V,tree,C.witness_affine_factory(t0,e,c),target)
print(name,"rows",len(M.rows),st,"replay",n,flush=True)
json.dump({"tree":tree},open("cornerX_%s.json"%name,"w"))
