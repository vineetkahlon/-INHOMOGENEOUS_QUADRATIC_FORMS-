"""Exact Farkas-tree replay (no LP, no floats)."""
from fractions import Fraction as F
def check_leaf(rows,V,cert):
    if any(x<0 for x in cert.values()): return False
    tot={v:F(0) for v in V}; rhs=F(0)
    for i,wt in cert.items():
        co,r=rows[i]
        for k,v in co.items(): tot[k]+=wt*v
        rhs+=wt*r
    return all(x==0 for x in tot.values()) and rhs<0

def check_tree(root,V,tree,witness_affine,target):
    """independent replay: witness_affine(label)->(g,h); exact."""
    tco,tc=target; n=[0,0]
    def rec(rows,node):
        n[0]+=1
        if "leaf" in node:
            n[1]+=1; cert={int(k):F(v) for k,v in node["leaf"].items()}
            assert all(0<=k<len(rows) for k in cert); 
            if not check_leaf(rows,V,cert): raise ValueError("bad leaf")
            return
        g,h=witness_affine(node["w"]); g={k:F(v) for k,v in g.items()}; h=F(h)
        r1=(g,-h); r2=({k:tco.get(k,F(0))-g.get(k,F(0)) for k in set(g)|set(tco)},h-tc)
        if len(node["c"])!=2: raise ValueError("missing child")
        rec(rows+[r1],node["c"][0]); rec(rows+[r2],node["c"][1])
    rec(list(root),tree); return n
