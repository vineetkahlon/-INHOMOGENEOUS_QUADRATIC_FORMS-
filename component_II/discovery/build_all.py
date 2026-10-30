import json, sys, os
from multiprocessing import Pool
from specs import all_families
from cover import build
FAMS=None
def job(name):
    fam=all_families()[name]
    path="certs/"+name.replace("/","__")+".json"
    if os.path.exists(path): return name,"exists"
    try:
        tree,cnt=build(fam)
        json.dump(tree,open(path,"w"),separators=(",",":")); return name,cnt
    except Exception as ex: return name,"FAIL "+str(ex)
if __name__=="__main__":
    names=sorted(all_families())
    with Pool(8) as p:
        for name,res in p.imap_unordered(job,names): print(name,res,flush=True)
