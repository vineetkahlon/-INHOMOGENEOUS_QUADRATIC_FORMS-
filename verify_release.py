"""Release-level check: (1) every file matches MANIFEST.json and no unlisted file exists;
(2) component II replays (verify_all.py); (3) component I replays with its own verify_end_to_end.py.
Standard library only.  Run from anywhere:  python3 verify_release.py"""
import hashlib, json, os, subprocess, sys
ROOT=os.path.dirname(os.path.abspath(__file__))
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda: f.read(1<<20), b""): h.update(b)
    return h.hexdigest()
def run(cmd,cwd):
    print("  $",cmd,"   (in %s)"%os.path.relpath(cwd,ROOT)); r=subprocess.run(cmd,shell=True,cwd=cwd,capture_output=True,text=True)
    return r
def main():
    man=json.load(open(os.path.join(ROOT,"MANIFEST.json")))["files"]
    present={}
    for d,dirs,fs in os.walk(ROOT):
        dirs[:]=[x for x in dirs if x!='.git']
        for f in fs:
            rel=os.path.relpath(os.path.join(d,f),ROOT)
            if rel in ("MANIFEST.json","SHA256SUMS") or "__pycache__" in rel or rel.endswith("_report.json") or rel.split(os.sep)[0]==".git": continue
            present[rel]=sha(os.path.join(ROOT,rel))
    bad=[k for k in man if present.get(k)!=man[k]]; extra=[k for k in present if k not in man]
    if bad or extra: sys.exit("manifest mismatch: changed/missing %s, unlisted %s"%(bad,extra))
    usage=json.load(open(os.path.join(ROOT,"component_I_usage.json")))
    for u in usage["used"]:
        for p,h in zip(u["path"].split(";"),u["sha256"]):
            if sha(os.path.join(ROOT,"component_I",p))!=h: sys.exit("used component-I tree differs from pinned hash: "+p)
    print("[1/3] manifest: %d files match; 25 used component-I trees match pinned hashes"%len(man))
    r=run("python3 -B verify_all.py > ../component_II_report.json",os.path.join(ROOT,"component_II"))
    rep=json.load(open(os.path.join(ROOT,"component_II_report.json")))
    if r.returncode or rep.get("status")!="ALL CHECKS PASSED": sys.exit("component II failed:\n"+r.stderr[-2000:])
    print("[2/3] component II:",rep["status"],"(%s s)"%rep["seconds"])
    r=run("python3 -S -B verify_end_to_end.py > ../component_I_report.json",os.path.join(ROOT,"component_I"))
    if r.returncode: sys.exit("component I verifier failed:\n"+r.stderr[-2000:])
    rep1=json.load(open(os.path.join(ROOT,"component_I_report.json")))
    ref=json.load(open(os.path.join(ROOT,"component_I_reference_replay_2026-10-29.json")))
    if rep1.get("status")!="PASS": sys.exit("component I report status is not PASS")
    for k in ("farkas_trees","nodes","terminal_contradictions"):
        if rep1["base"][k]!=ref["base"][k]: sys.exit("component I %s differs from reference: %s"%(k,rep1["base"][k]))
    if rep1["certificate_data_sha256"]!=ref["certificate_data_sha256"]: sys.exit("component I certificate hashes differ from reference")
    print("[3/3] component I: PASS, %d trees / %d nodes, all certificate hashes match the 2026-10-29 reference"%(rep1["base"]["farkas_trees"],rep1["base"]["nodes"]))
    print("RELEASE VERIFIED")
if __name__=="__main__": main()
