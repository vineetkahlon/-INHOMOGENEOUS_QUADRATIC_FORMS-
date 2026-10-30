"""Write MANIFEST.json and SHA256SUMS covering every file of the release (standard library only)."""
import hashlib, json, os, sys
ROOT=os.path.dirname(os.path.abspath(__file__))
SKIP={"MANIFEST.json","SHA256SUMS"}
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda: f.read(1<<20), b""): h.update(b)
    return h.hexdigest()
def main():
    if not os.path.exists(os.path.join(ROOT,"component_I","verify_end_to_end.py")):
        sys.exit("component_I/verify_end_to_end.py not found: copy the earlier archive into component_I/ first")
    if os.path.exists(os.path.join(ROOT,"component_I","PUT_ARCHIVE_HERE.md")):
        sys.exit("remove component_I/PUT_ARCHIVE_HERE.md after copying the archive")
    usage=json.load(open(os.path.join(ROOT,"component_I_usage.json")))
    missing=[u["case"] for u in usage["used"] if not u["path"]]
    if missing: sys.exit("fill the 'path' fields in component_I_usage.json for: %s"%missing)
    for u in usage["used"]:
        for p,h in zip(u["path"].split(";"),u["sha256"]):
            q=os.path.join(ROOT,"component_I",p.strip())
            if not os.path.exists(q): sys.exit("used tree not found: %s"%p)
            if sha(q)!=h: sys.exit("used tree differs from the pinned 2026-10-29 hash: %s"%p)
    files={}
    for d,dirs,fs in os.walk(ROOT):
        dirs[:]=[x for x in dirs if x!='.git']
        for f in fs:
            rel=os.path.relpath(os.path.join(d,f),ROOT)
            if rel in SKIP or "__pycache__" in rel or rel.endswith("_report.json") or rel.split(os.sep)[0]==".git": continue
            files[rel]=sha(os.path.join(ROOT,rel))
    json.dump({"files":dict(sorted(files.items())),"component_I_usage":usage},open(os.path.join(ROOT,"MANIFEST.json"),"w"),indent=1)
    with open(os.path.join(ROOT,"SHA256SUMS"),"w") as f:
        for k,v in sorted(files.items()): f.write("%s  %s\n"%(v,k))
    print("manifest written:",len(files),"files")
if __name__=="__main__": main()
