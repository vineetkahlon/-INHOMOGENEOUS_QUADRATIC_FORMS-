"""Master verifier for the consolidated constant-eight proof.

Standard library only. Mathematical domain inclusions and cited theorems
are not formalized by this program. Run without -O.
"""
import hashlib,json,subprocess,sys,time
from fractions import Fraction as F
from itertools import product
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def require(condition,message):
    if not condition:raise ValueError(message)

def scalar_audit():
    require(8**4>512*6,'Small-a Jackson threshold')
    for m in range(1,5):
        B=min(F(m+1,m+2),F(4,m+3))
        require(F(3*(m+3)**2,32*(m+1)*B*B)>1,'K lower bound')
    for m,K in ((4,2),(3,3),(2,5)):
        require(F((m+3)**4*(K+3)**3,(m+1)*(K+1))>4**7,'Jackson cutoff')
    require(F(8,3)<3,'K<=2 for m=1')
    for u in product((1,3,5,7),repeat=5):
        require((4*u[0]*u[1]-sum(x*x for x in u[1:]))%8==0,'Sharpness congruence')
    x=(F(3,2),)+(F(1,2),)*4
    require(x[0]*x[1]-sum(v*v for v in x[1:])/4==F(1,2),'Sharpness attainment')
    # Determinant of the 2x2 block is -1/4 and the remaining diagonal is -1/4.
    determinant=F(-1,4)*F(-1,4)**3
    require(determinant==F(1,256),'Sharpness determinant')
    require(F(1,2)**5/determinant==8,'Sharpness normalization')
    return {'classification':'passed','sharpness':'passed','odd_residue_tuples':4**5}

def expected_inventory():
    inventories={'large_a':{f'family_{i:02d}.json' for i in range(72)},
      'case21':{'low_plus.json','low_minus.json',*(f'diagonal_{i}.json' for i in range(3))},
      'case22':{*(f'L{i}.json' for i in (2,3,4)),*(f'family_{i}.json' for i in range(4))},
      'case11':{f'L{L}_{s}.json' for L in (1,2) for s in ('pp','pn')}}
    case12={f'L{L}_{s}.json' for L in (2,4) for s in ('pp','pn')}
    case12|={f'{r}_{s}.json' for r in ('L1_outer','L1_A_low','L1_A_high','L3_high') for s in ('pp','pn')}
    case12|={f'L3_low_{s}{t}.json' for s in ('pp','pn') for t in ('a','b')}
    case12|={f'four_{i}.json' for i in range(4)}|{f'resonance_{i}.json' for i in range(2)}
    inventories['case12']=case12
    return inventories

def run_json(relative):
    done=subprocess.run([sys.executable,'-S','-B',str(ROOT/relative)],cwd=ROOT,
                         text=True,capture_output=True)
    if done.returncode:
        raise RuntimeError(f'{relative} failed:\n{done.stdout}\n{done.stderr}')
    return json.loads(done.stdout)

def run():
    require(__debug__,'Run without optimization (-O): some component checks use assertions')
    start=time.monotonic();inventories=expected_inventory();hashes={};components={}
    for case,names in inventories.items():
        actual={p.name for p in (ROOT/case/'certificates').glob('*.json')}
        require(actual==names,f'Incorrect inventory for {case}')
        for name in sorted(names):
            p=ROOT/case/'certificates'/name
            hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
        print(f'Verifying {case}: {len(names)} certificates',file=sys.stderr,flush=True)
        components[case]=run_json(f'{case}/'+('verify.py' if case=='large_a' else 'verify_all.py'))
    totals=(sum(map(len,inventories.values())),sum(x['nodes'] for x in components.values()),sum(x['leaves'] for x in components.values()))
    require(totals==(110,50830,25470),'Unexpected certificate/node/leaf totals')
    zeros=sum(x.get('strict_zero_leaves',0) for x in components.values())
    require(zeros==381,'Unexpected strict zero-RHS count')
    audits={}
    for path in ('large_a/audit_identity.py','case21/audit_identities.py','case12/audit_identities.py'):
        print(f'Checking coordinate identities: {path}',file=sys.stderr,flush=True)
        audits[path]=run_json(path)
    return {'result':'All 110 certificates verified exactly.',
      'scope':'Finite certificates and scalar/coordinate audits for gamma8_consolidated.pdf. Analytic inclusions, classical inputs and BKR Theorems 4–5 are stated in the manuscript and are not formalized here.',
      'target':'0 < Q <= (8 |det Q|)^(1/5)',
      'certificate_count':totals[0],'nodes':totals[1],'leaves':totals[2],
      'strict_zero_rhs_leaves':zeros,'global_scalar_audit':scalar_audit(),
      'supplementary_audits':audits,'certificate_sha256':hashes,'components':components,
      'python_version':sys.version.split()[0],'elapsed_seconds':round(time.monotonic()-start,3)}

if __name__=='__main__':print(json.dumps(run(),indent=2))
