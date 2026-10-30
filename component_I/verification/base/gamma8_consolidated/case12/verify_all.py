"""Exact verification of the full (1,2) certificate cover; no third-party imports."""
import hashlib,json
from fractions import Fraction as F
from pathlib import Path
from model import d_bounds,SQRT43_UPPER
from binary_resonance_model import d_bounds as resonance_bounds
from verify_cover import verify as verify_ternary
from verify_four import verify as verify_four
from verify_resonance import verify as verify_resonance

def require(condition,message):
    if not condition:raise ValueError(message)

def run():
    root=Path(__file__).resolve().parent;records=[]
    require(F(413,1000)**3>F(587,1000)**5,'Invalid a bound')
    require(SQRT43_UPPER**2>F(4,3),'Invalid root bound')
    for L in (1,2,3,4):
        lo,hi=d_bounds(L)
        require(lo**2<=max(F(1,8),F(123,256*(L+1))),'Invalid lower d bound')
        if L>1:
            require(27*(L*L+3)**3*hi**5>=8192*(1-hi)**3,'Invalid upper d bound')
    require(F(29,50)**4<F(1,8),'Invalid four-family lower endpoint')
    require(F(7,12)**4<F(1,8),'Four-family target not below d')
    lo,hi=resonance_bounds()
    require(lo**3<F(1,16),'Invalid resonance lower d bound')
    require(hi**10>=(F(2,27))**3,'Invalid resonance upper d bound')

    cases=[]
    for L in (2,4):
        for sign in ('pp','pn'):cases.append((f'L{L}_{sign}',L,None,sign,(F(0),F(1,2))))
    for region,L in [('L1_outer',1),('L1_A_low',1),('L1_A_high',1),('L3_high',3)]:
        for sign in ('pp','pn'):cases.append((f'{region}_{sign}',L,region,sign,(F(0),F(1,2))))
    for sign in ('pp','pn'):
        for tag,lb in [('a',(F(0),F(1,4))),('b',(F(1,4),F(1,2)))]:
            cases.append((f'L3_low_{sign}{tag}',3,'L3_low',sign,lb))
    def load(name):
        path=root/'certificates'/f'{name}.json';raw=path.read_bytes();proof=json.loads(raw)
        require(not proof.get('pending'),'Pending proof branches')
        return proof,{'file':str(path.relative_to(root)),'sha256':hashlib.sha256(raw).hexdigest()}
    for name,L,region,sign,lb in cases:
        proof,record=load(name)
        expected=[(F(0),F(1,2)),(F(0),F(1,2)) if sign=='pp' else (F(-1,2),F(0)),lb]
        require(proof['L']==L and proof.get('region')==region,'Wrong ternary region')
        require([tuple(map(F,p)) for p in proof['shape_box']]==expected,'Incomplete shear cover')
        require('d_interval' not in proof,'Narrowed d interval')
        records.append({**record,'kind':'ternary',**verify_ternary(proof)})
    for family in range(4):
        proof,record=load(f'four_{family}')
        require(proof['family']==family,'Wrong four-family case')
        records.append({**record,'kind':'four-family',**verify_four(proof)})
    for family in range(2):
        proof,record=load(f'resonance_{family}')
        require(proof['family']==family,'Wrong binary resonance case')
        require(all(k not in proof for k in ('shape_box','d_interval','region')),'Narrowed resonance domain')
        records.append({**record,'kind':'binary resonance',**verify_resonance(proof)})
    return {'result':'All 22 certificates and scalar comparisons verified exactly.',
            'scope':'Reduced (m,K)=(1,2), with the analytic reductions in gamma8_consolidated.pdf, Section 9.',
            'nodes':sum(r['nodes'] for r in records),'leaves':sum(r['leaves'] for r in records),
            'strict_zero_leaves':sum(r.get('strict_zero_leaves',0) for r in records),
            'certificates':records}

if __name__=='__main__':print(json.dumps(run(),indent=2))
