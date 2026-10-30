"""Verify the complete seven-certificate cover, using only standard Python."""
import hashlib,json
from fractions import Fraction as F
from pathlib import Path
from verify_ternary import verify as verify_ternary
from verify_four import verify as verify_four
from model import d_bounds

def run():
    root=Path(__file__).resolve().parent
    # Scalar comparisons used in the analytic L=1 reduction.
    assert F(125,576)>F(32,75)**2
    assert F(18,35)**3<F(5,36)
    assert F(7,10)**4<F(5,4)**4/F(10)
    assert F(8,15)**4<F(1,12)
    assert F(1,12)>F(1,2)**4
    assert F(99,70)**2>2
    for L in (1,2,3,4):
        lo,hi=d_bounds(L)
        assert lo>0 and lo*lo<=F(125,288*(L+1))
        assert hi==F(64,25*(L+3))
    records=[]
    for L in (2,3,4):
        path=root/'certificates'/f'L{L}.json';raw=path.read_bytes();proof=json.loads(raw)
        assert proof['L']==L and proof['model_version']==3 and 'd_interval' not in proof
        assert not proof['pending']
        result=verify_ternary(proof)
        records.append({'file':str(path.relative_to(root)),'sha256':hashlib.sha256(raw).hexdigest(),**result})
    for family in range(4):
        path=root/'certificates'/f'family_{family}.json';raw=path.read_bytes();proof=json.loads(raw)
        assert proof['family']==family and not proof['pending']
        result=verify_four(proof)
        records.append({'file':str(path.relative_to(root)),'sha256':hashlib.sha256(raw).hexdigest(),**result})
    return {'result':'All seven certificates and scalar comparisons verified exactly.',
            'scope':'Reduced (m,K)=(2,2), using the BKR reduction and Macbeath lemma as mathematical inputs.',
            'nodes':sum(r['nodes'] for r in records),'leaves':sum(r['leaves'] for r in records),'certificates':records}

if __name__=='__main__':print(json.dumps(run(),indent=2))
