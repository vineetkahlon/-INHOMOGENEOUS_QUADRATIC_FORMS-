"""Verify the complete (1,1) cover using only exact standard-library arithmetic."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from model import d_bounds, SQRT43_UPPER
from verify_cover import verify


def require(condition, message):
    if not condition:
        raise ValueError(message)


def run():
    # Exact comparisons used to enclose the analytic domain.
    require(F(433,1000)**2 < F(3,16), 'Lower d endpoint is not outward')
    hi = F(629,1000)
    require(3*hi**3+2*hi-2 > 0, 'Upper d endpoint is not outward')
    require(F(413,1000)**3 > F(587,1000)**5, 'Invalid a upper bound')
    require(SQRT43_UPPER**2 > F(4,3), 'Invalid square-root upper bound')
    require(d_bounds(1) == (F(1,2),F(2,3)), 'Incorrect L=1 interval')
    require(d_bounds(2) == (F(433,1000),hi), 'Incorrect L=2 interval')

    root = Path(__file__).resolve().parent
    records = []
    for L in (1,2):
        for sign in ('pp','pn'):
            path = root/'certificates'/f'L{L}_{sign}.json'
            raw = path.read_bytes()
            proof = json.loads(raw)
            expected_shape = [(F(0),F(1,2)),
                              (F(0),F(1,2)) if sign=='pp' else (F(-1,2),F(0)),
                              (F(0),F(1,2))]
            shape = [tuple(map(F,pair)) for pair in proof['shape_box']]
            require(proof['L']==L and proof.get('model_version')==2, 'Wrong case')
            require(shape==expected_shape, 'Incomplete shear coverage')
            require('d_interval' not in proof, 'Narrowed d interval')
            require(not proof.get('pending'), 'Unfinished search')
            result = verify(proof)
            records.append({'file':str(path.relative_to(root)),
                            'sha256':hashlib.sha256(raw).hexdigest(), **result})
    # h >= 0 is enough by the explicitly proved symmetry
    # (x,h,k,beta1) -> (-x,-h,-k,-beta1), leaving lambda unchanged.
    return {
        'result':'All four certificates and scalar comparisons verified exactly.',
        'scope':'Reduced (m,K)=(1,1) at constant eight, using the mathematical domain and sign symmetry in gamma8_consolidated.pdf, Section 8.',
        'target':'0 < Q <= (8|det Q|)^(1/5)',
        'nodes':sum(r['nodes'] for r in records),
        'leaves':sum(r['leaves'] for r in records),
        'strict_zero_leaves':sum(r['strict_zero_leaves'] for r in records),
        'certificates':records}


if __name__=='__main__':
    print(json.dumps(run(),indent=2))
