"""Verify all five required certificates and scalar comparisons, exactly."""
from fractions import Fraction as F
from pathlib import Path
import json
from certificate import EXPECTED,verify

def scalar_checks():
    assert F(4419,10000)**2 < F(25,128)  # lower endpoint < 5 sqrt(2)/16
    assert F(1414214,1000000)**2 > 2
    assert F(6389,10000)**4 < F(1,6)
    assert F(4,7)**4 < F(1,6)
    assert F(3,5)**4 < F(1,6)  # delta >= (5/4) 6^(-1/4) > 3/4
    assert 3*F(3,4)**4==F(243,256)<1
    assert 2*F(7,10)**5+F(7,10)>1
    assert 2*F(11,16)**5<F(5,16)
    return 8

def main():
    directory=Path(__file__).resolve().parent/'certificates'
    if {p.name for p in directory.glob('*.json')}!=set(EXPECTED):raise ValueError('Wrong certificate inventory')
    results=[]
    for filename,domain in EXPECTED.items():
        proof=json.loads((directory/filename).read_text())
        if proof['domain']!=domain:raise ValueError('Wrong domain assignment')
        results.append(verify(proof))
    summary={'certificates':len(results),'nodes':sum(r['nodes'] for r in results),'leaves':sum(r['leaves'] for r in results),
             'scalar_checks':scalar_checks(),'target':'0<Q<d','verified_exactly':True,'results':results}
    print(json.dumps(summary,indent=2))
    return summary
if __name__=='__main__':main()
