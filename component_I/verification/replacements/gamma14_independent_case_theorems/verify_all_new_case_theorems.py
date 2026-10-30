#!/usr/bin/env python3
"""Run all exact checks for the FIVE imported (m,K) case replacements.

Scope: previous four case-replacement checks (including partial older L5/6)
plus unified (2,3) L=1..6 check. The old partial L5/6 check is intentionally
repeated as a cross-check; unique covering-rectangle counts exclude duplicate.

This does NOT verify gamma14.pdf's 110 original certificates or formalize
its canonical reductions and analytic inclusions.
"""
import pathlib,subprocess,sys,json
ROOT=pathlib.Path(__file__).resolve().parent
commands=[
    ('Prior (4,1),(3,1),(3,2),(2,4) and L5/6',ROOT/'previous_case_checks'/'verify_all_new.py'),
    ('Complete (2,3), L=1..6',ROOT/'case23_complete'/'verify_case23_complete.py')
]
for name,script in commands:
    proc=subprocess.run([sys.executable,'-S','-B',str(script)],cwd=script.parent,
                        text=True,capture_output=True)
    if proc.returncode:
        print(proc.stdout)
        print(proc.stderr,file=sys.stderr)
        raise RuntimeError('FAIL at '+name)
    print('PASS:',name)
    if 'Complete (2,3)' in name:
        out=json.loads(proc.stdout)
        assert out['status']=='PASS' and out['totals']['terminal_rectangles']==810
        print('  ',out['totals'])
    else:
        assert 'PASS: 4 verification runs' in proc.stdout
print('PASS ALL AVAILABLE CASE-THEOREM FINITE CHECKS')
print('Unique terminal rectangles in new BKR case replacements: 312 + 202 + 214 + 810 = 1538.')
print('LIMITATION: gamma14 original 110 certificates and analytic reductions NOT checked here.')
