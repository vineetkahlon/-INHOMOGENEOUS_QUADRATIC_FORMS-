#!/usr/bin/env python3
"""Reproduce exact finite verification of the case-replacement supplements.

Run: python3 -S -B verify_all_new.py
All certificate checks are exact rational with Python standard library only.
Proof acceptance does NOT rely on running the discovery code.
"""
import sys, subprocess, json, pathlib
ROOT=pathlib.Path(__file__).resolve().parent
PAIRS=[('verify_bkr_case31.py','case31_resonance_certificate.json'),
       ('verify_case32_cover.py','case32_certificate.json'),
       ('verify_case24_cover.py','case24_certificate.json'),
       ('verify_case23_L56.py','case23_L56_certificate.json')]
for script,cert in PAIRS:
    cmd=[sys.executable,'-S','-B',str(ROOT/'verification'/script),str(ROOT/'certificates'/cert)]
    p=subprocess.run(cmd,cwd=str(ROOT/'verification'),text=True,capture_output=True)
    if p.returncode:
        print(p.stdout);print(p.stderr,file=sys.stderr)
        raise RuntimeError('FAIL: '+script)
    if script=='verify_bkr_case31.py':
        assert 'PASS:' in p.stdout, 'missing case31 PASS marker'
    else:
        assert json.loads(p.stdout)['status']=='PASS', 'missing structured PASS marker'
    print('PASS',script,'<--',cert)
print('PASS: 4 verification runs; 4 of 5 imported cases independently replaced, plus 2 branches of the fifth.')
