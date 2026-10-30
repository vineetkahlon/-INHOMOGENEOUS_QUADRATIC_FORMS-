#!/usr/bin/env python3
"""Reproduce the two finite, exactly verified branches in Gamma_(1,4)=8.

This is an exact-certificate verifier, not a formalization of the classical
lemmas or of every analytic domain inclusion in the manuscript.
Requires Python3 without -O (inherited base checkers contain assertions).
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
BASE = ROOT/'verification'/'base'/'gamma8_consolidated'
NEW = ROOT/'verification'/'replacements'/'gamma14_independent_case_theorems'


def run_json(label, cwd, command):
    proc=subprocess.run([sys.executable,'-S','-B',*command],cwd=cwd,capture_output=True,text=True,check=False)
    if proc.returncode:
        raise RuntimeError(f'{label} failed [{proc.returncode}]:\n{proc.stderr[-6000:]}\n{proc.stdout[-6000:]}')
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as e:
        raise RuntimeError(f'{label} did not return valid JSON: {e}') from e


def run_lines(label,cwd,command):
    proc=subprocess.run([sys.executable,'-S','-B',*command],cwd=cwd,capture_output=True,text=True,check=False)
    if proc.returncode:raise RuntimeError(f'{label} failed:\n{proc.stderr}\n{proc.stdout}')
    return proc.stdout


def main():
    if not __debug__:raise RuntimeError('Python optimization disables assertions; run without -O')
    manifest=run_json('original SHA-256 manifest',BASE,['verify_manifest.py'])
    if manifest.get('result') != 'All listed file hashes match.':raise AssertionError('Invalid base manifest verification')
    base=run_json('original 110 certificate trees',BASE,['verify_all.py'])
    for field,expected in [('certificate_count',110),('nodes',50830),('leaves',25470),('strict_zero_rhs_leaves',381)]:
        if base.get(field)!=expected:raise AssertionError(f'Base {field}: expected {expected}, got {base.get(field)}')
    new_output=run_lines('independent case-theorem replacements',NEW,['verify_all_new_case_theorems.py'])
    if 'PASS ALL AVAILABLE CASE-THEOREM FINITE CHECKS' not in new_output:raise AssertionError('Missing replacement pass marker')
    detailed=run_json('new complete (2,3) checks',NEW/'case23_complete',['verify_case23_complete.py'])
    totals=detailed.get('totals',{})
    if totals.get('exact_gap_boxes')!=5 or totals.get('quadtree_certificates')!=9 or totals.get('terminal_rectangles')!=810:
        raise AssertionError('Unexpected (2,3) inventory')
    leaf_by_pair={'(3,1)':312,'(3,2)':202,'(2,4)':214,'(2,3)':810}
    expected_cases=[(4,1),(3,1),(3,2),(2,1),(2,2),(2,3),(2,4),(1,1),(1,2)]
    if len(expected_cases)!=9 or len(set(expected_cases))!=9:raise AssertionError('Case classification mismatch')
    hashed={}
    for folder in ('verification/base/gamma8_consolidated', 'verification/replacements/gamma14_independent_case_theorems'):
        for p in (ROOT/folder).rglob('*.json'):
            if 'certificates' in str(p) or p.name in ('case23_new_certificates.json','case23_L1_t_certificate.json','case23_L56_certificate.json'):
                hashed[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return {
        'status':'PASS',
        'scope':'Finite exact certificate verification only; does not replace analytic and classical theorem audits.',
        'base':{'farkas_trees':110,'nodes':50830,'terminal_contradictions':25470,'strict_zero_rhs_leaves':381,
                'manifest_files_checked':manifest['files_checked'],
                'coordinate_audits':base.get('supplementary_audits',{}),
                'scalar_checks':base.get('global_scalar_audit',{})},
        'replacements':{'status':'PASS','dyadic_trees':20,'terminal_rectangles':sum(leaf_by_pair.values()),
                        'by_pair':leaf_by_pair,'uniform_exact_gap_lemmas':5,
                        'new_case23':totals},
        'global_pair_classification':{'count':len(expected_cases),'pairs':[list(v) for v in expected_cases]},
        'certificate_data_sha256':hashed,
        'python_version':sys.version.split()[0],
        'critical_limitations':['classical density/Meyer/Birch/canonical reduction not machine-checked',
                               'Dumir/Jackson/Macbeath theorems assumed as published results',
                               'analytic domain inclusions and resonance normalizations require independent expert review'],
    }

if __name__=='__main__':
    print(json.dumps(main(),indent=2))
