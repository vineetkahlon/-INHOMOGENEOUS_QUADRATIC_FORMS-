ACTUAL PYTHON CODE AND CERTIFICATE DATA
======================================

This package contains the actual source files and certificate JSON files,
not just a SHA-256 manifest or verification report.

Version: the earlier October 26, 2026 proof package.
It is NOT the newer 28-Farkas-tree / 98-dyadic-cover release described in
manuscript gamma14-4.pdf.

CONTENTS
- 110 Farkas trees, with 50,830 nodes and 25,470 terminal contradictions.
- 20 distinct dyadic witness trees, with 1,538 terminal rectangles.
  Several trees are grouped together in individual JSON files; an older
  partial-case file is also retained as a duplicate cross-check.
- Actual Python model generators and exact checkers.
- Optional discovery programs and supporting documents retained from the
  original verification folders to preserve the internal hash manifests.

RUN
From this directory:
    python3 -S -B verify_end_to_end.py > full_verification.json
Do not use -O or -OO. Exact verification uses only Python's standard library.
Optional discovery programs can require other packages; they are not needed
for the verification command.

ACTUAL CERTIFICATE FILE EXAMPLES
verification/base/gamma8_consolidated/case22/certificates/L2.json
verification/base/gamma8_consolidated/case12/certificates/L2_pp.json
verification/base/gamma8_consolidated/large_a/certificates/family_31.json
verification/replacements/gamma14_independent_case_theorems/previous_case_checks/certificates/case31_resonance_certificate.json

The JSON files above contain tree nodes, witnesses, branching data and/or
rational infeasibility certificates. Files named SHA256SUMS or manifest
are supplemental integrity information, not substitutes for those data.

SCOPE
A successful run verifies the encoded finite computations, not every
analytic reduction or the classical theorems cited by the manuscript.

REPACKAGING
The verifier and verification-directory files were extracted byte-for-byte
from gamma14_end_to_end_verified_2026-10-26.zip. This root README and the
fresh verification output were added; no verifier or certificate was changed.
