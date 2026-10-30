# Constant eight: consolidated manuscript and verification bundle

Start with **gamma8_consolidated.pdf**. This 21-page manuscript contains the
complete assembled argument for `Gamma_(1,4) = 8`, with its remaining published
inputs stated explicitly. The LaTeX source is `gamma8_consolidated.tex`.
No earlier project note or code package is needed to read the new arguments
or run their exact verification.

The proof retains BKR Theorems 4–5 (the five reduced pairs in Sections 6–7)
as cited inputs, along with the classical reduction and gap lemmas. Their
proofs are not reproduced or machine-checked here. All replacement arguments
for the bound portions of BKR Sections 4–5 and 8–11 are integrated into this
manuscript. This is a computer-assisted proof draft for independent review,
not a proof-assistant formalization or a new classification of equality cases.

## Verify the complete certificate component

From this directory:

```bash
python3 -S -B verify_all.py > my_verification_summary.json
```

Use Python 3 without `-O` or `-OO`; the master command rejects optimized mode
because some component checkers use assertions. Only the standard library is
needed. Verification requires no network, NumPy, SciPy or LP solver.

The command checks the complete inventory and prescribed parameter domains,
all **110 certificates, 50,830 nodes and 25,470 terminal contradictions**,
the global scalar comparisons, and the sharpness example. It also runs the
three supplementary coordinate/domain and malformed-proof audits. There are
381 valid zero-right-hand-side contradictions using a strict inequality.
Those are exact endpoint checks, not numerical tolerances.

`verification_summary.json` records a successful run, the individual
certificate hashes, and all component results. Analytic domain inclusions
and the cited theorems are mathematical parts of the manuscript, not claims
that the software independently proves them.

## Layout

| Directory | Manuscript section | Certificates | Role |
|---|---:|---:|---|
| `large_a/` | 5 | 72 | Four fixed negative blocks and their forced shift classes |
| `case21/` | 6 | 5 | Entire (2,1) case, including the former Section 8 exception |
| `case22/` | 7 | 7 | Three ternary covers and four fixed families |
| `case11/` | 8 | 4 | Two L values and two shear signs, with the upper endpoint |
| `case12/` | 9 | 22 | Ternary covers, fixed families and coupled binary resonances |

The nonintegral-shift proof in Section 4 is analytic; its scalar checks are
included in the verification path. `proof_inventory.json` records this map,
the counts, targets and external dependencies.

Each component contains its exact model generators, certificates, checker,
and optional discovery programs. The master runs them in separate Python
processes to keep the local model imports unambiguous. The redundant old
exception-only certificate and the obsolete final-Section-11 certificate
are excluded.

## Inspect integrity and rebuild the manuscript

```bash
python3 -S -B verify_manifest.py
bash build_manuscript.sh
```

Rebuilding the manuscript requires a LaTeX installation with the packages
listed in its preamble. The manifest describes the delivered files; rebuilding
or regenerating files may change their hashes. It detects accidental changes,
but is not a digital signature or a substitute for certificate verification.

`provenance.json` records the source and copied-file hashes of the component
code and certificates. The only edits to copied verification code are prose
references pointing to the unified manuscript; the mathematical models and
certificate data are unchanged. The new master additionally checks the
small-minimum threshold and the complete global inventory.

## Optional discovery

Files named `build*.py` use NumPy and SciPy to discover rational certificates.
They are not run or imported during verification. For example, to regenerate
a small diagonal cover without altering the delivered data:

```bash
python3 case21/build.py --kind diagonal --family 2 --radius 2 --output /tmp/new_diagonal_cover.json
```

A new discovery tree may differ from the delivered tree. It must pass the
exact checker; an unsuccessful or time-limited numerical run proves nothing.

Source used by reference: Bhardwaj, Kathuria and Raka,
https://arxiv.org/abs/2402.18939, version 1.
