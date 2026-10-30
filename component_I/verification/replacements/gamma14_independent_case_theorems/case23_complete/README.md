# Gamma_(1,4)=8 — new independent replacement for (m,K)=(2,3)

**Scope**: This package gives a proposed independent computer-assisted proof of the *last BKR case-theorem dependency*, the canonical case (m,K)=(2,3). It does **not** by itself verify the original manuscript's 110 other rational certificates, the classical reduction, or the independent replacements for the four other imported pairs. Independent expert review is still required.

## Files

- `case23_independent_proof.pdf`, `.tex`: complete mathematical reduction through L=1,...,6.
- `verify_case23_complete.py`: Python-standard-library exact rational verifier. Rebuilds all five 2D residue-collision line arrangements, audits key rational bounds, checks every fixed witness on every dyadic rectangle, and checks tree completeness.
- `case23_L1_t_certificate.json`: t=1/2, first L=1 resonance (124 leaves).
- `case23_new_certificates.json`: L=1 at A=1/6 (2 models, 194 leaves), L=2 at t=1/4 (2 models, 164 leaves), L=3 at t=1/6 (2 models, 188 leaves).
- `case23_L56_certificate.json`: previous L=5,6 third-square H certificates (140 leaves).
- `build_case23_resonances.py`: OPTIONAL floating-point-aided discovery of six models. Not imported by verification.
- `verification_result.json`: an example full successful exact-checker result.
- `SHA256SUMS.txt`: hashes of final publication data, including proof and checker.

## Run

```bash
python3 -S -B verify_case23_complete.py > my_verification_result.json
```

Expected totals: five circular-gap models, nine quadtree certificates, 810 terminal rectangles, 267 internal nodes. Python 3.11+ standard library only. Numerical discovery and any optimization packages are unnecessary for verification.

### What the verifier guarantees

- Each arrangement includes all congruence/collision lines intersecting its exact rational rectangle and all pairwise intersections with the boundary. It evaluates the maximum circular gap exactly at every arrangement vertex.
- Each leaf's witness is valid over the **entire** dyadic rectangle and scalar mu interval (strictly positive and strictly below target).
- Every nonterminal quadtree node has exactly four children and the full unit square is covered.
- Scalar comparisons used to enclose the physical coefficient intervals are exact.

### What it does not guarantee automatically

The checker does not mechanically formalize all analytic steps in the attached paper: the canonical reduction, the application and resonant congruence classification of Dumir's/Macbeath's lemmas, the original full-case split, the identity of F/G/H with the original form, or the other 110 certificates. Those are described mathematically and require independent audit.

## Result inventory

| Branch | Method | Certificate leaves |
|---|---|---:|
| L=1, t=1/2 class (0,1/2) | G cover | 124 |
| L=1, t=1/2 class (1/2,0) | A-circle-gap/Macbeath, then two F covers | 106+88 |
| L=2 | Macbeath, then two G covers | 76+88 |
| L=3 | t-circle-gap/Macbeath, then two G covers | 124+64 |
| L=4 | Uniform nine-residue circular-gap theorem | 0 |
| L=5,6 | Previously verified H covers | 82+58 |
| **Total** | | **810** |

The earlier report for L=5,6 already identified a genuine obstruction to using the **last completed form H alone** at L=4. Retaining x5 in G resolves L=4 without contradiction.
