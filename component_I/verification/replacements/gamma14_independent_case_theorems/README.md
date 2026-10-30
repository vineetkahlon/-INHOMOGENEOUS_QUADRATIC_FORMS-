# The Gamma_(1,4)=8 case-theorem-independent research compendium

**Date:** 26 October 2026. **Status:** Complete proposed replacements for the five imported BKR case theorems, with exact finite verifiers. Independent mathematical review and reproduction of the *base* manuscript's 110 certificates are still needed before publication.

## Why this bundle exists

The 20-page `base/gamma14.pdf` proves the target Gamma_(1,4)=8 conditional on two imported specialized results (Bhardwaj–Kathuria–Raka Theorems 4 and 5) for five `(m,K)` pairs. The three supplements now provide proposed independent replacements for all five pairs. The base PDF has **not been edited in place**; read its Theorem 2.5 and references to imported BKR case theorems as replaced by the accompanying supplements. The longer original BKR case analysis is not needed by these proposed replacement arguments.

| Imported pair | Supplement | Independent finite verification |
|---|---|---|
| (4,1) | bkr_cases_41_31_note.pdf | Pure analytic replacement |
| (3,1) | bkr_cases_41_31_note.pdf | 312 dyadic rectangles |
| (3,2) | bkr_cases_32_24_note.pdf | 202 dyadic rectangles |
| (2,4) | bkr_cases_32_24_note.pdf | 214 dyadic rectangles |
| (2,3) | case23_independent_proof.pdf | 9 quadtree models, 810 rectangles + 5 exact circular-gap bounds |

The `(2,3)` coverage includes all L=1,...,6. The final-square obstruction for L=4 is overcome by retaining x5 in the prior square; the maximum circular gap of nine residues is exactly 1/2 on the entire relevant rational interval.

**Total unique extra terminal rectangles: 1,538**, in addition to the 110 original manuscript certificates that are *described*, but not included in this bundle.

## One-command exact verification of the new specialized-case replacements

```bash
python3 -S -B verify_all_new_case_theorems.py
```

This invokes the earlier case checks and the new full `(2,3)` verifier. It uses only the Python standard library. The checkers reconstruct all finite witness inequalities with exact rational arithmetic, check tree coverage, and verify stated scalar bounds. The circular-gap checker enumerates all rational arrangement vertices, not a sample grid. The previous `(2,3)` L=5,6 trees are deliberately run twice (once in each verifier) as a cross-check; the unique leaf inventory counts them only once.

**Important verification boundary:** These checkers do not automatically verify the foundational canonical reduction, the classical Dumir/Macbeath/Jackson lemmas, every analytic lifting and domain-inclusion argument in the supplements, or the 110 separate trees described in the base `gamma14.pdf` Appendix A. This is a research proof package for independent audit, not a formally checked or referee-approved result.

## Files

- `base/gamma14.pdf`: original 20-page base manuscript (unchanged).
- `supplements/`: three mathematical proof supplements (PDF and editable LaTeX sources).
- `previous_case_checks/`: certificate data, exact verifiers and optional discovery sources for the earlier four pairs plus `(2,3)` L=5,6.
- `case23_complete/`: full `(2,3)` theorem, new certificate data, complete exact verifier, and optional certificate generator.
- `gamma14_independent_case_theorems_compendium.pdf`: one navigable PDF with this status page, base manuscript and all three supplements.
- `verify_all_new_case_theorems.py`: master runner for all the available specialized-case finite checks.
- `SHA256SUMS.txt`: file-integrity hashes (not a digital signature).

## Recommended final steps

1. Have a mathematician check the canonical reduction, all Macbeath exception derivations, shear normalizations and Dumir liftings line by line.
2. Recover and independently run the **original 110 rational certificates** and model generators described in `gamma14.pdf`. They are a distinct, necessary component of its claimed global proof and are not shipped here.
3. Integrate the three supplements into the base LaTeX, remove the stale import of BKR Theorems 4/5, and revise the final theorem proof and case table.
4. Publish one immutable source/certificate archive and reproducible checksums with the final manuscript.
