# Γ₁,₄ = 8 — verification package, component II

Exact-arithmetic verification of every finite claim that is **new** in the manuscript
*Positive values of inhomogeneous quadratic forms of type (1,4): the constant Γ₁,₄ = 8*.

## Run

```
python3 verify_all.py > inventory.json
```

Python 3 standard library only (`fractions`, `json`, `hashlib`). numpy and scipy are never imported by the
verifier. Running it with both modules blocked gives the same result. Runtime is about 15 s. The script ends with
`"status": "ALL CHECKS PASSED"` or exits with an error naming the failing item.

## What is checked

| Paper | Item | Checker |
|---|---|---|
| Lemma 2.9 | two-piece gap bound (160 + 229 arrangement vertices) | `verify_twopiece.py` |
| Prop. 3.3, Table 2 | 72 large-a dyadic covers | `verify_all.py` → `cover.check` |
| Thm 4.3, Tables 4–5 | (3,1) interval tests, residue gaps, n=2 residues for the 8 classes | `verify_31.py` |
| Thm 4.3 | 1 (3,1) dyadic cover (the all-half class) | `cover.check` |
| Lemma 5.2 | uniform gap rows U1–U5 (exact line arrangement) | `verify_gaps.py` |
| Lemma 5.3, Table 6 | 11 shared resonant covers R0a–R9 | `cover.check` |
| §6.2 | 4 (2,2) L=1 covers | `cover.check` |
| §6.4 | (2,2) L=2,3 exact t-partitions (listed in `PARTITIONS.md`) | `verify_partitions.py` |
| §6.4, Table 7 | 6 (2,2) resonant covers (one with pinch-ball leaves) | `cover.check` |
| Lemma 2.15, Thm 6.1 | ladder data at the three pinch points: rungs, polynomial rung families, ratio monotonicity, handoff inequalities, radii | `verify_pinch.py` |
| §6.4, App. A.4 | 3 Farkas trees at the exact pinch points, roots regenerated | `verify_corners_exact.py` |
| Thm 8.3 | 4 (1,2) L=1 covers | `cover.check` |
| Supplement (six variables) | pruning thresholds, (3,1) scalars, mod-2 residues | `verify_sixvar.py` |
| throughout | 214 scalar comparisons and congruences from the text | `verify_scalars.py` |

`specs.py` is the single source of truth for every fixed family. `verify_all.py` refuses to run if
the certificate directory and the specification differ.

## Files

- `specs.py`, `families.py`: family definitions (witness constants, μ-intervals, targets, neighbourhoods).
- `cover.py`: dyadic-cover builder and the exact checker `check`.
- `gap.py`, `partition.py`: exact maximal circular gap (Prop. 2.8) and t-partition helpers.
- `corner_model.py`, `farkas_check.py`: corner root polytopes (App. A.4) and LP-free Farkas replay.
- `certs/`: all certificates (JSON). `inventory.json`: the output of the last full run (per-certificate counts and hashes).
- `PARTITIONS.md`: the exact t-partitions of §6.4.
- `retired/`: the seven (3,1) covers made unnecessary by the n=2 residue argument (kept for reference, not checked).
- `discovery/`: the tools that *found* the certificates (`corner22.py` and `farkas.py` use scipy's LP).
  They are not needed for verification.

## Component I

25 Farkas trees for (2,1), (1,1) and (1,2) (23,753 nodes, 11,889 leaves) are retained unchanged from the
archive of the earlier unified argument (`verify_end_to_end.py` there). Superseded and no longer used: the 72
large-a trees, the seven (2,2) trees, the four (1,2) fixed-family trees, and the two (1,2) L=4 ternary trees
`case12/certificates/L4_pp.json` and `L4_pn.json` (now covered by Lemma 2.10, the high-quotient closure).

**Before publication**, merge component I into this release with a manifest that lists exactly the 25 retained
trees and their hashes, so that the master command checks both components and rejects missing or extra files.

## Trust base

Exact rational arithmetic. Floating point appears only as a seed for rational cube-root lower bounds,
which are then certified by exact cubing, and in the discovery tools. The inclusion of the true parameter
sets in the certified domains is proved in the manuscript and is not machine-checked.
