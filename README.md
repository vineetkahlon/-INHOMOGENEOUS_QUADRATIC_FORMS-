# Γ₁,₄ = 8 — complete verification release

Finite data for *Positive values of inhomogeneous quadratic forms of type (1,4): the constant Γ₁,₄ = 8*.
DOI: [to be added]   ·   Paper: [arXiv id / version to be added] (the manuscript is not included in this repository)

## Contents

| Folder | What | Verifier |
|---|---|---|
| `component_I/` | the earlier unified-argument archive, with its manuscript drafts removed (see `component_I/STRIPPED_FILES.json`); certificates, model generators and checkers are unchanged; the paper uses 25 of its Farkas trees (23,753 nodes) | `verify_end_to_end.py` (its own) |
| `component_II/` | new data: 98 dyadic covers, 3 Farkas trees at the pinch points, gap/partition/ladder/scalar checks | `verify_all.py` |
| `component_I_usage.json` | the 25 component-I trees the paper uses (by lemma), with paths and pinned SHA-256 hashes; superseded files | `make_manifest.py`, `verify_release.py` |
| `component_I_reference_replay_2026-10-29.json` | reference report of a full replay of component I (110 trees, 50,830 nodes, all certificate hashes) | compared by `verify_release.py` |
| `MANIFEST.json`, `SHA256SUMS` | hashes of every file | `verify_release.py` |

## Verify

    python3 verify_release.py

checks every hash against the manifest (and rejects unlisted files) and the 25 used component-I trees
against their pinned hashes, replays component II, then replays component I and requires its report to
match the 2026-10-29 reference (tree and node counts, and the hash of every certificate file). Python 3 standard library only; no floating-point LP is used for acceptance.



