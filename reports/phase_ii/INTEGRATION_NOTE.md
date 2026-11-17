# Phase III integration of archived Phase II evidence

PR #4 was merged with owner authorization on 2026-11-17. Merge commit
`9f721db1d498444996df13f54cb6a478554ca83b` has the same tree as audited
`835263e526053d2043daf04eeb2212863489a047`. GitHub showed no applicable branch
rules, no required approvals, no blocking reviews and successful exact-head CI
and hosted reproduction before the merge. No protections were changed.

The integration branch starts at that merge and cherry-picks, with source trailers,
the four unique Phase II commits `f71f2f8`, `cbff26a`, `e191937`, and `b001260`.
The historical source branch remains unchanged at
`b0012604ee9462d5064d09f622ad37dbc11e26f5`. There were no content conflicts.
The 28 transferred files comprise 15 reports, 10 evidence files, two scripts and
one test file. Reports describing PR #4 as open are historical snapshots.
Their results retain their actual generating commits and execution links.

Integration verification used the existing locked Python environment. `pip check`,
`compileall`, Ruff and 27 data-independent tests passed (28 full-data tests deselected).
`scripts/verify_phase_ii.py` against the retained exact-head hosted bundle verified
45 outer partitions, 5,310 predictions and 1,465 metric/count/interval records.
The existing scientific pipeline and locks are unchanged, so an expensive raw-CEL
rerun is not needed for this evidence-and-prototype integration.

The prototype has synthetic regression evidence. It has no acquisition-matched
empirical validation. Issues #6 and #7 remain open. This integration requires
review and is not authorized for automatic merge.
