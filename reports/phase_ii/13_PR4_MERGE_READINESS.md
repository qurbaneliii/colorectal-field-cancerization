# PR #4 merge readiness

**READY FOR OWNER-APPROVED MERGE.** PR #4 is open, out of draft and mergeable into main. The owner retains the merge decision. This Codex audit carries no independent human approval.

| Gate | Result |
| --- | --- |
| Executive revision | `835263e526053d2043daf04eeb2212863489a047` |
| Hosted current-GEO pipeline | PASS · 20 stages and 54 tests |
| Push and PR CI | PASS |
| Artifacts and hashes | PASS · 196 generated files, 10 manuscript tables, 18 manuscript sections and six model cards checked |
| Scientific evidence replay | PASS · 45 outer folds, 5,310 held-out predictions, 1,465 metric, count and interval records |
| Critical defects | No confirmed P0 remains. SVM calibration P1 corrected, tested and rerun. |
| Repository payload | PASS · no raw archives, CELs, model binaries or Parquet files in the remediation. Locks and both historical TAR pins unchanged. |

The mean outer-fold SVM macro-F1 changed from 0.961987 to 0.970519 for Task A, 0.997558 to 1.000000 for Task B, and 0.989811 to 0.988137 for Task C. The corrected scores are comparator results. The primary Elastic Net and random forest fold metrics, split assignments, Task C lock and external evaluations reproduced within numerical tolerance.

The Task B acquisition alias remains a P1 scientific limitation. Historical byte recovery, inductive preprocessing and matched-cohort validation remain research tasks. Publication and clinical readiness have separate evidence requirements.

[PR #4](https://github.com/qurbaneliii/colorectal-field-cancerization/pull/4) · [Exact-HEAD workflow](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/43846818826) · [Artifact bundle](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/43846818826/artifacts/14221831791) · [CI](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/43846822896)
