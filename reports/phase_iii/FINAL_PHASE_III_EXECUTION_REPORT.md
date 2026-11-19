# Phase III execution report

This report separates completed engineering, measured conditional experiments and
unresolved independent validation. No Codex review is independent human approval.

| Workstream | Status | Main evidence | Remaining dependency |
|---|---|---|---|
| A. PR #4 merge | PASS | Merge `9f721db1d498444996df13f54cb6a478554ca83b`; main CI 44026475308 | None for authorized merge |
| B. Phase II integration | PASS | PR #8, `56d7cf0dfc575094e9a1c40d621f6816c97530e0`; CI 44026875052 | Owner review and acceptance; remains unmerged |
| C. Frozen preprocessing | PARTIAL | PR #10; measured fixed-panel experiment and independent replay | Full nested discovery, stringent sensitivity and established-method comparison remain NOT RUN |
| D. Independent validation | PARTIAL / BLOCKED | PR #9, `9249314143cae30c4cf5c398fa73036616edc7e7`; CI 44028635372 | Institutional access, approved compute, eligible data and assay transport |

## Integration evidence

The owner-authorized merge occurred at 2026-11-17 05:46:47 UTC. Audited PR HEAD
`835263e526053d2043daf04eeb2212863489a047` is the second parent of the merge and
the merged tree is identical. The previous main was
`c2a36a1d2f80fd25e1565b543c909eebef64ca09`. There were no applicable branch rules,
required approvals, blocking reviews or unresolved review threads at merge time.
No protection was bypassed. Exact-head hosted evidence remains run 43846818826
(20 stages / 54 tests). Main's post-merge CI passed.

Four unique Phase II commits from `b0012604ee9462d5064d09f622ad37dbc11e26f5`
were cherry-picked with source trailers onto the merge. There were no content
conflicts. The original branch is preserved. Integration transferred 28 files
plus an integration note. Dependencies, Ruff, syntax, 27 unit/contract tests,
synthetic smoke and independent replay of 45 partitions, 5,310 predictions and
1,465 records passed. Historical evidence retains its generating revisions.

## Research execution

The frozen experiment's exact analysis revision is recorded in
`results/phase_iii/frozen_preprocessing/execution.json`. Source TAR pins, compressed
CEL member digests, baseline matrix hashes, mapping hashes, patient partitions,
reference arrays and numeric models are retained. The R environment is 4.5.1
with preserved dependencies; Python uses the existing lock. Primary processing,
historical hashes and published estimates were not changed.

The new method is a minimal training-quantile/median-PM comparator. It is not
fRMA and changes both normalization and summarization. Fixed historical panels
and prespecified L2 logistic regression make this a conditional sensitivity,
not an unbiased reproduction of complete nested feature discovery. All 246
development CELs and 233 external endpoint CELs were processed independently.
Detailed measured results, invariance, calibration flags and paired donor
intervals are in the frozen-preprocessing reports and machine-readable outputs.

No external outcome selected parameters or thresholds. GSE41258 remains previously
exposed retrospective evidence. The experimental two-gene transport refit and its
historical reference cutoff 0.41 do not replace the primary Elastic Net result.

The initial comparison at `b00e552` failed when the archived external baseline's
390 rows were joined to 233 eligible endpoint records. The error stopped the run.
Commit `cb407ac0069200fce57e9a046ee667069b6ebf8e` applies the same eligibility/order
to both arms, adds a regression test and reruns from the unchanged CEL cache.
The failed attempt is retained as failure evidence, not counted as a successful
scientific run. The corrected branch passes 32 data-independent tests, with 28
full-data tests deselected. Ruff and the locked dependency check pass.

Issue #6 preparation implements schema checks, explicit exclusions, donor
consistency, source-ordered biopsy choice, source preservation and synthetic-only
execution by default. The suite passed 36 tests; 10 targeted tests passed after
final eligibility edits. A synthetic CLI selected four donors. Eight reports
cover public-source research, eligibility, access, analysis, platform mapping,
readiness and limitations. Their branch depends on PR #8, as does PR #10.

## Research boundaries

Historical TAR equality remains FAIL and mismatch attribution remains BLOCKED
under Issue #5. Missing historic bytes were not replaced. Task B case status and
acquisition remain aliased. No preprocessing change resolves that limitation.
JGAD000785 requires qualified investigator/institutional details, ethics and
institutional authorization, Type I approval, approved compute and eligible
donor/sample records. No restricted data were accessed or application submitted.
SNORA12 coverage and compatible microarray-to-RNA-seq transport remain unresolved.
Clinical validity, independent human review and author declarations remain open.

## GitHub evidence

- [Merged audit PR #4](https://github.com/qurbaneliii/colorectal-field-cancerization/pull/4)
- [Main CI](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/44026475308)
- [Phase II integration PR #8](https://github.com/qurbaneliii/colorectal-field-cancerization/pull/8)
- [Integration CI](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/44026875052)
- [Controlled-access preparation PR #9](https://github.com/qurbaneliii/colorectal-field-cancerization/pull/9)
- [Preparation CI](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/44028635372)
- [Frozen real-CEL experiment PR #10](https://github.com/qurbaneliii/colorectal-field-cancerization/pull/10)

Only PR #4 was authorized for merge. New PRs remain open for owner review.
