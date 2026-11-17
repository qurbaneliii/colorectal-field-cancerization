# Independent validation status

| Deliverable | Status | Evidence or dependency |
|---|---|---|
| Public-source research and proposed protocol | PASS | Reports 01–05 and source register |
| Eligibility implementation | PASS | Explicit schema, exclusions and donor rules |
| Synthetic tests | PASS | 36 data-independent tests; 28 full-data tests deselected; metadata CLI selected four synthetic donors |
| JGAD000785 authorization | BLOCKED | Investigator, institution, ethics, approval and approved compute |
| Assay/model portability | BLOCKED | Full panel detection and locked microarray-to-RNA-seq representation |
| Independent empirical Task B validation | NOT RUN | Approved eligible data and compatible locked model |

Issue #6 remains open. This preparation branch depends on the Phase II integration
PR #8. It will be reviewed separately and retargeted to main after that dependency
merges. Authorization to merge PR #4 does not authorize merging this research PR.

Locked `pip check`, compileall and Ruff also passed. These software results do
not constitute empirical evidence about the controlled cohort. The test fixture
and CLI manifest are retained under `results/phase_iii/independent_validation/`.
