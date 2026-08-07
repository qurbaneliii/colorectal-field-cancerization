# Publication-readiness acceptance report

This report records acceptance gates, not scientific narrative. PASS requires an
executed, inspected, non-empty artifact; code presence alone is insufficient.

| Gate | Status | Executed evidence | Remaining action |
|---|---|---|---|
| A. Data integrity | PASS | Metadata/CEL audit, 246 GSE44076 arrays, 390 GSE41258 arrays, 98 complete pairs, canonical external patient-tissue table | None |
| B. Raw-array QC | PASS | All-array densities, custom diagnostics, two arrayQualityMetrics reports, GPL96 NUSE, exclusion sensitivity | None |
| C. Differential expression | PASS | Full-rank U0/U1/U2 and patient-fixed-effect limma outputs with BH FDR and provenance | None |
| D. Confounding and biology | PASS | Demographic, subgroup, MCP-counter, stress, QC, preprocessing, trajectory, and tier-specific enrichment sensitivities | Independent cell-resolved confirmation remains a scientific limitation |
| E. Modeling and external boundaries | PASS | Grouped nested selection; Task B primary; Task C secondary; locked GSE44076 threshold; exact model versus transport-refit distinction | Independent Task B cohort remains unavailable |
| F. Uncertainty | PASS | Separate fold and aggregated-OOF estimands, 1,000 group bootstraps, task-specific 1,000-permutation tests, calibration | Prospective uncertainty remains unavailable |
| G. Reproducibility | PARTIAL | Locked Python/R files, 40 passing local tests, R verification, and a passing baseline-branch CI run | Current-revision CI awaits an authenticated push; Docker is host-blocked until the Windows Linux-container engine/WSL service is available |
| H. Manuscript | PASS | Complete sections, verified references, no citation placeholders, expanded results/discussion, full legends, ten tables, vector/raster figures | Supply author, affiliation, and funding metadata |

Overall decision: **PARTIAL publication readiness**. The scientific and reporting
gates are complete, but current-revision remote CI and independent container
execution remain blocked at the host. This does not invalidate the executed
analyses; it prevents a full reproducibility PASS.
