# FINAL PHASE II REPORT

[Executive summary](01_EXECUTIVE_SUMMARY.md) · [Code review](02_INDEPENDENT_STYLE_CODE_REVIEW.md) · [Merge gate](13_PR4_MERGE_READINESS.md)

The independent check reproduced primary Task C ROC-AUC 0.982894, macro-F1 0.828155, sensitivity 1 and specificity 0.557692 from 233 arrays across 190 patients. The confirmed SVM calibration defect was internal to the comparator and has regression coverage. Task B case/acquisition rank deficiency leaves the cancer-specific interpretation unresolved.

[Reproducibility](11_REPRODUCIBILITY_ACCEPTANCE.md) · [Dataset registry](07_EXTERNAL_DATASET_REGISTRY.md) · [External protocol](09_EXTERNAL_VALIDATION_PROTOCOL.md) · [Publication](12_PUBLICATION_READINESS.md) · [Blockers](14_UNRESOLVED_RESEARCH_BLOCKERS.md)

**READY FOR OWNER-APPROVED MERGE.** Both CI checks and the full hosted raw-CEL run passed for the corrected integration head. The frozen external prototype passed synthetic tests and remains separate. [Exact verification and recommendation](13_PR4_MERGE_READINESS.md)
