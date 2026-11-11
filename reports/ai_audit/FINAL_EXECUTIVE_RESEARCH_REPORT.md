# Independent colorectal field-cancerization audit

**Primary conclusion.** GSE44076 contains a strong transcriptomic association
between the two normal-mucosa cohorts. Controls were biopsied at colonoscopy.
Cancer-associated mucosa came from surgical margins at least 10 cm from the
tumor. This study design cannot establish how much of the signal is
cancer-specific. [Original GSE44076 study](https://doi.org/10.1186/1476-4598-13-46)

| Question | Audited evidence |
|---|---|
| Do the cohorts and pairings match GEO? | Yes. GSE44076 has 246 arrays from 148 groups, including 98 tumor/mucosa pairs. GSE41258 has 233 eligible arrays from 190 neoplasia patients, including 43 pairs. |
| Does Task B persist with fold-local selection and confounding challenges? | On the deposited matrix, fold-selected macro-F1 was 1.000, stress-only 0.992 and composition-proxy-only 0.941. After marker/stress exclusion and fold-local residualization, macro-F1 was 0.843 (95% conditional donor bootstrap interval 0.777–0.901). |
| Does the Task B score travel? | The separate four-gene GSE4107 refit achieved ROC-AUC 0.992 in 22 donors and macro-F1 0.313 at the fixed threshold. Both cohorts have the same acquisition alias. |
| What does GSE41258 support? | The two-gene Task C refit achieved replayed ROC-AUC 0.983, macro-F1 0.828, sensitivity 1.000 and specificity 0.558 at 0.41. This is tumor versus normal colon in a neoplasia cohort. |
| Does a more deeply nested Task C threshold improve development performance? | A separate deposited rank policy selected cutoffs from 0.39 to 0.50 across 15 folds. Mean outer macro-F1 was 0.990 under both the inner-selected and default policies. |
| Do the limma effects agree across expression sources? | The four deposited-expression refits correlated above 0.99998 with the archived raw-CEL effect estimates on the same samples. The acquisition confounding applies to both routes. |
| Is there a better Task B cohort? | Ikuta et al. sampled healthy and cancer-associated mucosa endoscopically. JGAD000785 requires controlled access. |

The branch archives the older 19/14-gene experiment, verifies the active
external result against saved predictions, restores GEO scan-date metadata,
binds the threshold lock to its input bytes and corrects probability column
mapping. Current manuscript numerical claims were preserved and their
methods and scope clarified.

**Reproduction status.** Both current GEO archives passed complete integrity
and metadata checks. Both differ from the historical TAR hashes. The original
TARs or historical CEL digests are required to attribute those differences.
Fresh RMA reproduced all eight archived limma outputs to numerical precision.
The grouped nested CV reproduced all 135 fold macro-F1 values and
hyperparameter selections. A fresh Task C refit conditional on the archived
panel reproduced the 0.41 lock and the primary external operating point,
calibration and bootstrap interval endpoints. The continuation records input
hashes, stage-level execution and remaining gates.
[Raw-CEL reproduction continuation](15_raw_reproduction_continuation.md)

## Audit package

- [Repository forensic audit](01_repository_forensic_audit.md)
- [Baseline reproduction](02_baseline_reproduction.md)
- [Data quality](03_data_quality_report.md)
- [Statistical genomics](04_statistical_genomics_audit.md)
- [Signature robustness](05_signature_robustness_report.md)
- [Machine-learning validation](06_ml_validation_audit.md)
- [Cross-platform transport](07_cross_platform_transport_report.md)
- [External cohort research](08_external_dataset_research.md)
- [Artifact consistency](09_artifact_consistency_report.md)
- [Experiment register](10_experiment_register.md)
- [Implementation changelog](11_implementation_changelog.md)
- [Before and after evaluation](12_before_after_evaluation.md)
- [Remaining limitations](13_remaining_limitations.md)
- [Future research roadmap](14_future_research_roadmap.md)
- [Raw-CEL reproduction continuation](15_raw_reproduction_continuation.md)
