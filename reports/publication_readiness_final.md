# Publication-readiness final report

## 1. Executive summary

Final readiness decision: **PARTIAL**. Scientific claims remain restricted
to retrospective biomarker discovery; no clinical-readiness claim is made.

## 2. Repository changes

Raw platform assertions, provenance separation, limma, enrichment,
configuration-driven grouped nested CV, group-calibrated SVM, fold-internal
compact panels, patient-balanced external validation, CI, locks, tests, and
manuscript generators were added or corrected.

## 3. Raw-data audit

GSE44076 contains 246 arrays and GSE41258 contains 390 arrays. Detailed
archive, mapping, and checksum evidence is linked below.

## 4. RMA and QC

Raw-CEL completion status: True. No PCA-only exclusion is permitted.

## 5. Sample exclusions

Deterministic metadata exclusions and technical candidates are recorded in
`data/metadata/sample_exclusion_log.csv` and the metadata audit.

## 6. Final differential expression

| comparison               | direction   |   tested_genes |   significant_genes |
|:-------------------------|:------------|---------------:|--------------------:|
| adjacent_vs_healthy      | down        |          11934 |                 659 |
| adjacent_vs_healthy      | up          |           6556 |                1003 |
| tumor_vs_adjacent_paired | down        |           9393 |                2025 |
| tumor_vs_adjacent_paired | up          |           9097 |                2477 |
| tumor_vs_healthy         | down        |          10401 |                2239 |
| tumor_vs_healthy         | up          |           8089 |                2764 |

## 7. Field-cancerization results

Candidate tables prioritize raw-CEL limma, sample consistency, sensitivity
agreement, and uncertainty; three group means are not treated as progression.

## 8. Enrichment

Direction- and contrast-specific GO, Reactome, and ranked results use the tested
gene universe. Enrichment is associative, not causal.

## 9. Model methodology

Elastic Net was prespecified; SVM and Random Forest are comparators.

## 10. Leakage controls

Outer, inner, calibration, feature-selection, compact-panel, and threshold
operations preserve patient/donor boundaries and training-only fitting.

## 11. Internal validation

| task                     |   mean_macro_f1 |   median_macro_f1 |   sd_macro_f1 |   fold_minimum |   fold_maximum |
|:-------------------------|----------------:|------------------:|--------------:|---------------:|---------------:|
| task_a_three_class       |        0.941281 |          0.951178 |     0.02413   |       0.89418  |       0.975118 |
| task_b_field_effect      |        0.992628 |          1        |     0.0152627 |       0.962677 |       1        |
| task_c_tumor_vs_adjacent |        0.984713 |          1        |     0.0188782 |       0.947222 |       1        |

## 12. Feature stability

Selection counts use unique `(repeat, outer_fold)` keys and report sign and
coefficient variability.

## 13. Compact-panel selection

Task C signature (5 genes): FOXQ1, CEMIP, ETV4, GTF2IRD1, PACC1.

## 14. Calibration

Brier score, log loss, calibration curves, and binary calibration intercept and
slope are generated where estimable.

## 15. Permutation tests

Task A is explicitly conditional; Task B is group-level; Task C is paired
within patient. Resolution is controlled by the configured iteration count.

## 16. External validation

| representation                | evaluation_set    | is_primary   |   arrays |   unique_patients | analysis_provenance   |   f1_macro |   f1_weighted |   balanced_accuracy |   log_loss |   precision_adjacent_normal |   recall_adjacent_normal |   f1_adjacent_normal |   precision_tumor |   recall_tumor |   f1_tumor |   roc_auc |   pr_auc |   brier_score |   sensitivity |   specificity |   negative_predictive_value |   positive_predictive_value |   calibration_intercept |   calibration_slope |
|:------------------------------|:------------------|:-------------|---------:|------------------:|:----------------------|-----------:|--------------:|--------------------:|-----------:|----------------------------:|-------------------------:|---------------------:|------------------:|---------------:|-----------:|----------:|---------:|--------------:|--------------:|--------------:|----------------------------:|----------------------------:|------------------------:|--------------------:|
| within_sample_percentile_rank | canonical_patient | True         |      190 |               190 | raw_cel_rma           |   0.791895 |      0.874549 |            0.916667 |   0.411845 |                    0.509091 |                        1 |             0.674699 |                 1 |       0.833333 |   0.909091 |  0.996914 | 0.999447 |        0.1249 |      0.833333 |             1 |                    0.509091 |                           1 |                     nan |                 nan |

## 17. Sensitivity analyses

Raw-versus-deposited, compact size, representation, threshold, and
canonical-patient versus all-array results are preserved separately.

## 18. Reproducibility

Python and R locks, Docker restore, normal CI, manual full-data CI, Make targets,
test evidence, and artifact checksums are included.

## 19. Remaining limitations

See `manuscript/limitations.md`; independent Task B field-effect validation and
prospective biological/clinical validation remain required.

## 20. Final readiness decision

**PARTIAL**

## Final evidence table

| Requirement              | Status   | Evidence file                                            | Verification command                                                                     | Key result                                         | Remaining action                                   |
|:-------------------------|:---------|:---------------------------------------------------------|:-----------------------------------------------------------------------------------------|:---------------------------------------------------|:---------------------------------------------------|
| A - Data integrity       | PASS     | reports/cel_archive_audit.md                             | Rscript R/01_extract_and_audit_cel.R                                                     | Required evidence exists and was inspected         | None                                               |
| B - Raw preprocessing    | PASS     | data/metadata/GSE44076_platform_validation.csv           | Rscript R/02_preprocess_gse44076.R; Rscript R/03_preprocess_gse41258.R                   | Required evidence exists and was inspected         | None                                               |
| C - Biological inference | PASS     | results/tables/supplementary_full_raw_cel_de_results.csv | Rscript R/04_differential_expression.R; Rscript R/05_functional_enrichment.R             | Required evidence exists and was inspected         | None                                               |
| D - Predictive validity  | PASS     | results\metrics\compact_panel_fold_metrics.csv           | python scripts/run_compact_panel_analysis.py --provenance raw_cel_rma                    | Required evidence exists and was inspected         | None                                               |
| E - External validity    | PASS     | results\metrics\external_validation_metrics.csv          | python scripts/run_external_validation.py --provenance raw_cel_rma                       | Required evidence exists and was inspected         | None                                               |
| F - Reproducibility      | PARTIAL  | reports/docker_verification.md                           | python -m pip check; Rscript R/verify_environment.R; python -m pytest -q; docker build . | One or more required evidence artifacts are absent | Complete and rerun the listed verification command |
| G - Publication outputs  | PASS     | manuscript/manuscript.md                                 | python scripts/build_manuscript_outputs.py                                               | Required evidence exists and was inspected         | None                                               |
