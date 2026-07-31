# Model selection rationale

Elastic Net is the prespecified primary model because the study is a sparse,
interpretable biomarker-discovery analysis. Linear SVM and Random Forest are
benchmarks, not candidates selected by their external-cohort behavior. All
models used the same repeated nested patient/donor-group folds, and all learned
preprocessing and calibration occurred inside training data.

| task                     | selected_model   |   elastic_net_mean_macro_f1 |   best_comparator_point_estimate | decision_rule                                                                                         | analysis_provenance   |
|:-------------------------|:-----------------|----------------------------:|---------------------------------:|:------------------------------------------------------------------------------------------------------|:----------------------|
| task_a_three_class       | elastic_net      |                    0.941281 |                         0.961987 | Elastic Net was prespecified as the primary scientific model; comparators are sensitivity benchmarks. | raw_cel_rma           |
| task_b_field_effect      | elastic_net      |                    0.992628 |                         0.997558 | Elastic Net was prespecified as the primary scientific model; comparators are sensitivity benchmarks. | raw_cel_rma           |
| task_c_tumor_vs_adjacent | elastic_net      |                    0.984713 |                         0.989811 | Elastic Net was prespecified as the primary scientific model; comparators are sensitivity benchmarks. | raw_cel_rma           |

Expression provenance: `raw_cel_rma`. Compact-panel performance and the final
locked Task C artifact are produced separately by
`scripts/run_compact_panel_analysis.py`.
