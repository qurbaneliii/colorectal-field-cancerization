# Task B field-effect model report

Task B (healthy versus adjacent-normal) is the primary predictive endpoint.
Every outer split preserved donor/patient groups. Within each outer-training
partition, feature ranking and panel-size evaluation were repeated in grouped
inner folds; the smallest panel satisfying the configured macro-F1 tolerance,
balanced-accuracy floor, and calibration-loss limit was selected before the
outer test fold was evaluated once.

- Deterministic outer repeats: 3
- Outer folds per repeat: 5
- Mean nested-policy outer-fold macro-F1: 1.0000
- Final full-development panel: CLC, DYNC1H1, FOS, VIP, SNORA12
- Panel-size origin: GSE44076 grouped inner CV only
- Independent healthy-versus-adjacent validation: unavailable
- Intended use: retrospective biomarker discovery; not clinically ready

Patient/donor-bootstrap intervals of aggregated repeated OOF predictions:

| task                | estimand                                     | metric            |      mean |   median |   standard_deviation |   ci_lower |   ci_upper |   bootstrap_iterations | analysis_provenance   |
|:--------------------|:---------------------------------------------|:------------------|----------:|---------:|---------------------:|-----------:|-----------:|-----------------------:|:----------------------|
| task_b_field_effect | patient_bootstrap_of_aggregated_repeated_oof | roc_auc           | 1         | 1        |          3.70074e-17 |  1         |  1         |                   1000 | raw_cel_rma           |
| task_b_field_effect | patient_bootstrap_of_aggregated_repeated_oof | pr_auc            | 1         | 1        |          8.66834e-17 |  1         |  1         |                   1000 | raw_cel_rma           |
| task_b_field_effect | patient_bootstrap_of_aggregated_repeated_oof | f1_macro          | 1         | 1        |          0           |  1         |  1         |                   1000 | raw_cel_rma           |
| task_b_field_effect | patient_bootstrap_of_aggregated_repeated_oof | f1_weighted       | 1         | 1        |          0           |  1         |  1         |                   1000 | raw_cel_rma           |
| task_b_field_effect | patient_bootstrap_of_aggregated_repeated_oof | balanced_accuracy | 1         | 1        |          0           |  1         |  1         |                   1000 | raw_cel_rma           |
| task_b_field_effect | patient_bootstrap_of_aggregated_repeated_oof | brier_score       | 0.0267885 | 0.026674 |          0.00215564  |  0.0228925 |  0.0312873 |                   1000 | raw_cel_rma           |
| task_b_field_effect | patient_bootstrap_of_aggregated_repeated_oof | log_loss          | 0.164543  | 0.164377 |          0.0068046   |  0.152586  |  0.178205  |                   1000 | raw_cel_rma           |
