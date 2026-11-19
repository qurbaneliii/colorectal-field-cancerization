# Measured model comparison

These are conditional historical-panel L2 logistic-regression comparisons. The same model settings are used in both arms. They are distinct from the published nested Elastic Net estimates.

Mean across the 15 paired held-out folds per task

| task                     | method    |   f1_macro |   balanced_accuracy |   roc_auc |   pr_auc |   sensitivity |   specificity |   brier_score |   log_loss |
|:-------------------------|:----------|-----------:|--------------------:|----------:|---------:|--------------:|--------------:|--------------:|-----------:|
| task_b_field_effect      | frozen    |   1.000000 |            1.000000 |  1.000000 | 1.000000 |      1.000000 |      1.000000 |      0.001027 |   0.016775 |
| task_b_field_effect      | joint_rma |   1.000000 |            1.000000 |  1.000000 | 1.000000 |      1.000000 |      1.000000 |      0.001103 |   0.016696 |
| task_c_tumor_vs_adjacent | frozen    |   0.994909 |            0.994912 |  0.999315 | 0.999363 |      0.989825 |      1.000000 |      0.006184 |   0.033589 |
| task_c_tumor_vs_adjacent | joint_rma |   0.994909 |            0.994912 |  0.999500 | 0.999538 |      0.989825 |      1.000000 |      0.006200 |   0.034022 |

The paired donor-bootstrap intervals in `paired_conditional_intervals.csv` use repeat-averaged predictions, not the mean-fold estimand above. They do not cover feature discovery or cross-study uncertainty.

Repeat-averaged metrics corresponding to that interval estimand

| task                     | method    |   f1_macro |   balanced_accuracy |   roc_auc |
|:-------------------------|:----------|-----------:|--------------------:|----------:|
| task_b_field_effect      | frozen    |   1.000000 |            1.000000 |  1.000000 |
| task_b_field_effect      | joint_rma |   1.000000 |            1.000000 |  1.000000 |
| task_c_tumor_vs_adjacent | frozen    |   0.994898 |            0.994898 |  0.999688 |
| task_c_tumor_vs_adjacent | joint_rma |   0.994898 |            0.994898 |  0.999688 |

Calibration fits and convergence flags are retained for each fold; separation can prevent reliable slope estimation. Brier score and log loss remain reported.
