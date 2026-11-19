# Measured platform transport

The shared annotation universe contains 12,039 genes. CEMIP and ETV4 are refit on development ranks; FOXQ1 is absent externally. The experimental frozen reference is learned on 196 development arrays.

| method    |   threshold |   f1_macro |   balanced_accuracy |   roc_auc |   pr_auc |   sensitivity |   specificity |   brier_score |   log_loss |   tn |   fp |   fn |   tp |   calibration_intercept |   calibration_slope | calibration_converged   |   mean_probability |   observed_fraction |
|:----------|------------:|-----------:|--------------------:|----------:|---------:|--------------:|--------------:|--------------:|-----------:|-----:|-----:|-----:|-----:|------------------------:|--------------------:|:------------------------|-------------------:|--------------------:|
| frozen    |    0.410000 |   0.949071 |            0.936783 |  0.982894 | 0.993594 |      0.988950 |      0.884615 |      0.027588 |   0.121138 |   46 |    6 |    2 |  179 |               -0.270738 |            1.046918 | True                    |           0.786387 |            0.776824 |
| joint_rma |    0.410000 |   0.854774 |            0.807692 |  0.982682 | 0.993460 |      1.000000 |      0.615385 |      0.049353 |   0.188770 |   32 |   20 |    0 |  181 |               -2.114153 |            1.046613 | True                    |           0.851285 |            0.776824 |

Both arms above use a new prespecified L2 refit and retain 0.41 as a historical reference cutoff. The primary archived Elastic Net transport result remains ROC-AUC 0.982894, macro-F1 0.828155, sensitivity 1.0 and specificity 0.557692. Deposited-expression specificity 0.942 remains a separate sensitivity. Differences between the new joint-RMA arm and that historical model include classifier policy and cannot be attributed to normalization.

No external outcome selected the method, parameters, panel or cutoff. GSE41258 is previously exposed retrospective evidence. Rank transport does not establish calibrated probabilities or clinical validity.
