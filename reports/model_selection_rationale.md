# Model selection rationale

The prespecified priority order was leakage safety, external validity, balanced
performance, adjacent-normal recall, calibration, stability, compactness,
interpretability, and computational simplicity. Elastic Net was retained for
all three tasks. It was evaluated in the same repeated nested patient-group
cross-validation as the linear SVM and Random Forest; no point estimate alone
determined the decision.

| task                     | selected_model   |   elastic_net_mean_macro_f1 |   best_point_estimate | decision_rule                                                                                                                    |
|:-------------------------|:-----------------|----------------------------:|----------------------:|:---------------------------------------------------------------------------------------------------------------------------------|
| task_a_three_class       | elastic_net      |                    0.952751 |              0.962039 | Elastic Net retained for leakage safety, parsimony, stability, and interpretability; point estimates are not the sole criterion. |
| task_b_field_effect      | elastic_net      |                    0.992605 |              1        | Elastic Net retained for leakage safety, parsimony, stability, and interpretability; point estimates are not the sole criterion. |
| task_c_tumor_vs_adjacent | elastic_net      |                    0.987218 |              0.987218 | Elastic Net retained for leakage safety, parsimony, stability, and interpretability; point estimates are not the sole criterion. |

Task C locked signature size: 19 genes. The model artifact
was serialized and reloaded successfully. External-cohort labels were not used
in this fit or any hyperparameter decision.
