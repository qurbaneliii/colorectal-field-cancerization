# Task B confounding sensitivity

Task B is the primary healthy-versus-adjacent endpoint. All molecular
sensitivity analyses use the five-gene panel locked by grouped inner CV. The
covariate residualizers are fitted on each outer-training fold and applied to
the held-out fold, preventing outcome or test-row leakage. The demographic-only
model uses age, sex, and left/right location. The matching analysis uses
one-to-one, without-replacement same-sex matching within a
10-year age caliper.

- Locked panel: CLC, DYNC1H1, FOS, VIP, SNORA12
- Primary fixed-panel aggregated OOF macro-F1: 1.0000
- Matched subset: 88 samples (44 per class)
- Custom-QC exclusion sensitivity: 147 samples; reselected panel
  (5 genes): CLC, DYNC1H1, FOS, VIP, APOLD1
- Independent Task B validation: unavailable
- Interpretation: internally validated association model; not clinically ready

## Aggregated repeated OOF estimates

| scenario                                 |   n_samples |   n_patients | panel_genes_locked   |   f1_macro |   f1_weighted |   balanced_accuracy |   log_loss |   precision_adjacent_normal |   recall_adjacent_normal |   f1_adjacent_normal |   precision_healthy |   recall_healthy |   f1_healthy |   roc_auc |   pr_auc |   brier_score |   sensitivity |   specificity |   negative_predictive_value |   positive_predictive_value |   calibration_intercept |   calibration_slope |
|:-----------------------------------------|------------:|-------------:|:---------------------|-----------:|--------------:|--------------------:|-----------:|----------------------------:|-------------------------:|---------------------:|--------------------:|-----------------:|-------------:|----------:|---------:|--------------:|--------------:|--------------:|----------------------------:|----------------------------:|------------------------:|--------------------:|
| age_sex_location_residualized            |         148 |          148 | True                 |   0.977234 |      0.979679 |            0.974898 |   0.206645 |                    0.979798 |                 0.989796 |             0.984772 |            0.979592 |             0.96 |     0.969697 |  0.993469 | 0.990688 |     0.0461847 |          0.96 |      0.989796 |                    0.979798 |                    0.979592 |              -0.182279  |             3.41526 |
| age_sex_matched_fixed_panel              |          88 |           88 | True                 |   1        |      1        |            1        |   0.202858 |                    1        |                 1        |             1        |            1        |             1    |     1        |  1        | 1        |     0.0373315 |          1    |      1        |                    1        |                    1        |              -0.0502956 |            10.2422  |
| borderline_sample_excluded_fixed_panel   |         147 |          147 | True                 |   1        |      1        |            1        |   0.131913 |                    1        |                 1        |             1        |            1        |             1    |     1        |  1        | 1        |     0.0184442 |          1    |      1        |                    1        |                    1        |              -0.45916   |             7.71161 |
| composition_and_demographic_residualized |         148 |          148 | True                 |   0.718149 |      0.736678 |            0.742041 |   0.566695 |                    0.8625   |                 0.704082 |             0.775281 |            0.573529 |             0.78 |     0.661017 |  0.785918 | 0.652525 |     0.192512  |          0.78 |      0.704082 |                    0.8625   |                    0.573529 |              -0.685796  |             1.46224 |
| demographic_only                         |         148 |          148 | False                |   0.626036 |      0.663733 |            0.627347 |   0.642498 |                    0.75     |                 0.734694 |             0.742268 |            0.5      |             0.52 |     0.509804 |  0.682041 | 0.578544 |     0.225326  |          0.52 |      0.734694 |                    0.75     |                    0.5      |              -0.632357  |             1.78235 |
| primary_fixed_panel                      |         148 |          148 | True                 |   1        |      1        |            1        |   0.131182 |                    1        |                 1        |             1        |            1        |             1    |     1        |  1        | 1        |     0.0182343 |          1    |      1        |                    1        |                    1        |              -0.404288  |             7.6561  |

Gene coefficients and selection frequencies for the confounding scenarios are
reported separately. Those scenarios keep the primary panel locked, so
coefficient shrinkage to zero is their predefined gene-level sensitivity
measure. The QC-only analysis additionally repeats full-development grouped
inner-CV panel selection after omitting the borderline array; it is reported as
a sensitivity panel and never replaces the primary model.
