# Predictive validation audit

## Task B | Cohort association and sensitivity

The original nested compact-panel policy has 15 outer folds over three repeats
and reports mean macro-F1, balanced accuracy and ROC-AUC of 1.000. Grouped
training/test separation, fold-local feature filtering and inner tuning are
visible in `src/modeling/nested_cv.py` and `scripts/run_compact_panel_analysis.py`.
The inputs already contain information from cohort-wide quantile normalization
and multiarray summarization. The study estimates discrimination of the
retrospectively processed specimens. Deployable single-sample performance is
unmeasured.

The archived Task B confounding analysis used the full-development panel and
its hyperparameters, while residualization and composition PCA were fold-local.
Reported macro-F1 was 1.000 for the fixed panel, 0.626 for demographics, 0.977
after demographic residualization and 0.850 after demographics and composition.
The age/sex matched subset retains the original acquisition difference.

The new GSE44076 sensitivity used the exact GEO-deposited matrix, 18,490
mapped gene summaries and repeated grouped CV. Panels of up to five genes
were reselected inside each training fold with fixed Elastic Net parameters.
The sample-aggregated repeated out-of-fold results were:

| Scenario | Macro-F1 | ROC-AUC | Interpretation |
|---|---:|---:|---|
| Demographics alone | 0.629 | 0.690 | Measured case mix carries some signal |
| Marker-expression composition proxies alone | 0.941 | 0.987 | A substantial signal is shared with bulk composition proxies |
| Curated stress panel alone | 0.992 | 0.998 | Preanalytical sensitivity is plausible |
| Fold-selected expression | 1.000 | 1.000 | Separation persists with held-out donors |
| Reselect after removing known markers and stress genes | 1.000 | 1.000 | The signal extends beyond the two curated lists |
| Reselect after removing those genes and residualizing | 0.843 | 0.932 | Separation remains after this adjustment |

The last estimate's 95% donor bootstrap interval was 0.777–0.901. The paired
macro-F1 change relative to fold-selected expression was −0.157
(95% interval −0.223 to −0.099). This is a conditional resampling interval of
stored out-of-fold predictions. It does not include uncertainty from training
new models or establish a causal decomposition of the signal. Removing marker
genes and fitting proxy residualization are distinct interventions. Correlated
genes and unmeasured acquisition effects can persist.

## Task C | Selection and threshold estimands

The old threshold lock used a full-development signature and hyperparameters.
Its grouped refits avoid direct test-patient fitting, but selection used the
full development cohort. The reported 0.995 development balanced accuracy was
also calculated where the threshold was selected.

`scripts/run_development_audit.py` evaluated a separate three-gene rank policy
using GSE44076 deposited expression. Each threshold-fit cross-fit had its own
grouped tuning and training-only gene selection. Outer held-out patients
contributed to evaluation after all choices were made. Platform overlap used
the published annotation maps only. GSE41258 expression and outcomes were
absent from this experiment.

Thresholds ranged from 0.39 to 0.50 across 15 outer training sets. Both the
inner-selected and fixed 0.50 operating policies had the same outer-fold
mean macro-F1 (0.9898), balanced accuracy (0.9898), ROC-AUC (1.0000),
Brier (0.0187) and log loss (0.1173). There was no observed operating gain
under this policy on these samples. The refits and preprocessing differ
from the historical two-gene raw-CEL transport model. The new analysis
does not supersede its external threshold.

## Permutation and clinical boundaries

The archived permutation nulls preserve patient groups but reuse
hyperparameters and splits selected under the observed outcomes. They refer
to a conditional diagnostic, with Task A fixing healthy labels. Task B's
class–protocol alias violates exchangeability for a cancer-specific null.
A complete procedure would regenerate all data-dependent splits and
selection under each permutation and define an exchangeable design first.
Retrospective discrimination does not establish clinical utility.
See [Varma and Simon](https://doi.org/10.1186/1471-2105-7-91) and
[TRIPOD+AI](https://doi.org/10.1136/bmj-2023-078378).
