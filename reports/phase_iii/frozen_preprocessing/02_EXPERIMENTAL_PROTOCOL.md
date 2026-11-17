# Locked experimental protocol

The machine-readable protocol was committed before model results were inspected
at `007d2c9`. It specifies 4,097 quantiles, seed 44076, historical fixed Task B/C
panels, balanced L2 logistic regression with C=1, a 0.5 internal cutoff and the
historical 0.41 external reference cutoff. There is no parameter search or
validation-label cutoff selection. [Protocol](../../../experiments/frozen_reference/protocol.json)

All three archived repeats and five outer folds per task are reused. For every
partition, the quantile target and scaler are fit only on its training arrays.
The model uses a fixed panel and fixed hyperparameters, so there are no inner
validation folds or reused inner references. This is a conditional panel
experiment. Historical panel discovery used the development cohort, so these
estimates are not unbiased de novo feature-selection performance.

The same classifier is fit on archived joint-RMA values for a paired baseline.
This baseline is distinct from the published nested Elastic Net pipeline. Only
the frozen arm processes held-out CELs without the rest of the evaluation cohort.
Shared manufacturer annotation is fixed independently of outcomes.

Invariance tolerances are 1e-12 absolute and relative. A raw-PM check compares
independent background correction with multicolumn processing for first/last
members of each tissue class. Every internal held-out panel is transformed alone
and in a batch using the same immutable reference; probabilities and decisions
are checked. Paired 1,000-draw donor-cluster bootstrap intervals use repeat-averaged
OOF scores and remain conditional on the fitted models and historical panel.

External transport refits the two covered genes CEMIP/ETV4 on all eligible
development patients, using within-sample percentile ranks over a fixed shared
gene universe. FOXQ1 remains unavailable on the external platform. Scores are
written before outcome metric computation. GSE41258 is previously exposed
retrospective evidence. The 0.41 cutoff is a historical reference, not a newly
validated operating point for this experimental L2 model.
