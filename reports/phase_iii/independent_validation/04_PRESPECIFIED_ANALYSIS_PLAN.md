# Proposed prespecified independent analysis

Status: a reviewable protocol before access, not a completed preregistration or
an empirical validation. Freeze a version and institutional approval before
restricted outcomes or expression are inspected. Hypothesis: a development-locked
Task B score discriminates healthy from cancer-associated normal mucosa under
comparable acquisition and site conditions. Causality is not the estimand.

Use the primary eligibility in report 02. The primary metric is donor-level ROC
AUC with 2,000 stratified donor bootstrap draws (seed 44076). Secondary metrics are
PR AUC, macro-F1, balanced accuracy, sensitivity, specificity, Brier score and
log loss at the development-locked cutoff. Evaluate calibration intercept/slope
descriptively; do not recalibrate the primary test on validation outcomes. Record
nonconvergence and wide intervals. Intervals are conditional on the locked model.

Freeze gene identifiers, assay coverage, normalization, rank universe, classifier,
class order and threshold before scoring. CLC, DYNC1H1, FOS, VIP and SNORA12 are the
historical Task B panel. Missing genes or incompatible quantification block the
original-panel transport claim. No external-label gene selection, threshold
search or model comparison selects a winning primary model.

Age and sex are required for the primary eligible set. Record all exclusions and
missingness by cohort/site. Secondary logistic association models include fixed
score, age and sex, with center/batch only if estimable. Check matrix rank, class
overlap and separation first. Stage, distance and left/right-site contrasts are
secondary descriptive analyses. Use donor grouping for repeated samples. Apply
Benjamini-Hochberg within the prespecified secondary association family; a sole
primary AUC estimate has no multiple-model selection claim. No prospective power
claim is justified until eligible donor counts are known.

Inspect batch and covariate structure without outcome-driven correction. No
joint development/test ComBat. Technical replicate agreement and label-free
sample-alone invariance are justified technical controls. Shuffled labels can
test software behavior but cannot resolve acquisition confounding. Retain a log
of every evaluation and protocol amendment. Poor discrimination, missing
coverage or infeasibility must be reported without replacing the primary test.
