# Introduction

Field cancerization describes molecular alteration extending beyond a visible
neoplasm into apparently non-neoplastic tissue. Genomic and transcriptomic
differences have been demonstrated in colorectal mucosa surrounding tumors
[1], and GSE44076 was specifically generated to compare healthy-donor mucosa,
tumor-adjacent mucosa, and matched colorectal tumors [2]. Tumor-adjacent normal
tissue is therefore a biologically informative compartment rather than a
universally exchangeable healthy control [3]. Pan-cancer analyses likewise show
that adjacent tissue can contain inflammation, stromal response, and other
signals distinct from both healthy tissue and tumor [4].

Three questions must remain separate. The primary biological and predictive
question is whether healthy-donor mucosa differs from adjacent-normal mucosa
(Task B). A supporting three-class analysis asks whether healthy, adjacent, and
tumor tissue states can be discriminated (Task A), without interpreting group
means as a longitudinal sequence. A secondary question asks whether a
GSE44076-derived tumor-versus-adjacent classifier can transfer across Affymetrix
platforms to GSE41258 tumor and normal-colon specimens (Task C). GSE41258 cannot
independently validate healthy-versus-adjacent field cancerization because it
does not contain the required healthy and adjacent-normal endpoint definition.

Microarray biomarker discovery creates several opportunities for optimism:
normalization may be dataset-specific, multiple samples from one patient may
cross validation folds, feature and panel selection may use test outcomes, and
thresholds may be tuned on external labels. Nested cross-validation is required
when tuning is part of the model-building procedure [10], and grouped resampling
is necessary when observations share a patient or donor [11]. This study was
therefore designed around raw-data provenance, prespecified covariate and
composition sensitivity, training-only transformations, grouped nested
validation, task-specific permutation tests, and explicit separation of the
primary model from the cross-platform transport refit.

We aimed to identify candidate field-cancerization genes with graded evidence;
quantify demographic, tissue-composition, stress, and QC sensitivity; derive an
interpretable internally validated Task B classifier; and evaluate a distinct
Task C tumor-normal signal across platforms while preserving honest boundaries
on calibration, transport, and clinical readiness.
