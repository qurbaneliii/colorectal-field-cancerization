# External validation report

This is a secondary cross-platform tumor-normal evaluation, not independent
validation of the healthy-versus-adjacent field-effect task. The exact
3-gene Task C model could not be deployed because
FOXQ1 was
absent on GPL96. The feature specification, Elastic Net family,
hyperparameters, within-sample rank representation, and threshold were locked
using GSE44076. A separate transport model was refitted on GSE44076 using the
2 common signature genes and evaluated once on GSE41258.
External labels were never used for fitting, selection, transformation, or
threshold choice.

- External labels retained verbatim: `normal_colon`, `primary_tumor`
- Locked threshold: 0.41 (GSE44076 grouped repeated OOF)
- Default threshold sensitivity: 0.50
- Primary set: 233 arrays from 190 patients
- Primary ROC-AUC: 0.9829
- Primary PR-AUC: 0.9935
- Primary balanced accuracy: 0.7788
- Primary macro-F1: 0.8282
- Primary NPV: 1.0000
- Primary Brier score: 0.0533
- Primary log loss: 0.1950

Patient-cluster bootstrap intervals, the all-array analysis, the optional
one-array-per-patient analysis, and the failed or degraded training-z-score
transport are reported as distinct estimands. High ranking discrimination does
not by itself establish threshold or probability calibration transport. This
retrospective biomarker-discovery model is not clinically ready.
