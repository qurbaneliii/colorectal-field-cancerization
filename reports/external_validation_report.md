# External validation report

The GSE41258 endpoint was primary tumor versus normal colon in patients with
colorectal neoplasia. Task B requires a separate cohort with cancer-free controls
and cancer-associated normal mucosa collected under comparable protocols.

The full GSE44076 signature contained FOXQ1, CEMIP, ETV4.
The cross-platform model was refitted on GSE44076 with
CEMIP, ETV4. The threshold was selected on
grouped GSE44076 out-of-fold scores with the full-development panel and
hyperparameters held fixed. GSE41258 outcomes were used for evaluation.
The cohort has been examined in earlier project analyses.

- Locked threshold 0.41
- Default threshold sensitivity 0.50
- Canonical set of 233 arrays from 190 patients
- ROC-AUC 0.9829 and average precision 0.9935
- Balanced accuracy 0.7788 and macro-F1 0.8282
- Sensitivity 1.0000 and specificity 0.5577
- Brier score 0.0533 and log loss 0.1950
- Calibration intercept -2.0661 and slope 2.6198

Patient-cluster bootstrap intervals and the sensitivity estimands are in the
external validation CSVs. Calibration coefficients describe score transport in
this cohort. They were not applied to the model or the locked predictions.
