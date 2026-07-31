# External validation report

The prespecified primary transfer representation was within-sample percentile
rank because it is computed independently within each sample, does not require
external labels, reduces reliance on platform-specific intensity scale, and
preserves relative gene ordering. Training-derived z-scores are a sensitivity
analysis. GSE41258 labels were not used for feature, panel, model,
transformation, threshold, or hyperparameter selection.

- Eligible arrays: 233 (181 primary tumor; 52 normal colon)
- Eligible unique patients: 190
- Primary deterministic one-array-per-patient subset: 190
- Original Task C signature: 5 genes
- Cross-platform common signature: 4 genes
- Excluded signature genes: FOXQ1
- Locked threshold: 0.5

Primary point estimates: ROC-AUC 0.9969, PR-AUC
0.9994, balanced accuracy
0.9167, macro F1
0.7919, Brier score
0.1249, and log loss
0.4118. Patient-cluster bootstrap intervals and the
all-array sensitivity analysis are saved in `results/metrics`.

This validates only the tumor-versus-normal-colon component. GSE41258 is not an
independent validation cohort for cancer-free healthy versus tumor-adjacent
field cancerization.
