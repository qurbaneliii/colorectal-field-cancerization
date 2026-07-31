# Methods

GSE44076 (GPL13667) was the sole development cohort (50 healthy donors, 98
adjacent-normal samples, and 98 paired tumors). GSE41258 (GPL96) was restricted
to author-included Primary Tumor and Normal Colon arrays after deterministic
removal of metastases, non-colon tissues, polyps, cell lines, ambiguous records,
and `_ez`/`rehyb` technical-replicate candidates.

Raw CEL files were processed independently by platform using oligo RMA for
HG-U219 and affy RMA for HG-U133A. CEL chip/CDF compatibility, GEO identifiers,
sample counts, finite log2-scale output, and probe-set counts were asserted.
Probes required unambiguous symbol and Entrez mappings; multiple probes per
gene were aggregated by the per-sample median. Raw and normalized
distributions, RLE, MA diagnostics, PCA, correlation, clustering, robust
sample-median/IQR diagnostics, and arrayQualityMetrics were generated. PCA
separation alone was not an exclusion criterion.

Primary differential expression used limma on raw-CEL RMA expression. Unpaired
designs tested adjacent-normal versus healthy and tumor versus healthy. A
patient fixed-effect design tested paired tumor versus adjacent-normal. All
coefficient names and design rank were checked programmatically, and BH FDR was
applied. Field candidates required FDR < 0.05, |log2FC| >= 0.5, at least 70%
sample-direction consistency, and processed-matrix sign agreement. Enrichment
used the tested-gene universe and separated direction/contrast-specific sets.

Three predictive tasks were prespecified. Elastic Net was the primary model;
group-calibrated linear SVM and Random Forest were comparators. Validation used
5-fold outer and 4-fold inner StratifiedGroupKFold across 3 deterministic repeats. Every learned filter, feature selection, scale, hyperparameter, calibration, and compact panel was fitted inside training data. Compact panels of [5, 10, 15, 20, 30] genes and the fold-specific full stable signature were evaluated. Patient/donor bootstrap intervals quantified uncertainty.

The final Task C model and threshold were locked using GSE44076 only. The
prespecified primary cross-platform representation was within-sample percentile
rank over the label-independent common-gene universe. GSE41258 point estimates
used a deterministic one-array-per-patient subset; all arrays were retained for
patient-cluster bootstrap sensitivity. No external label informed model,
feature, representation, threshold, or hyperparameter selection. Primary seed:
44076.
