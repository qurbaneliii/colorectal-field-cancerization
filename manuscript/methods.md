# Methods

## Datasets and eligibility

GSE44076 (GPL13667) was the sole development cohort: 50 cancer-free healthy
donors, 98 tumor-adjacent normal samples, and 98 tumors from the same 98
patients. GSE41258 (GPL96) was normalized and annotated independently. Its main
external subset was restricted before analysis to author-included Primary Tumor
and Normal Colon arrays; metastases, non-colon normal tissues, polyps, cell
lines, ambiguous samples, and titles matching `_ez` or `rehyb` were excluded.

## Preprocessing and quality control

The executed Python sensitivity route used the normalized expression deposited
in the GEO series matrices. Official GEO platform annotation was fetched
independently; probes with absent or ambiguous gene-symbol mapping were removed,
and multiple unambiguously mapped probes were aggregated per sample by the
median. The two platforms were never jointly normalized. Distribution,
PCA, correlation, missingness, and robust sample-median outlier-candidate checks
were generated; separation alone never triggered exclusion.

Raw-CEL RMA scripts use `oligo`/`pd.hg.u219`/`hgu219.db` for GSE44076 and
`affy`/`hgu133a.db` for GSE41258, followed by label-independent median
probe-to-gene aggregation, arrayQualityMetrics, limma, and clusterProfiler.
They were not executed in this environment because Rscript was unavailable.
Consequently, inferential differential-expression tables in this run are
explicitly labeled processed-matrix Welch or paired-t sensitivity analyses and
must not be represented as the final limma analysis.

## Statistical and trajectory analysis

Adjacent-normal versus healthy and tumor versus healthy sensitivity contrasts
used Welch tests. Tumor versus adjacent-normal used patient-matched paired
tests. Benjamini-Hochberg correction was applied. Field candidates required
FDR < 0.05, absolute mean log2 difference ≥
0.5, and at least 70% sample-direction consistency.
Trajectory categories were assigned by explicit rules in `src/biology.py` with
a 0.15 log2 tolerance.

## Predictive modeling and leakage prevention

Three objectives were analyzed separately: healthy/adjacent/tumor,
healthy/adjacent, and adjacent/tumor. Repeated 5-fold
outer and 4-fold inner StratifiedGroupKFold
splits used 2 deterministic seeds. Patient
pairs were indivisible; healthy donors had unique groups. Zero-variance
removal, training-variance filtering, univariate supervised selection,
standardization, tuning, and fitting all occurred inside the training pipeline.
Elastic Net logistic regression was compared with a linear SVM and Random
Forest. Macro F1 was primary. Out-of-fold predictions and every assignment were
retained. No SMOTE or neural network was used.

Feature stability was the outer-fold selection frequency plus coefficient-sign
consistency. Elastic Net was prespecified as preferred when performance was
within uncertainty of a complex comparator. The Task C panel and threshold
were locked from GSE44076. GSE41258 labels were never used for tuning. External
evaluation compared training-derived scaling with within-sample percentile
ranks computed over the common-gene universe. Confidence intervals used
patient/donor-group bootstrapping.

## Software and seeds

Python 3.12.10 on Windows-11-10.0.26200-SP0; scikit-learn
1.8.0. Primary seed: 44076;
additional deterministic seeds: [44076, 44077].
