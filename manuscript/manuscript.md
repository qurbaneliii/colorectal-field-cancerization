# Introduction

Colorectal carcinogenesis can be accompanied by molecular alterations beyond
the histologically apparent tumor, a concept commonly described as field
cancerization [REF]. Consequently, tumor-adjacent mucosa should not be assumed
to be equivalent to mucosa from cancer-free donors [REF]. Transcriptome-wide
profiling can identify candidate field alterations, but small-sample,
high-dimensional analyses are vulnerable to patient leakage, unstable feature
selection, platform effects, and optimistic validation [REF].

This study tests whether histologically normal tumor-adjacent colon differs
from genuinely healthy colon mucosa, characterizes healthy-to-adjacent-to-tumor
expression patterns without asserting biological progression, and evaluates
compact, interpretable signatures under nested patient-aware validation. A
separate external cohort tests only the tumor-versus-normal component.


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


# Results

## Cohorts and raw analysis

| accession   | platform   | tissue_class    |   arrays |   unique_patients_or_donors | role                             |
|:------------|:-----------|:----------------|---------:|----------------------------:|:---------------------------------|
| GSE44076    | GPL13667   | healthy         |       50 |                          50 | development                      |
| GSE44076    | GPL13667   | adjacent_normal |       98 |                          98 | development                      |
| GSE44076    | GPL13667   | tumor           |       98 |                          98 | development                      |
| GSE41258    | GPL96      | normal_colon    |       52 |                          52 | external tumor/normal validation |
| GSE41258    | GPL96      | primary_tumor   |      181 |                         181 | external tumor/normal validation |

Raw-CEL limma summary:

| comparison               | direction   |   tested_genes |   significant_genes |
|:-------------------------|:------------|---------------:|--------------------:|
| adjacent_vs_healthy      | down        |          11934 |                 659 |
| adjacent_vs_healthy      | up          |           6556 |                1003 |
| tumor_vs_adjacent_paired | down        |           9393 |                2025 |
| tumor_vs_adjacent_paired | up          |           9097 |                2477 |
| tumor_vs_healthy         | down        |          10401 |                2239 |
| tumor_vs_healthy         | up          |           8089 |                2764 |

## Internal predictive validation

| task                     |   mean_macro_f1 |   median_macro_f1 |   sd_macro_f1 |   fold_minimum |   fold_maximum |
|:-------------------------|----------------:|------------------:|--------------:|---------------:|---------------:|
| task_a_three_class       |        0.941281 |          0.951178 |     0.02413   |       0.89418  |       0.975118 |
| task_b_field_effect      |        0.992628 |          1        |     0.0152627 |       0.962677 |       1        |
| task_c_tumor_vs_adjacent |        0.984713 |          1        |     0.0188782 |       0.947222 |       1        |

The final compact Task C signature contained 5 genes:
FOXQ1, CEMIP, ETV4, GTF2IRD1, PACC1.

## External tumor-versus-normal validation

| representation                | evaluation_set    | is_primary   |   arrays |   unique_patients | analysis_provenance   |   f1_macro |   f1_weighted |   balanced_accuracy |   log_loss |   precision_adjacent_normal |   recall_adjacent_normal |   f1_adjacent_normal |   precision_tumor |   recall_tumor |   f1_tumor |   roc_auc |   pr_auc |   brier_score |   sensitivity |   specificity |   negative_predictive_value |   positive_predictive_value |   calibration_intercept |   calibration_slope |
|:------------------------------|:------------------|:-------------|---------:|------------------:|:----------------------|-----------:|--------------:|--------------------:|-----------:|----------------------------:|-------------------------:|---------------------:|------------------:|---------------:|-----------:|----------:|---------:|--------------:|--------------:|--------------:|----------------------------:|----------------------------:|------------------------:|--------------------:|
| within_sample_percentile_rank | canonical_patient | True         |      190 |               190 | raw_cel_rma           |   0.791895 |      0.874549 |            0.916667 |   0.411845 |                    0.509091 |                        1 |             0.674699 |                 1 |       0.833333 |   0.909091 |  0.996914 | 0.999447 |        0.1249 |      0.833333 |             1 |                    0.509091 |                           1 |                     nan |                 nan |

This external result does not validate cancer-free healthy versus
tumor-adjacent field cancerization.


# Discussion

The analysis separates the field-effect question from generic tumor detection.
Raw-CEL limma evidence and sample-level consistency provide the primary basis
for candidate field genes, while machine-learning coefficients contribute
predictive stability rather than defining biology. Patient-aware nesting limits
optimism from paired samples and feature selection.

Within-sample ranks offer a biologically neutral transfer representation: they
require no external labels, reduce dependence on platform intensity scale, and
preserve relative ordering. Strong discrimination should still be considered
separately from calibration and cannot establish clinical utility. Independent
cancer-free healthy and tumor-adjacent cohorts, tissue-composition analyses,
and prospective assay validation remain necessary.


# Limitations

- Retrospective public microarray cohorts cannot establish prospective clinical utility.
- GPL13667 and GPL96 differ in probe design and distribution.
- GSE41258 validates tumor versus normal colon, not the primary healthy-versus-adjacent field effect.
- GSE44076 is enriched for stage II and microsatellite-stable disease.
- Bulk expression conflates epithelial state with immune/stromal composition.
- Explicit batch covariates are incomplete and may remain confounded with phenotype.
- High-dimensional selection can remain unstable despite nested grouped validation.
- Differential expression and enrichment support association, not causal mechanism.
- The model is a candidate biomarker-discovery model, not a clinically ready diagnostic.


# Data and code availability

Expression data are available from NCBI GEO under GSE44076 and GSE41258. Raw
archives and large generated matrices are intentionally excluded from Git.
Checksums, deterministic metadata, compact result tables, environment locks,
commands, and model documentation are included in the repository. Users must
place the five named GEO files under `data/raw/<accession>/` and run `make all`.


# Figure legends

1. Study design and leakage-safe analysis flow.
2. Sample inclusion flow with array and patient counts.
3. Raw and post-RMA array QC distributions and diagnostics.
4. Raw-CEL RMA PCA, correlation, and hierarchical clustering.
5. Raw-CEL limma adjacent-normal versus healthy volcano plot.
6. Patient-fixed-effect tumor versus adjacent-normal volcano plot.
7. Raw-versus-deposited log2FC concordance.
8. Final field-candidate heatmap and group distributions.
9. Repeated nested-CV model comparison and three-class OOF confusion matrix.
10. Task B and Task C OOF ROC, precision-recall, and calibration.
11. Repeat-aware stable Elastic Net coefficients and compact-panel curve.
12. Patient-balanced GSE41258 ROC/PR, calibration, and confusion matrix.
13. Direction- and contrast-specific GO/Reactome/ranked enrichment results.
