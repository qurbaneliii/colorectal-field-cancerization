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
