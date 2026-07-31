# Results

## Cohorts and processing

The primary audit recovered exactly 246 arrays: 50 healthy, 98 adjacent-normal,
and 98 tumor, with 98 complete patient pairs and no count discrepancy. Official
platform annotation and deterministic aggregation yielded 19,040
GSE44076 genes. GSE41258 yielded
13,299
genes; 11,336
symbols were common.

The processed-matrix sensitivity analysis identified 1,654
adjacent-versus-healthy candidates meeting the prespecified FDR, effect, and
direction-consistency criteria. These are provisional until raw-CEL
Bioconductor preprocessing and limma are executed.

## Internal validation

Elastic Net outer-fold mean macro F1 values were:

| task                     | model       |   f1_macro |
|:-------------------------|:------------|-----------:|
| task_a_three_class       | elastic_net |   0.952751 |
| task_b_field_effect      | elastic_net |   0.992605 |
| task_c_tumor_vs_adjacent | elastic_net |   0.987218 |

All reported predictions were generated for untouched patient-group outer
folds. Adjacent-normal recall, per-class metrics, calibration scores, and all
fold-level distributions are in `results/metrics/nested_cv_metrics.csv`.

## Signature and external validation

The locked Task C stable signature contained
19
genes before cross-platform intersection. GSE41258 main validation contained
233 eligible arrays.

| representation                |   roc_auc |   pr_auc |   balanced_accuracy |   f1_macro |   sensitivity |   specificity |
|:------------------------------|----------:|---------:|--------------------:|-----------:|--------------:|--------------:|
| training_zscore               |  0.992775 | 0.997055 |            0.870166 |   0.769768 |      0.740331 |      1        |
| within_sample_percentile_rank |  0.99575  | 0.998633 |            0.987622 |   0.987622 |      0.994475 |      0.980769 |

This validates only the tumor-versus-normal-colon component. GSE41258 does not
provide an independent cancer-free healthy-versus-adjacent validation cohort.
