# Final analysis report

## 1. Executive summary

The repository now implements a patient-aware, feature-leakage-safe hybrid
R/Python field-cancerization workflow. Primary counts and all 98 pairs were
verified exactly. The Python deposited-matrix route, nested validation,
signature locking, and external validation were executed. Raw-CEL Bioconductor
stages remain blocked solely because R is absent and are not claimed complete.

## 2. Scientific question

Whether histologically normal tumor-adjacent colon carries a reproducible
transcriptomic field effect relative to genuinely healthy mucosa, and whether
compact explainable signatures distinguish healthy, adjacent, and tumor states
without patient leakage.

## 3. Dataset audit and sample inclusion

| accession   | platform   | tissue_class    |   included_samples | role                |
|:------------|:-----------|:----------------|-------------------:|:--------------------|
| GSE44076    | GPL13667   | healthy         |                 50 | development         |
| GSE44076    | GPL13667   | adjacent_normal |                 98 | development         |
| GSE44076    | GPL13667   | tumor           |                 98 | development         |
| GSE41258    | GPL96      | normal_colon    |                 52 | external validation |
| GSE41258    | GPL96      | primary_tumor   |                181 | external validation |

Detailed deterministic exclusions are in `reports/metadata_audit.md`.

## 4. Preprocessing and quality control

The executed route used independently deposited normalized GEO matrices,
official platform annotations, unambiguous probe mapping, and per-sample median
aggregation. It produced 300-DPI PNG plus SVG/PDF QC outputs. No outlier was
automatically removed. Raw CEL files were checksummed but not modified.

## 5. Differential expression and field findings

| comparison                      | direction   |   tested_genes |   significant_genes |
|:--------------------------------|:------------|---------------:|--------------------:|
| adjacent_normal_vs_healthy      | down        |          12298 |                 673 |
| adjacent_normal_vs_healthy      | up          |           6742 |                1027 |
| tumor_vs_adjacent_normal_paired | down        |           9657 |                2079 |
| tumor_vs_adjacent_normal_paired | up          |           9383 |                2566 |
| tumor_vs_healthy                | down        |          10711 |                2302 |
| tumor_vs_healthy                | up          |           8329 |                2844 |

Field candidates passing the provisional processed-matrix rules: 1654.
Final limma inference is pending R/Bioconductor execution.

## 6. Model methodology and leakage prevention

All patient/donor groups were disjoint in every fold. Every learned predictive
operation was fitted inside the inner/outer training boundary. Three tasks,
Elastic Net, linear SVM, and Random Forest were evaluated; Elastic Net was
selected under the prespecified interpretability/stability rule.

## 7. Internal validation, stability, and explainability

| task                     | model       |   f1_macro |
|:-------------------------|:------------|-----------:|
| task_a_three_class       | elastic_net |   0.952751 |
| task_b_field_effect      | elastic_net |   0.992605 |
| task_c_tumor_vs_adjacent | elastic_net |   0.987218 |

Fold predictions, metrics, assignments, coefficient stability, compact panels,
and the reload-tested locked model are saved under `results/` and `models/`.

## 8. External validation

| representation                |   roc_auc |   pr_auc |   balanced_accuracy |   f1_macro |   sensitivity |   specificity |
|:------------------------------|----------:|---------:|--------------------:|-----------:|--------------:|--------------:|
| training_zscore               |  0.992775 | 0.997055 |            0.870166 |   0.769768 |      0.740331 |      1        |
| within_sample_percentile_rank |  0.99575  | 0.998633 |            0.987622 |   0.987622 |      0.994475 |      0.980769 |

The external claim is restricted to tumor versus normal colon.

## 9. Statistical uncertainty and permutation

Patient/donor-group bootstrap confidence intervals are stored in
`results/metrics/`. Group-preserving permutation results:

| task                     |   observed_mean_macro_f1 |   null_mean |   null_standard_deviation |   null_95_percentile |   permutation_p_value |   iterations | permutation_unit               |
|:-------------------------|-------------------------:|------------:|--------------------------:|---------------------:|----------------------:|-------------:|:-------------------------------|
| task_a_three_class       |                 0.952917 |    0.664983 |                 0.02622   |             0.702991 |            0.00990099 |          100 | patient/donor group preserving |
| task_b_field_effect      |                 1        |    0.471956 |                 0.0475771 |             0.562994 |            0.00990099 |          100 | patient/donor group preserving |
| task_c_tumor_vs_adjacent |                 0.989444 |    0.492204 |                 0.0574963 |             0.58449  |            0.00990099 |          100 | patient/donor group preserving |

## 10. Functional enrichment

Correct-background, BH-adjusted GO enrichment is implemented in
`R/05_functional_enrichment.R`; it was not executed because R is unavailable.
No pathway mechanism is fabricated.

## 11. Reproducibility

Exact executed command sequence:

```powershell
.\.venv\Scripts\python scripts/run_data_audit.py
.\.venv\Scripts\python scripts/run_processed_matrix_pipeline.py
.\.venv\Scripts\python scripts/run_modeling.py
.\.venv\Scripts\python scripts/run_permutation_tests.py
.\.venv\Scripts\python scripts/run_external_validation.py
.\.venv\Scripts\python scripts/build_manuscript_outputs.py
.\.venv\Scripts\python -m pytest -q
```

Complete route after R is installed: `make all`.

The committed R lock is intentionally a bootstrap lock because R was absent in
this execution. Run `R/00_install_packages.R` and commit the resolved
package-filled `renv.lock` before the manuscript environment is frozen.

## 12. Remaining scientific risks and next steps

Execute raw-CEL RMA/QC/limma/enrichment in the pinned Docker/R environment,
compare deposited-matrix sensitivity results with raw-CEL results, validate the
field-effect panel in an independent cancer-free/adjacent cohort, assess
cell-composition and molecular-subtype confounding, and perform prospective
assay validation before any clinical claim.
