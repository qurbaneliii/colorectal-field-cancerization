# Compact-panel selection rationale

Candidate panel construction was repeated independently inside every outer
training fold. Genes were ranked from the inner-selected Elastic Net fitted to
outer-training data only; each restricted model was then evaluated on the
untouched outer test fold. No global panel was used to estimate internal
performance.

| task                     |   selected_panel_size |   full_signature_macro_f1 |   selected_panel_macro_f1 | rule                                                                      | analysis_provenance   |
|:-------------------------|----------------------:|--------------------------:|--------------------------:|:--------------------------------------------------------------------------|:----------------------|
| task_a_three_class       |                     5 |                  0.943201 |                  0.944373 | Smallest panel within 0.020 macro F1 of the full fold-specific signature. | raw_cel_rma           |
| task_b_field_effect      |                     5 |                  0.992628 |                  1        | Smallest panel within 0.020 macro F1 of the full fold-specific signature. | raw_cel_rma           |
| task_c_tumor_vs_adjacent |                     5 |                  0.984713 |                  0.994909 | Smallest panel within 0.020 macro F1 of the full fold-specific signature. | raw_cel_rma           |

The final Task C gene list was derived only after panel-size performance was
estimated, by aggregating training-fold rankings across all repeat/fold keys.
Selected Task C genes: FOXQ1, CEMIP, ETV4, GTF2IRD1, PACC1.
