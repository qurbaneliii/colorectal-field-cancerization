# Figure legends

1. **Study design and pipeline.** Independent cohort preprocessing, biological
   analysis, patient-group nested cross-validation, signature locking, and the
   restricted external-validation claim.
2. **Sample-count flowchart.** Deposited and eligible samples for both GEO
   accessions, including exact primary tissue counts.
3. **Preprocessing QC.** Deposited normalized expression distributions,
   platform-specific PCA, and sample-correlation heatmaps.
4. **Primary PCA and clustering.** GSE44076 PCA and top-variable-gene heatmap
   across healthy, adjacent-normal, and tumor tissues.
5. **Adjacent-versus-healthy volcano.** Processed-matrix sensitivity contrast;
   red points pass FDR and effect-size thresholds. This is not the pending limma result.
6. **Paired tumor-versus-adjacent plot.** Top patient-matched processed-matrix
   contrasts; positive values indicate greater tumor expression.
7. **Field-signature heatmap and trajectories.** Top candidates after FDR,
   effect-size, and direction-consistency filtering.
8. **Nested-CV performance.** Untouched outer-fold macro F1 across models and
   prediction tasks.
9. **Out-of-fold confusion matrix, ROC/PR, and calibration.** Repeated
   predictions were averaged per sample before visualization.
10. **Stable coefficients.** Median standardized Elastic Net coefficients for
    frequently selected genes.
11. **External validation.** Locked Task C model evaluated in eligible GSE41258
    Primary Tumor and Normal Colon arrays using two cross-platform representations.
