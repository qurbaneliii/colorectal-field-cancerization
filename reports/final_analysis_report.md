# Final scientific analysis report

## Scientific question and data

This retrospective study separates the primary healthy-versus-adjacent field
question (Task B), supporting three-class discrimination (Task A), and secondary
tumor-normal cross-platform transfer (Task C). GSE44076 contributes 246 arrays;
GSE41258 contributes 233 canonical external arrays from
190 patients.

## Executed findings

- QC: all-array densities and arrayQualityMetrics were executed for both cohorts;
  no primary exclusion changed. GPL96 NUSE was computed with affyPLM.
- Differential expression: U0/U1/U2 detected
  1,662/
  1,580/
  1,521 genes;
  paired tumor-adjacent analysis detected
  4,502.
- Field evidence: 101 high-confidence and 1,377
  provisional genes; 74 high-confidence genes were composition robust.
- Biology: ECM, collagen/proteoglycan turnover, vascular/stromal activation, and
  provisional metabolic attenuation were the principal nonredundant themes.
- Task A: aggregated OOF Elastic Net macro-F1 0.959 (0.929-0.982).
- Task B: CLC, DYNC1H1, FOS, VIP, SNORA12; nested-policy mean outer-fold macro-F1
  1.000; aggregated OOF ROC-AUC
  1.000 (1.000-1.000). No independent external Task B cohort exists.
- Task B confounding: macro-F1 was
  0.977
  after demographic residualization and
  0.718
  after composition-plus-demographic residualization.
- Task C: FOXQ1, CEMIP, ETV4; every gene is strictly stable. Threshold
  0.41 came only from grouped GSE44076 OOF rank
  predictions. Only CEMIP, ETV4 transferred to GPL96.
- External Task C: ROC-AUC 0.985 (0.957-1.000), macro-F1
  0.831 (0.755-0.889), specificity
  0.561 (0.417-0.692); strong ranking coexisted with
  incomplete threshold and calibration transport.

## Interpretation and limitations

The data support an internally validated, cohort-specific field-effect classifier
and a distinct externally evaluated tumor-normal transfer model. Bulk-tissue
composition, residual cohort confounding, missing external Task B validation,
platform-specific gene loss, and retrospective design preclude clinical or
causal claims. Full methodological detail and uncertainty are in the manuscript;
machine-readable evidence is indexed by `reports/result_artifact_manifest.csv`.
