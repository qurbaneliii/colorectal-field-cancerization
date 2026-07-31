# Limitations

- These are retrospective public microarray cohorts, not prospectively collected
  clinical-validation samples.
- GPL13667 and GPL96 differ in probe design and measurement distribution;
  within-sample ranks reduce but do not eliminate platform shift.
- GSE41258 validates tumor versus normal colon, not the healthy
  cancer-free-versus-adjacent field effect.
- Stage and molecular-subtype coverage is limited; GSE44076 is enriched for
  stage II, microsatellite-stable disease.
- Bulk-tissue expression may reflect cell-composition changes as well as
  epithelial field biology.
- No prospective or clinical-utility validation was performed.
- High-dimensional small-sample selection remains vulnerable to instability,
  even with nested grouped validation and stability thresholds.
- The executed environment lacked R, so raw-CEL RMA, final limma inference,
  arrayQualityMetrics, and enrichment remain pending.
- The signature is a biomarker-discovery candidate, not a clinically validated
  or clinically ready diagnostic.
