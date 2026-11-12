# Remaining limitations

- The primary Task B cohort cannot identify a cancer-specific field signal
  separately from acquisition. The small independent GSE4107 cohort repeats
  the same protocol split.
- The closest acquisition-matched candidate requires approved JGA access.
  Cross-site and stage overlap, mapping and cohort provenance need auditing
  before any locked evaluation.
- The effect of frozen-reference preprocessing on held-out raw-array scores
  is unknown. Both raw-CEL RMA and deposited normalized-expression analyses
  inherit cohort-wide preprocessing.
- GSE41258 is an exposed neoplasia cohort. The exact three-gene model lacks
  coverage on GPL96. The two-gene refit has good ranking and weak transport of
  the operating threshold. External calibration is measured on a particular
  case mix and carries substantial patient-level uncertainty.
- Marker-based composition scores and stress-gene ablations cannot distinguish
  cell-state changes from handling and tissue mixture. Biological mediation
  and confounding remain entangled.
- Repeated grouped folds share training patients. Bootstrap intervals of
  averaged OOF scores are conditional on the fitted models. Training and
  selection uncertainty and transport calibration require new data.
- Historical TAR byte identity requires the original inputs or historical CEL
  fingerprints. The current GEO archives support a separately pinned raw-CEL
  run. Clinical utility remains open.
  [Raw-CEL continuation](15_raw_reproduction_continuation.md)
