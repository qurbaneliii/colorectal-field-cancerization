# Before and after evaluation

The archived raw-CEL pipeline and the new deposited-matrix experiments use
different preprocessing and, in some cases, different panel policies.
Comparisons across those routes establish sensitivity, not matched changes
in predictive performance.

| Dimension | Original raw-CEL snapshot | Audited result | Status |
|---|---|---|---|
| Task A nested compact macro-F1 | 0.9482 | Checked archived summary | Full re-estimation pending |
| Task B nested compact macro-F1 | 1.000 | 1.000 in separate deposited fold-selected sensitivity | Cohort association reproduced in a separate route |
| Task B composition robustness | Fixed-panel adjusted 0.8497 | Fold-reselected sensitivity 0.8431 with interval 0.7769–0.9008 | Comparable qualitative finding, different estimands |
| Task B stability | CLC/FOS/VIP compact selections 15/12/14 of 15; DYNC1H1 7 | Four named genes selected 15 of 15 in unadjusted deposited sensitivity | Same-cohort robustness, different mapping and policy |
| Task B external | Primary model internal only | GSE4107 four-gene refit ROC-AUC 0.9917, macro-F1 0.3125 at fixed 0.50 | New confounded cohort-association sensitivity |
| Task C nested compact macro-F1 | 0.9898 | 0.9898 in a separate nested rank policy | Distinct estimands |
| Task C external ROC-AUC | 0.9829 | 0.9829 from locked prediction replay | Historical estimate verified |
| Task C external macro-F1 | 0.8282 | 0.8282 from locked prediction replay | Historical estimate verified |
| Task C external specificity | 0.5577 | 0.5577 at historical lock, 0.5246 equal patient weighting | New weighting sensitivity |
| Task C external calibration | Brier 0.0533; intercept −2.066, slope 2.620 | Replayed; calibration-in-the-large −0.492 | Diagnostic only |
| Data integrity | GEO scan dates dropped | Source dates present on 233 of 233 eligible external arrays | Verified correction |
| Leakage and selection | Cohort-wide RMA and development-conditional threshold OOF | Explicit estimands and tested nested threshold workflow | Raw preprocessing comparison pending |
| Reproducibility and artifacts | Contradictory active 19/14- and 3/2-gene reports | Historical artifacts retained and active metrics replay contract in CI | Verified artifact correction |
