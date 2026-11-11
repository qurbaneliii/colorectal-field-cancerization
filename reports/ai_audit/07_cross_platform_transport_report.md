# Task C cross-platform transport

The canonical GSE41258 replay used 233 saved predictions from 190 patients.
No model was fitted and no operating point was optimized on the replayed
labels. The analysis reproduces all 12 metric rows in the canonical results,
including named sensitivity estimands. GSE41258 has already informed the
repository's research history.

| Representation and threshold | ROC-AUC | Average precision | Macro-F1 | Balanced accuracy | Sensitivity | Specificity | Brier |
|---|---:|---:|---:|---:|---:|---:|---:|
| Percentile rank, 0.41 development lock | 0.9829 | 0.9935 | 0.8282 | 0.7788 | 1.0000 | 0.5577 | 0.0533 |
| Percentile rank, 0.50 sensitivity | 0.9829 | 0.9935 | 0.8898 | 0.8530 | 0.9945 | 0.7115 | 0.0533 |
| Training z-score, 0.41 lock applied as sensitivity | 0.9822 | 0.9928 | 0.4372 | 0.5000 | 1.0000 | 0.0000 | 0.1962 |

At 0.41, the rank model had 29 true negatives, 23 false positives, zero false
negatives and 181 true positives. At 0.50, it had 37 true negatives and one false
negative. With equal total weight per patient, locked specificity was 0.5246
and sensitivity was 1.0. Forty-two of 43 paired patients had a higher tumor
score with both representations.

The rank score averaged 0.8234 versus a tumor proportion of 0.7768. The joint
calibration intercept and slope were −2.066 and 2.620. Calibration-in-the-large,
estimated with slope fixed to one, was −0.492. The rank Brier score beat the
0.1734 prevalence-only benchmark. The archived 95% patient-cluster intervals
were 0.417–0.692 for specificity, −10.768 to −1.392 for the joint intercept,
and 1.875–12.632 for the slope. With training z-scores, even the minimum
normal-colon probability was 0.8833. Average probability was 0.9795 and the Brier
score was worse than the prevalence benchmark. Both representations retained
strong ranking. The absolute-expression model placed the useful ranking
entirely above either decision threshold.

Separate GEO-deposited expression diagnostics showed two mapped probes per
gene on HG-U219 and one per gene on HG-U133A. For CEMIP and ETV4 the mean
absolute shift was 128.7 and 93.0 primary SDs, compared with 1.08 and 0.69 SD
for within-array ranks. GEO records describe PLIER summaries for GSE41258
and RMA normalization for GSE44076. The deposited scales differ markedly and the
external expression-only diagnostic uses all 390 arrays, including samples
outside the Task C endpoint. These values cannot quantify the raw-RMA platform
effect or distinguish probe effects from expression scale and tissue mix.
The locked threshold remains the historical primary operating point.

The score and weighting replays are in
`results/ai_audit/external_prediction_replay.csv` and
`results/ai_audit/external_score_quantiles.csv`. Archived patient-cluster
intervals are in `results/metrics/external_validation_bootstrap_ci.csv`.
The independent 1,000-draw patient-cluster replay reproduced all 13 primary
intervals to within 1.5 × 10⁻¹⁴ in endpoint differences. See
`results/ai_audit/external_rank_cluster_bootstrap_replay.csv`.

A separate fresh raw-CEL refit with the archived panel and hyperparameters
reselected 0.41 using GSE44076 only. The new GSE41258 predictions reproduced
the primary metrics and all 13 patient-cluster interval endpoints.
[Raw-CEL continuation](15_raw_reproduction_continuation.md)

Calibration diagnostics follow [Van Calster et al.](https://doi.org/10.1186/s12916-019-1466-7).
