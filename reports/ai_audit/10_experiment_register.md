# Experiment register

| ID | Hypothesis | Baseline and change | Partition and seed | Result and decision |
|---|---|---|---|---|
| E01 | Cohort metadata and eligibility can be reconstructed from source records | Baseline tracked metadata versus independent GEO text parsing and clinical supplement | Full metadata, no model seed | Counts and compared fields agree after preserving source dates and normalizing missing sentinels |
| E02 | Current external metrics describe the archived locked predictions | Recompute 12 metric rows and compare cards, threshold table and publication tables | Saved GSE41258 prediction sets, no model fitting | Reconciled current JSON/Markdown; archived 19/14-gene artifacts |
| E03 | Demographics, composition and stress can account for some cohort discrimination | Five-gene fold-local selection and prespecified proxy, removal and residualization ablations | Three by five patient-grouped outer folds, seeds 44076–44078 | Perfect discrimination persists after gene-list exclusion. Residualized F1 0.843 and paired change −0.157 |
| E04 | A threshold selected within grouped training partitions improves an outer rank policy | Same fitted outer rank models scored with inner-selected cutoff and fixed 0.50 | Three by five outer folds, four threshold folds and four tuning folds, seeds 44076–44078 | Threshold 0.39–0.50; equal outer operating performance |
| E05 | GSE41258 ranking is robust to patient weighting while absolute-score levels shift | Replay 0.41 and prespecified 0.50, give each patient total weight one, and reproduce 1,000 bootstrap draws | All canonical archived predictions, bootstrap seed 44076 | ROC-AUC remains near 0.98. Rank specificity changes from 0.558 to 0.525 under patient weighting. All 13 primary intervals reproduce |
| E06 | Gene ranks compress cross-platform shifts | Compare the same two gene names using the verified deposited matrices and the 12,039-gene common annotation universe | Primary Task C development arrays, all external arrays with labels withheld | Large absolute shifts shrink after ranking. Assay scale and case mix remain confounded |
| E07 | A frozen Task B association model can transport to another early-onset cohort | Fit published mapped panel on GSE44076 deposited ranks. Evaluate GSE4107 once at 0.50 | Entire GSE44076 Task B development set, 22 external donors, seed 44076 | ROC-AUC 0.992, macro-F1 0.313. No evidence for operating-point transport |
| E09 | Raw-archive DE effects recur using deposited expression and a matched design | Refit limma U0/U1/U2 and paired models on the same mapped gene universe. Compare effect estimates and test an explicit minimum effect with treat | 148 unpaired or 196 paired arrays, deterministic | Effect correlations above 0.99998. U1 retained 1,578 archived hits with the same direction |

The machine-readable results under `results/ai_audit/` include source and
configuration digests, thresholds and fold-specific parameters. E03/E04 ran
in 756 seconds on the four-logical-processor Windows host. E07 ran in 65
seconds. E03 is sensitive to cohort-wide preprocessing and proxy choice. E04
is sensitive to compact panel instability and remains a development estimate.
E05 is conditional on the archived predictions and cohort. E07 has limited
precision and preserves the acquisition alias. The frozen-reference raw
preprocessing comparison and complete null refits require raw expression
reproduction and are pending.
