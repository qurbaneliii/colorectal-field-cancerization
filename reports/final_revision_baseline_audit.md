# Final revision baseline audit

Audit date: 2026-08-01 (Asia/Baku)  
Baseline commit: `448c060a6433a24c7618f0274af60f25e9453fa7`  
Working branch: `feat/final-publication-methodology` (created from `origin/main`)  
Remote default branch: `main`  
Prior pull request: [#1](https://github.com/qurbaneliii/colorectal-field-cancerization/pull/1), merged  
Latest evidenced prior CI: workflow run `30650313294`, completed successfully for commit `a5a53d79f1fda4c9caeccb958ea01a34be1a2f46`

## Repository status before methodological edits

The latest remote `main` was inspected and checked out before this report was
created. The pre-existing checkout contained only the untracked `tmp/` folder;
it is unrelated to this revision and was left untouched. Running the required
data audit refreshed timestamps, the audited commit SHA, and manifest hashes in
three tracked evidence files. No scientific source file had been edited at the
time of this baseline assessment.

The repository is an existing hybrid R/Python analysis with raw-CEL
preprocessing, processed-matrix sensitivity analysis, grouped nested
cross-validation, enrichment, external validation, figures, reports, and a
manuscript scaffold. These valid components will be extended in place.

## Available data and generated artifacts

All five required GEO inputs are present locally and ignored by Git:

| Accession | Input | Bytes |
|---|---|---:|
| GSE44076 | `GSE44076_RAW.tar` | 532,213,760 |
| GSE44076 | `GSE44076_series_matrix.txt.gz` | 32,596,137 |
| GSE41258 | `GSE41258_RAW.tar` | 1,378,918,400 |
| GSE41258 | `GSE41258_series_matrix.txt.gz` | 17,537,584 |
| GSE41258 | `GSE41258_clinical_data.txt.gz` | 3,135 |

Local ignored artifacts include 636 extracted CEL files, both raw-CEL RMA
expression matrices, deposited-matrix expression matrices, six model files,
and full result tables. Approximate local footprint is 10.1 GB under `data/`,
281 MB under `results/`, and 0.22 MB under `models/`.

The refreshed data audit completed with exit code 0:

```powershell
.\.venv-publication\Scripts\python.exe scripts/run_data_audit.py
```

It confirmed 246 GSE44076 arrays (50 healthy, 98 adjacent-normal, 98 tumor),
98 complete tumor-adjacent patient pairs, 390 GSE41258 arrays, and 233 eligible
external tumor/normal arrays from 190 patients, with no expected-count
discrepancy.

## Current raw-CEL and QC status

Committed evidence reports successful platform-specific raw-CEL processing:
49,386 HG-U219 probe sets for GSE44076 and 22,283 HG-U133A probe sets for
GSE41258, followed by deterministic gene-level median aggregation. Local RMA
expression and ExpressionSet artifacts are present.

R 4.5.1 and Bioconductor 3.21 are installed at `D:\R\R-4.5.1`; the executable
is not on `PATH`. The locked R verification passed with exit code 0 when the
absolute executable was used, but `renv` reported that the project is
out-of-sync. `arrayQualityMetrics` 3.64.0 and `affyPLM` 1.84.0 are installed.

The current QC implementation has publication-critical deficiencies:

- pre- and post-RMA density figures use only the first array;
- no GSE41258 arrayQualityMetrics output exists, and the only GSE44076 HTML
  directory is local, ignored, and lacks a committed machine-readable summary;
- manuscript text nevertheless states that arrayQualityMetrics was generated;
- GSE41258 NUSE has not been produced;
- one GSE44076 sample has one extreme IQR flag, but no sample has the required
  two independent severe failures; no exclusion was applied and no explicit
  exclusion sensitivity report exists.

## Current biological inference

Raw-CEL limma tables contain 18,490 genes per contrast. At FDR < 0.05 and
|log2FC| >= 0.5, the unadjusted adjacent-versus-healthy analysis reports 1,662
genes (1,003 up and 659 down), and the paired tumor-versus-adjacent model
reports 4,502 genes (2,477 up and 2,025 down).

The primary field-effect inference is currently `expression ~ tissue` only.
Age and sex are available and associated with tissue group, but no
age/sex-adjusted limma model, covariate design-rank audit, adjusted/unadjusted
concordance, composition sensitivity, stress-gene audit, or evidence-tiered
field signature exists. The current 1,437-gene candidate table is therefore a
screening universe, not a defensible final field signature.

## Current predictive outputs

The existing repeated grouped nested CV uses three deterministic outer repeats
and keeps patients/donors grouped. Mean outer-fold Elastic Net macro-F1 is
0.9413 for Task A, 0.9926 for Task B, and 0.9847 for Task C. Task B does not yet
have its own final serialized model, compact signature, model card, nested
panel-size selection, or demographic/residualization sensitivity analysis.

The current compact-panel workflow constructs features within each outer
training fold, but chooses the global panel size using outer-test performance.
That remaining selection optimism prevents an unbiased final panel claim.

The current Task C five-gene compact consensus panel is FOXQ1, CEMIP, ETV4,
GTF2IRD1, and PACC1. In panel-specific selection, PACC1 appears in 8/15 folds
(0.533), below the configured 0.65 strict stability threshold, so the complete
five-gene panel cannot be labeled strictly stable. The general Elastic Net
stability table reports PACC1 in 10/15 folds (0.667), illustrating why the
stability estimand must be named explicitly.

## Current external validation

The current external workflow uses a deterministic one-array-per-patient set
(190 arrays), discarding one tissue from 43 patients with both normal-colon and
primary-tumor arrays. It rewrites the external biological labels to
`adjacent_normal` and `tumor`. The primary rank-transfer point estimates are
ROC-AUC 0.9969 and macro-F1 0.7919, while the training-z-score sensitivity
performs poorly. Four of five signature genes are present on GPL96; FOXQ1 is
absent. The workflow refits a four-gene transport model but does not yet expose
separate primary and transport artifacts or describe the refit completely.

The threshold is fixed at 0.5 rather than selected and locked from GSE44076
repeated OOF predictions. Patient-cluster bootstrap code exists, but the
primary point estimate is not based on the required canonical one-array-per-
patient-by-tissue set.

## Baseline tests and environments

All commands used the existing `.venv-publication` environment.

| Command | Result | Exit code |
|---|---|---:|
| `python -m pytest -q -m "not full_data"` | 6 passed, 17 deselected | 0 |
| `python -m pytest -q` | 23 passed | 0 |
| `python scripts/run_synthetic_smoke_test.py` | PASS; 96 predictions | 0 |
| `python -m pip check` | no broken requirements | 0 |
| `python -m ruff check src scripts tests` | all checks passed | 0 |
| `python -m compileall -q src scripts tests` | completed | 0 |
| `D:\R\R-4.5.1\bin\Rscript.exe R/verify_environment.R` | passed; renv out-of-sync warning | 0 |
| `docker info` | Linux engine pipe unavailable | 1 |

Docker verification is BLOCKED because the Docker Desktop Linux engine cannot
be reached at `npipe:////./pipe/dockerDesktopLinuxEngine`. WSL also reports
`Wsl/0x80070422` (required service disabled or unavailable). No Docker build or
run success is claimed.

## Current manuscript and report status

The manuscript has Introduction, Methods, Results, Discussion, Limitations,
Data and Code Availability, and Figure Legends files, but it remains a
generated scaffold. `[REF]` placeholders remain; required title page, abstract,
keywords, conclusion, declarations, references, acknowledgments, and
supplementary-materials sections are absent. The Methods falsely state that
arrayQualityMetrics was generated. `reports/final_analysis_report.md` is
identical to the readiness report, and multiple existing acceptance gates are
marked PASS based on file existence rather than executed evidence.

## Differences between claims and evidence

1. Cohort-wide density is claimed, but only the first array is plotted.
2. arrayQualityMetrics is claimed, but complete two-cohort evidence and summary
   tables are absent.
3. Task B is described as validated but lacks a distinct final workflow,
   compact signature, model card, and confounding sensitivity.
4. Compact-panel nesting is overstated because panel size is selected from
   outer-test results.
5. The five-gene Task C panel is described without clearly separating strict
   stability from compact-panel consensus frequency.
6. External terminology and sampling do not preserve `normal_colon` and
   `primary_tumor`, and the current primary set discards paired tissues.
7. The exact serialized five-gene model is not transportable when FOXQ1 is
   absent; a refit occurs but artifacts and wording do not fully distinguish it.
8. The default 0.5 threshold is not a primary-cohort OOF-derived locked
   threshold.
9. Performance estimands and their uncertainty are not separated consistently.
10. Tissue composition and preanalytical stress remain unassessed.
11. The manuscript has unverifiable placeholders and unsupported QC wording.
12. Publication/readiness PASS gates are inferred from artifact existence.

## Publication-critical issues to resolve

- correct all-array density plots and execute or honestly delimit AQM/NUSE;
- add full-rank U0/U1/U2 limma models, covariate audit, concordance, QC
  sensitivity, and evidence tiers;
- implement a validated microarray-compatible composition sensitivity and a
  documented immediate-early/stress audit;
- make Task B the distinct primary model with fold-local residualization,
  demographic/balancing sensitivities, fully nested panel-size selection,
  stability, calibration, and a model card;
- make Task C panel selection fully nested, name consensus versus strict
  stability accurately, derive its threshold exclusively from GSE44076 OOF
  predictions, preserve external labels, and use canonical patient-by-tissue
  evaluation with patient-cluster uncertainty;
- separate primary and transport model artifacts and report refitting exactly;
- separate outer-fold, repeat-level, aggregated-OOF, and bootstrap estimands;
- complete enrichment, signature intersection, publication figures/tables,
  manuscript, references, tests, and acceptance evidence;
- restore R lock consistency, retain Docker as BLOCKED until it genuinely runs,
  and verify the new GitHub CI run after push.
