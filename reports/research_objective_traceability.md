# Research objective traceability

Baseline commit: `70998d6ffb537ddf720adb46ff217c373edce9fb`

## Primary biological objective: healthy versus adjacent-normal field effect

| Traceability element | Repository mapping |
|---|---|
| Objective | Determine whether histologically normal tumor-adjacent colorectal mucosa differs reproducibly from cancer-free healthy colon mucosa. |
| Required evidence | Raw-CEL preprocessing; metadata integrity; covariate-aware differential expression; robust field-gene tiers; composition and stress sensitivities; leakage-safe Task B prediction; feature stability; enrichment; uncertainty. |
| Implementation | `R/02_preprocess_gse44076.R`; `R/04_differential_expression.R`; `R/06_tissue_composition_sensitivity.R`; `scripts/run_covariate_sensitivity.py`; `scripts/run_stress_gene_sensitivity.py`; `scripts/run_modeling.py`; `scripts/run_compact_panel_analysis.py`; `scripts/run_task_b_confounding_sensitivity.py`. |
| Evidence | 246 aligned arrays; U0/U1/U2 limma tables; 101 high-confidence and 1,377 provisional genes; 74/101 composition-robust Tier 1 genes; five-gene Task B panel; nested-policy outer-fold macro-F1 1.000; demographic-only macro-F1 0.626; age/sex/location-residualized macro-F1 0.977. |
| Status | PARTIAL at baseline. The central evidence exists and is strong, but Task B composition residualization uses ten correlated scores rather than the prespecified limited PC representation, and the required per-gene concordance artifact is missing. |

## Supporting three-state tissue objective: healthy versus adjacent-normal versus tumor

| Traceability element | Repository mapping |
|---|---|
| Objective | Establish that the three tissue states are distinguishable without treating this as the primary endpoint. |
| Required evidence | Three-class grouped nested CV; no patient leakage; calibrated probability metrics; supporting biological trajectory analyses; conditional permutation null. |
| Implementation | `scripts/run_modeling.py`; `scripts/run_compact_panel_analysis.py`; `src/modeling/*`; `scripts/run_permutation_tests.py`; trajectory tables. |
| Evidence | Nested-policy outer-fold macro-F1 0.948; aggregated OOF Elastic Net interval reported separately; conditional 1,000-permutation p=0.000999. |
| Status | PASS. Task A is consistently supporting rather than primary. |

## Secondary tumor-transition objective: tumor versus adjacent-normal

| Traceability element | Repository mapping |
|---|---|
| Objective | Derive an interpretable tumor-versus-adjacent signature for secondary translational evaluation. |
| Required evidence | Patient-aware paired DE; grouped nested Task C models; fully nested panel selection; strict stability terminology; primary-only threshold locking. |
| Implementation | `R/04_differential_expression.R`; `scripts/run_modeling.py`; `scripts/run_compact_panel_analysis.py`; `scripts/select_task_c_threshold.py`. |
| Evidence | 4,502 paired DE genes; three-gene Task C panel FOXQ1/CEMIP/ETV4; all genes strict-stability frequency and sign consistency 1.0; nested-policy outer-fold macro-F1 0.990; locked threshold 0.41. |
| Status | PASS scientifically; PARTIAL artifact contract because canonical Task C filenames are absent at baseline. |

## External-transfer objective: GSE44076 to GSE41258

| Traceability element | Repository mapping |
|---|---|
| Objective | Evaluate a locked tumor-versus-normal signature across HG-U219 and HG-U133A without presenting it as external validation of healthy-versus-adjacent field cancerization. |
| Required evidence | Independent platform normalization; label-independent common-gene intersection; primary-data-only threshold; distinct transport refit; original biological labels; canonical patient-tissue evaluation; patient-cluster uncertainty; calibration and platform-shift sensitivities. |
| Implementation | `R/03_preprocess_gse41258.R`; `src/data/cross_platform.py`; `scripts/select_task_c_threshold.py`; `scripts/run_external_validation.py`. |
| Evidence | 233 eligible/canonical arrays from 190 patients; 43 paired-tissue patients; FOXQ1 absent and two genes transported; ROC-AUC 0.9829, macro-F1 0.8282, specificity 0.5577, sensitivity 1.0, Brier 0.0533; 1,000 patient-cluster bootstrap intervals. |
| Status | PASS. Evidence is strong for discrimination but only moderate for threshold/probability transport; it does not validate Task B. |

## Reproducibility objective: same code, data, and environment to reproducible results

| Traceability element | Repository mapping |
|---|---|
| Objective | Make every principal claim reproducible from versioned code, immutable input checksums, locked dependencies, deterministic seeds, tests, and auditable artifacts. |
| Required evidence | Python and R locks; analysis configuration; input hashes; executable full route; unit/full-data tests; synthetic smoke; CI; Docker execution; model/result provenance. |
| Implementation | `requirements-lock.txt`; `renv.lock`; `config/analysis.yaml`; `Makefile`; `Dockerfile`; `.github/workflows/*`; tests; artifact manifest. |
| Evidence | Python dependency/compile/lint checks pass; 37 tests pass; R environment check passes; synthetic smoke passes; GitHub CI run 30693185908 passes. Docker engine is unavailable. Model cards record an older generation commit and lack a source hash/dirty-state field. |
| Status | PARTIAL. Local and CI reproducibility are evidenced; Docker is BLOCKED and artifact provenance requires regeneration after code stabilization. |

## Objective boundary

The repository supports two separate claims:

1. Primary: internally reproducible transcriptomic differences and predictive
   information distinguish cancer-free healthy mucosa from tumor-adjacent
   histologically normal mucosa in GSE44076.
2. Secondary: a locked tumor-transition signature transfers with strong ranking
   discrimination but incomplete threshold/calibration transport from GSE44076
   to GSE41258.

GSE41258 does not contain a cancer-free healthy-versus-adjacent design and is
therefore not external validation of the primary field-effect endpoint.
