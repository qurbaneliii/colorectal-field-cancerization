# Transcriptomic field cancerization in colorectal cancer

Publication-oriented reanalysis of **GSE44076** with raw Affymetrix CEL
files. **GSE41258** is used only for locked tumor-versus-normal-colon external
validation; it cannot independently validate the healthy-versus-adjacent field-effect
task.

The executed primary route is:

- R 4.5.1 / Bioconductor 3.21 raw-CEL RMA and platform-specific annotation;
- limma inference, including patient-fixed-effect paired contrasts;
- patient/donor-grouped repeated nested cross-validation in Python;
- fold-internal compact-panel selection and group-aware SVM calibration;
- GSE44076-locked rank threshold and patient-clustered cross-platform evaluation;
- 1,000-iteration task-specific permutation tests and grouped bootstrap intervals.

This is retrospective biomarker discovery, not a clinical-readiness claim. Three tissue
group means are not interpreted as longitudinal progression.

## Executed result snapshot

- GSE44076: 246 HG-U219 arrays (50 healthy, 98 adjacent, 98 tumor; 98 complete pairs).
- GSE41258: 390 HG-U133A arrays audited; 233 tumor/normal arrays from 190 patients are
  eligible for external validation.
- Raw RMA assertions: 49,386 GSE44076 and 22,283 GSE41258 probe sets.
- High-confidence field signature: **101 genes**; 74 retain
  composition-adjusted support and 3 overlap the curated immediate-early/stress set.
- Primary Task B compact panel: **CLC, DYNC1H1, FOS, VIP, SNORA12**. The primary
  raw-CEL estimate is internal. Validation under matched acquisition protocols
  remains outstanding.
- Task B confounding challenge: macro-F1 was 0.977 after fold-local age/sex/location
  residualization and 0.850 after demographics plus two composition PCs learned
  inside each outer-training fold (demographic-only baseline 0.626).
- Secondary Task C compact panel: **FOXQ1, CEMIP, ETV4**. All three meet the strict
  stability rule; only **CEMIP** and **ETV4** are present on GPL96.
- External Task C primary analysis: 233 canonical patient-tissue arrays from 190
  patients; patient-clustered ROC AUC 0.9829 and macro-F1 0.8282 at the
  GSE44076-locked threshold of 0.41. The two-gene transport refit is distinct from the
  exact three-gene primary model.
- Reproduction and scientific readiness are documented in the independent audit.

See [the final readiness report](reports/publication_readiness_final.md),
[the manuscript](manuscript/manuscript.md), and
[the artifact manifest](reports/result_artifact_manifest.csv).
The audit trail starts with
[objective alignment](reports/objective_alignment_baseline_audit.md) and the
[research traceability matrix](reports/research_objective_traceability.md).
The [raw-CEL reproduction continuation](reports/ai_audit/15_raw_reproduction_continuation.md)
records the historical and current GEO archive identities, fresh runs, and
remaining reproducibility blockers.

## Reproduce

The supported interpreter versions are Python 3.12 and R 4.5.1. On Windows, invoke the
installed R executable explicitly if `Rscript` is not on `PATH`.

```powershell
py -3.12 -m venv .venv-publication
.\.venv-publication\Scripts\python -m pip install -r requirements-lock.txt
.\.venv-publication\Scripts\python -m pip install --no-deps -e .

D:\R\R-4.5.1\bin\Rscript.exe R/00_install_packages.R
D:\R\R-4.5.1\bin\Rscript.exe R/verify_environment.R
.\.venv-publication\Scripts\python scripts/run_data_audit.py
D:\R\R-4.5.1\bin\Rscript.exe R/01_extract_and_audit_cel.R
D:\R\R-4.5.1\bin\Rscript.exe R/02_preprocess_gse44076.R
D:\R\R-4.5.1\bin\Rscript.exe R/03_preprocess_gse41258.R
.\.venv-publication\Scripts\python scripts/convert_r_outputs.py
D:\R\R-4.5.1\bin\Rscript.exe R/04_differential_expression.R
.\.venv-publication\Scripts\python scripts/run_processed_matrix_pipeline.py
.\.venv-publication\Scripts\python scripts/run_raw_vs_processed_sensitivity.py
.\.venv-publication\Scripts\python scripts/run_covariate_sensitivity.py
D:\R\R-4.5.1\bin\Rscript.exe R/06_tissue_composition_sensitivity.R
.\.venv-publication\Scripts\python scripts/run_stress_gene_sensitivity.py
.\.venv-publication\Scripts\python scripts/run_modeling.py --provenance raw_cel_rma
.\.venv-publication\Scripts\python scripts/run_compact_panel_analysis.py --provenance raw_cel_rma
D:\R\R-4.5.1\bin\Rscript.exe R/05_functional_enrichment.R
.\.venv-publication\Scripts\python scripts/run_task_b_confounding_sensitivity.py
.\.venv-publication\Scripts\python scripts/run_permutation_tests.py --provenance raw_cel_rma
.\.venv-publication\Scripts\python scripts/select_task_c_threshold.py
.\.venv-publication\Scripts\python scripts/run_external_validation.py --provenance raw_cel_rma
.\.venv-publication\Scripts\python scripts/build_publication_tables.py
.\.venv-publication\Scripts\python scripts/build_publication_figures.py
.\.venv-publication\Scripts\python scripts/build_manuscript_outputs.py
.\.venv-publication\Scripts\python -m pip check
.\.venv-publication\Scripts\ruff check src scripts tests
.\.venv-publication\Scripts\python -m pytest -q
```

`make all` expresses the same full route on systems with GNU Make. `make smoke` runs a
no-data grouped-CV smoke test. `make sensitivity-all` reproduces the explicitly labeled
GEO-deposited-series-matrix sensitivity route.

Stage verified archives in an isolated checkout with
`scripts/stage_geo_inputs.py`. It reuses a checksum-verified cache and defaults
to the historical manifest. The current GEO pins and the explicit opt-in
command are documented in the continuation audit.

The executed QC bundle includes cohort-wide pre/post distributions, RLE, MA, PCA,
sample correlation, hierarchical clustering, objective multimetric exclusion flags,
and `arrayQualityMetrics` HTML reports for both cohorts. GPL96 also has affyPLM NUSE;
no surrogate NUSE is used for the oligo platform.

## Provenance and safeguards

Raw inputs under `data/raw/` are never modified. Their paths, sizes, SHA-256 hashes,
archive-member counts, GEO mappings, and platform checks are recorded under
`data/metadata/`. Raw-CEL and deposited-matrix artifacts use separate filenames and
result namespaces.

Model splits preserve patient and donor groups. Nested CV fits variance filtering,
univariate selection, scaling, tuning, calibration and panel ranking within training
partitions. RMA was fit on the entire cohort before those splits. The historical
Task C threshold uses grouped development scores conditional on a panel and
hyperparameters selected on all development samples. The external scores and
threshold were frozen for evaluation.

The locked environments are `requirements-lock.txt` and `renv.lock`. CI runs lint,
no-data tests, and a synthetic grouped-CV smoke test. The manual full-data workflow
supports a self-hosted runner with local inputs or `github-hosted` with the
`current-geo` pins. Hosted runs retain the input hashes, stage logs, generated
publication outputs and their checksums. The optional `analysis_sha` input replays
an exact prior commit with matching pins and dependency locks.

The [independent scientific audit](reports/ai_audit/FINAL_EXECUTIVE_RESEARCH_REPORT.md)
documents sampling-protocol confounding, validation estimands and real-data
sensitivities. Run `python scripts/verify_external_results.py` to reconcile the
active external reports with saved predictions.
[Final acceptance code review](reports/ai_audit/16_final_acceptance_review.md)

## Data availability

The raw and deposited expression inputs are public through NCBI GEO (GSE44076 and
GSE41258). Large raw/intermediate matrices and fitted models are intentionally ignored by
Git; the checksum manifest and generation commands make their provenance auditable.
