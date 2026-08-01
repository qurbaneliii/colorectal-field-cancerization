# Transcriptomic field cancerization in colorectal cancer

Publication-oriented, leakage-safe reanalysis of **GSE44076** with raw Affymetrix CEL
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
- Primary Task B compact panel: **CLC, DYNC1H1, FOS, VIP, SNORA12**. It is internally
  validated only; no independent healthy-versus-adjacent cohort is available.
- Secondary Task C compact panel: **FOXQ1, CEMIP, ETV4**. All three meet the strict
  stability rule; only **CEMIP** and **ETV4** are present on GPL96.
- External Task C primary analysis: 233 canonical patient-tissue arrays from 190
  patients; patient-clustered ROC AUC 0.9829 and macro-F1 0.8282 at the
  GSE44076-locked threshold of 0.41. The two-gene transport refit is distinct from the
  exact three-gene primary model.
- Final readiness: **PARTIAL** until the Docker image is built against a running Docker
  engine; all scientific gates completed locally.

See [the final readiness report](reports/publication_readiness_final.md),
[the manuscript](manuscript/manuscript.md), and
[the artifact manifest](reports/result_artifact_manifest.csv).

## Reproduce

The supported interpreter versions are Python 3.12 and R 4.5.x. On Windows, invoke the
installed R executable explicitly if `Rscript` is not on `PATH`.

```powershell
py -3.12 -m venv .venv-publication
.\.venv-publication\Scripts\python -m pip install -r requirements-lock.txt
.\.venv-publication\Scripts\python -m pip install --no-deps -e .

D:\R\R-4.5.1\bin\Rscript.exe R/00_install_packages.R
D:\R\R-4.5.1\bin\Rscript.exe R/verify_environment.R
D:\R\R-4.5.1\bin\Rscript.exe R/01_extract_and_audit_cel.R
D:\R\R-4.5.1\bin\Rscript.exe R/02_preprocess_gse44076.R
D:\R\R-4.5.1\bin\Rscript.exe R/03_preprocess_gse41258.R
.\.venv-publication\Scripts\python scripts/convert_r_outputs.py
D:\R\R-4.5.1\bin\Rscript.exe R/04_differential_expression.R
.\.venv-publication\Scripts\python scripts/run_raw_vs_processed_sensitivity.py
.\.venv-publication\Scripts\python scripts/run_covariate_sensitivity.py
D:\R\R-4.5.1\bin\Rscript.exe R/06_tissue_composition_sensitivity.R
.\.venv-publication\Scripts\python scripts/run_stress_gene_sensitivity.py
D:\R\R-4.5.1\bin\Rscript.exe R/05_functional_enrichment.R
.\.venv-publication\Scripts\python scripts/run_modeling.py --provenance raw_cel_rma
.\.venv-publication\Scripts\python scripts/run_compact_panel_analysis.py --provenance raw_cel_rma
.\.venv-publication\Scripts\python scripts/run_task_b_confounding_sensitivity.py
.\.venv-publication\Scripts\python scripts/run_permutation_tests.py --provenance raw_cel_rma
.\.venv-publication\Scripts\python scripts/select_task_c_threshold.py
.\.venv-publication\Scripts\python scripts/run_external_validation.py --provenance raw_cel_rma
.\.venv-publication\Scripts\python scripts/build_publication_figures.py
.\.venv-publication\Scripts\python scripts/build_publication_tables.py
.\.venv-publication\Scripts\python scripts/build_manuscript_outputs.py
.\.venv-publication\Scripts\python -m pip check
.\.venv-publication\Scripts\ruff check src scripts tests
.\.venv-publication\Scripts\python -m pytest -q
```

`make all` expresses the same full route on systems with GNU Make. `make smoke` runs a
no-data grouped-CV smoke test. `make sensitivity-all` reproduces the explicitly labeled
GEO-deposited-series-matrix sensitivity route.

The executed QC bundle includes cohort-wide pre/post distributions, RLE, MA, PCA,
sample correlation, hierarchical clustering, objective multimetric exclusion flags,
and `arrayQualityMetrics` HTML reports for both cohorts. GPL96 also has affyPLM NUSE;
no surrogate NUSE is used for the oligo platform.

## Provenance and safeguards

Raw inputs under `data/raw/` are never modified. Their paths, sizes, SHA-256 hashes,
archive-member counts, GEO mappings, and platform checks are recorded under
`data/metadata/`. Raw-CEL and deposited-matrix artifacts use separate filenames and
result namespaces.

All model splits preserve patient/donor groups. Variance filtering, univariate selection,
scaling, tuning, calibration, panel ranking, and threshold use training data only. No
SMOTE, deep learning, external-label tuning, or silent mock fallback is used.

The locked environments are `requirements-lock.txt` and `renv.lock`. CI runs lint,
no-data tests, and a synthetic grouped-CV smoke test; the manual full-data workflow is
provided for a self-hosted runner with the local GEO inputs.

## Data availability

The raw and deposited expression inputs are public through NCBI GEO (GSE44076 and
GSE41258). Large raw/intermediate matrices and fitted models are intentionally ignored by
Git; the checksum manifest and generation commands make their provenance auditable.
