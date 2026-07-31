# Publication-readiness baseline

Audit date: 2026-07-31 (Asia/Baku)  
Baseline commit: `62cce06e252cfc01582ab8bc1664d8c1cc037834`  
Working branch: `feat/publication-ready-field-cancerization`  
Starting branch: `main`  
Working tree before the audit report was written: clean

## Repository tree and implementation inventory

The repository contains 93 tracked files and preserves a hybrid R/Python
implementation. The relevant tree at baseline was:

```text
.
|-- config/                 analysis and path YAML
|-- data/
|   |-- raw/                five ignored GEO inputs
|   |-- interim/            ignored sklearn cache; no extracted CEL tree
|   |-- processed/          six ignored deposited-matrix artifacts
|   `-- metadata/           committed manifests, sample audits, mappings
|-- manuscript/             methods, results, limitations, figure legends
|-- models/                 two ignored Task C joblib artifacts
|-- R/                      install, extraction, RMA, limma, enrichment scripts
|-- reports/                six committed baseline reports
|-- results/
|   |-- differential_expression/  four ignored processed-matrix tables
|   |-- figures/            45 ignored PNG/PDF/SVG files
|   |-- metrics/            nine ignored modeling/validation tables
|   `-- tables/             20 committed tables
|-- scripts/                seven Python entry points
|-- src/                    data, biology, modeling, interpretation, plotting
`-- tests/                  15 collected tests across six files
```

The repository has one commit. No existing changes were present, and no valid
work was removed or overwritten during this audit.

## Input data

All five expected files are available locally under ignored `data/raw/` paths.
The raw files were read for audit only and were not modified.

| Accession | File | Size (bytes) | SHA-256 | Archive members |
|---|---|---:|---|---:|
| GSE44076 | `GSE44076_RAW.tar` | 532,213,760 | `8db70acb37d0fc85e85d090459312c8e97960f3cdcd19d1998ee03fb6be5f522` | 246 |
| GSE44076 | `GSE44076_series_matrix.txt.gz` | 32,596,137 | `096db6c83516d01e39a8c1e3b84b2b3759d9edb133e0ce40807251a69b2671bb` | not applicable |
| GSE41258 | `GSE41258_RAW.tar` | 1,378,918,400 | `965cc38ccf8cf5f803de2bf3e15da0232fd629196cd5c9e75f415a106902576d` | 390 |
| GSE41258 | `GSE41258_series_matrix.txt.gz` | 17,537,584 | `95931c3bb5fef20e4e056e4e3b1c354aa9272a7e1c60e1fd54e2849b30bb582a` | not applicable |
| GSE41258 | `GSE41258_clinical_data.txt.gz` | 3,135 | `6ed90de8810b609d175c34b43443a64b6f973c9bf8aff85b7e666e16cbdb8194` | not applicable |

Both TAR member-name sets are unique. All 246 and 390 GEO sample-to-CEL names
match archive members. A streamed header inspection identified `HG-U219` in
the first GSE44076 Calvin CEL and `HG-U133A` in the first GSE41258 CEL. This is
preliminary evidence only; every CEL and the platform/CDF compatibility still
require the corrected integration-stage audit.

## Cohort audit

- GSE44076: 246 included arrays; 50 healthy, 98 adjacent-normal, 98 tumor;
  148 unique donor/patient groups; 98 complete tumor-adjacent pairs; no
  duplicate pair entries.
- GSE41258: 233 eligible arrays; 181 primary tumors and 52 normal colon arrays;
  190 unique included patients; 43 patients with both tissues, 138 tumor-only,
  9 normal-only; 43 patients have more than one eligible array; 8 technical
  replicate candidates are excluded.

## Available processed matrices and local artifacts

The following local matrices exist but are ignored by Git and were generated
from GEO deposited series matrices, not raw-CEL RMA:

| Artifact | Rows | Samples | Observed range | Baseline interpretation |
|---|---:|---:|---:|---|
| `GSE44076_probe_expression.parquet` | 49,386 | 246 | not recalculated here | deposited matrix |
| `GSE44076_gene_expression.parquet` | 19,040 | 246 | 1.4418 to 13.57575 | plausibly log2-scale |
| `GSE41258_probe_expression.parquet` | 22,283 | 390 | not recalculated here | deposited matrix |
| `GSE41258_gene_expression.parquet` | 13,299 | 390 | 0.351 to 528,000 | not log2-scale |

There are 11,336 common gene symbols. The GSE41258 scale disproves the current
documentation's blanket description of deposited values as normalized log2
expression and requires explicit transformation/provenance handling.

Two ignored model artifacts are present:

- `models/task_c_locked_elastic_net.joblib`: 19-gene Task C signature, 196
  GSE44076 training arrays, deposited-matrix provenance.
- `models/task_c_locked_rank_elastic_net.joblib`: 14 common genes and a
  GSE44076-only fitted rank model.

Ignored local metrics include repeated nested-CV predictions and metrics,
fold assignments, Elastic Net coefficients, bootstrap tables, permutation
tests, and external predictions. Committed compact summaries include Tables
1-5, feature stability, signatures, exploratory differential-expression
summaries, metadata, PCA scores, and processed-matrix field candidates.

## Baseline test execution

Exact command:

```powershell
.\.venv\Scripts\python -m pytest -q
```

| Item | Result |
|---|---:|
| Collected | 15 |
| Passed | 15 |
| Failed | 0 |
| Skipped | 0 |
| Warnings | 0 |
| Exit code | 0 |
| Runtime | 44.18 seconds |

These tests do not cover several mandatory publication gates: real CEL
platform mismatch failures, inner-fold leakage, configuration-grid equality,
correct repeat-aware stability denominators, compact-panel nesting,
patient-balanced external evaluation, threshold provenance, or a synthetic
end-to-end pipeline.

## Runtime environments

- Python: 3.12.10 in `.venv`; pip 25.0.1; scikit-learn 1.8.0.
- `python -m pip check`: exit 1. Conflict: `opentelemetry-proto 1.41.1`
  requires `protobuf>=5,<7`, while the environment contains `protobuf 4.25.9`.
- Rscript: unavailable on PATH; exit 1.
- GNU Make: unavailable on PATH; exit 1.
- Docker client: present, but the Docker Desktop Linux engine is unavailable;
  `docker info` exits 1.
- Disk: D: has approximately 939 GB free; C: has approximately 2.27 GB free.
- `renv.lock`: bootstrap only, with an empty `Packages` object; not a resolved
  publication environment.

## Completed stages at baseline

- Raw input discovery, checksum capture, metadata parsing, cohort counts, and
  sample-to-archive-name matching.
- Deposited-series-matrix ingestion, label-independent probe aggregation,
  exploratory Welch/paired-t sensitivity analysis, and exploratory figures.
- Two-repeat grouped nested CV for three tasks, comparator fitting, stored OOF
  predictions, preliminary feature stability, model serialization, and
  array-level GSE41258 evaluation.
- A limited unit-test suite and partial manuscript/report scaffolding.

## Blocked or incomplete stages at baseline

- Full CEL extraction and per-file integrity/readability verification.
- Genuine all-array chip/platform/CDF validation.
- Raw-CEL RMA, full raw/normalized QC, arrayQualityMetrics, RLE/NUSE, objective
  exclusion evidence, and provenance-separated raw matrices.
- Raw-CEL limma, paired-design validation, raw-versus-deposited concordance,
  final field candidates, enrichment, and pathway figures.
- Leakage-safe compact-panel evaluation and final raw-derived model locking.
- Patient-balanced external primary analysis and all-array cluster bootstrap.
- Calibration intercept/slope, correctly labeled Task A permutation null,
  sufficient final permutation resolution, and required sensitivity analyses.
- Resolved R/Python locks, restorable Docker build, CI, artifact manifest,
  complete manuscript, and all required publication tables/figures/reports.

## Baseline methodological risks

1. `assert_platform()` only prints an expectation and cannot reject a mismatch.
2. GSE41258 forces `cdfname = "hgu133acdf"` without demonstrating the detected
   CEL annotation/CDF environment and probe-set compatibility.
3. R output filenames can collide semantically with deposited-matrix products;
   provenance is not enforced in matrix filenames or sidecars.
4. The limma script builds an unpaired coefficient name by string
   concatenation and does not emit the required complete schema or provenance.
5. YAML and Python contain different hyperparameter grids and only two outer
   repeats; compact panel sizes omit 5 and 15 genes.
6. LinearSVC decision scores are sigmoid/softmax transformed and reported as
   probabilities, invalidating its log-loss/Brier/calibration comparisons.
7. Feature stability counts coefficient rows using only `outer_fold` size and
   does not explicitly count unique `(repeat, outer_fold)` selections.
8. Candidate compact panels are global lists and have no untouched outer-fold
   performance estimate.
9. GSE41258's 233 arrays represent 190 patients; current point estimates
   overweight patients with both tissues and do not implement a canonical
   patient-level primary set.
10. The within-sample-rank representation is described as prespecified only
    after both variants were evaluated; the final protocol must lock its
    scientific rationale independently of external labels.
11. Current Task A permutation preserves/conditionally alters label structure
    but is not labeled as a conditional null test.
12. The committed reports overstate completeness relative to the absent raw
    pipeline, unresolved environments, missing CI, and incomplete manuscript.

## Baseline readiness decision

**BLOCKED for publication readiness.** The existing processed-matrix work is a
valuable, explicitly exploratory sensitivity analysis. The raw-CEL biological
inference gates and several predictive-validity/reproducibility gates have not
been executed or satisfied in this environment.
