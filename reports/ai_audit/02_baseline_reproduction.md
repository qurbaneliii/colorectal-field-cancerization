# Baseline reproduction

Verified upstream baseline `c2a36a1d2f80fd25e1565b543c909eebef64ca09`.
Baseline commands used a detached worktree at that commit. The audit branch
kept the deposited, archived-prediction and raw-CEL evidence distinct.
Execution took place on 11 November 2026 on Windows x64 with four logical
processors and 7.75 GiB physical RAM. Python 3.12.14 ran from an isolated
virtual environment with `requirements-lock.txt`. R 4.5.1 ran from an
audit-local installation. The Docker engine was unavailable.

## Verified inputs

The original raw archives, two deposited matrices and GSE41258 clinical
supplement are published in GEO. The tables below distinguish expected and
observed SHA-256 values. The full compressed raw inputs require about 1.91 GB
before extraction, QC or model intermediates. Raw downloads were staged in
the local audit cache without overwriting any repository input.

| Input | Bytes expected | Verification |
|---|---:|---|
| GSE44076 series matrix | 32,596,137 | Exact match to tracked SHA-256 |
| GSE41258 series matrix | 17,537,584 | Exact match to tracked SHA-256 |
| GSE41258 clinical supplement | 3,135 | Exact match to tracked SHA-256 |
| GSE44076 raw archive | 532,213,760 | Full archive downloaded. All 246 member names matched the inventory and all HG-U219 CEL gzip streams passed CRC checks. The archive digest differs from the historical manifest. |
| GSE41258 raw archive | 1,378,918,400 | Full verified input unavailable after download timeouts. The public endpoint and reported size were confirmed. |

The GSE44076 observed SHA-256 was
`059abdaa76f731c9d1af116875bb60663fb875c296bbf01e353de48eeed0a030`.
The historical manifest expects
`8db70acb37d0fc85e85d090459312c8e97960f3cdcd19d1998ee03fb6be5f522`.
Matching member names and readable CEL files establish archive integrity.
They do not establish byte identity with the baseline input. The current
archive was preserved outside the primary route.
[GSE44076 raw archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44076/suppl/GSE44076_RAW.tar),
[GSE41258 raw archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41258/suppl/GSE41258_RAW.tar)

## Execution evidence

| Command | Result | Duration and exit |
|---|---|---|
| `python -m pip install -r requirements-lock.txt` and `python -m pip check` | Locked Python dependencies installed. Dependency check passed. | Exit 0 |
| `python -m pytest -q -m "not full_data"` | Baseline 12 passed. Revised 22 passed. | Exit 0 |
| `python -m pytest -q -m full_data --tb=no` | 15 passed, 1 skipped, 12 failed on prerequisites including absent raw archives, processed matrices, fitted models or generated fold/prediction artifacts. | 21.12 s, exit 1 |
| `python scripts/run_synthetic_smoke_test.py` | Grouped synthetic execution produced 96 predictions. | Exit 0 |
| `Rscript --vanilla -e "renv::restore(prompt = FALSE)"` | Downloaded 204 packages. Restore failed when pinned source packages required the Windows build toolchain. Missing `sh` and compiler tools were reported. | Exit 1 |
| R parsing of `R/*.R` | All nine baseline R scripts parsed. | Exit 0 |
| `python scripts/run_processed_matrix_pipeline.py` at baseline | Deposited exploratory route built 19,040 primary and 13,299 external genes, with 11,336 shared. It emitted 1,654 exploratory field candidates using the baseline GEO mapping. | 242.60 s, exit 0 |
| `Rscript --vanilla scripts/audit_deposited_limma.R <mapped gene matrix CSV> results/ai_audit` | Four limma designs ran with the cached locked limma 3.64.3 and statmod 1.5.2. Input was GSE44076 deposited expression mapped using the pinned primary annotation and the repository median aggregation. | 12.69 s for the R fit, exit 0 |
| `python scripts/audit_geo_metadata.py --gse44076 <brief> --gse41258 <brief> --gse41258-clinical <supplement>` | Reconciled GEO fields, groups, exclusions and clinical joins. | Exit 0 |
| `python scripts/verify_external_results.py` | Recomputed locked external metrics from published scores, checked reference metadata and reconciled the active report. | Exit 0 |
| `python scripts/run_development_audit.py --series-matrix <GSE44076 matrix>` | Executed the real-data Task B ablations and nested Task C policy. Uses pinned primary-route annotation over the deposited matrix. | 756.42 s, exit 0 |
| `python scripts/audit_external_transport.py --bootstrap-replay --primary-matrix <GSE44076 matrix> --external-matrix <GSE41258 matrix>` | Replayed saved scores, reproduced 1,000 patient-cluster bootstrap draws for the canonical rank evaluation and measured deposited-expression drift. | 43.89 s, exit 0 |
| `python scripts/run_task_b_transfer_sensitivity.py --primary-matrix <GSE44076 matrix> --external-matrix <GSE4107 matrix> --external-annotation <GPL570 annotation> --external-geo-text <GSE4107 brief>` | Executed locked early-onset cohort sensitivity. | 64.61 s, exit 0 |
| `python scripts/select_task_c_threshold.py` and `python scripts/run_external_validation.py --provenance geo_deposited_series_matrix` in a scratch worktree | Exercised the revised lock, model and report generation using verified deposited matrices and the archived feature specification. The scratch config used 10 bootstrap replicates for execution verification. | 47.93 s including input preparation, exit 0 |

The archived raw-CEL QC, original raw-RMA limma fit, composition, enrichment,
full nested modeling, permutations and full external refit were blocked.
The canonical external score replay and the separate real-data sensitivities
provide narrower evidence. The synthetic smoke test exercises software only.
Evidence hashes, splits, seeds, thresholds and outcome tables are in
`results/ai_audit/`.
