# Test execution report

## Publication environment

- Python: 3.12.10 in `.venv-publication`
- Dependency integrity: `python -m pip check` passed
- Direct imports: all locked runtime and test dependencies passed
- Static compilation: `python -m compileall -q src scripts tests` passed
- Lint: `ruff check src scripts tests` passed
- Unit/no-data suite: 6 passed, 17 deselected
- Synthetic grouped-CV smoke test: PASS; 96 predictions, 8 metric rows,
  192 assignment rows, and 40 coefficient rows
- Full suite after all generated artifacts: **23 passed, 0 failed** in 13.10 s
- R: 4.5.1 / Bioconductor 3.21; `R/verify_environment.R` passed

## Full-data execution

- CEL archive audit: PASS (246 GSE44076; 390 GSE41258)
- GSE44076 raw RMA/QC: PASS (49,386 probe sets)
- GSE41258 raw RMA/QC: PASS (22,283 probe sets)
- limma and enrichment: PASS
- grouped repeated nested CV: PASS (135 fold-metric rows; 5,310 held-out predictions)
- compact panels and model reload: PASS
- 1,000-iteration permutation tests: PASS for all three tasks
- external validation and 1,000-iteration grouped bootstrap: PASS, with
  calibration explicitly NA when separation prevented estimation
- manuscript/artifact build: PASS (63 PNG figures, 43 CSV tables, 309 manifest
  rows at the recorded build point)

Docker build/runtime verification is BLOCKED by the disabled host WSL service;
see `reports/docker_verification.md`.
