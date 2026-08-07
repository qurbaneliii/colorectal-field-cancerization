# Test execution report

## Publication environment

- Python: 3.12.10 in `.venv-publication`
- Dependency integrity: `.venv-publication\\Scripts\\python.exe -m pip check`
  passed with no broken requirements.
- Static compilation: `.venv-publication\\Scripts\\python.exe -m compileall -q
  src scripts tests` passed.
- Lint: `.venv-publication\\Scripts\\ruff.exe check src scripts tests` passed.
- Baseline full suite before this revision: 37 passed, 0 failed.
- Final unit/no-data suite after contract expansion: 12 passed, 28 deselected.
- Synthetic grouped-CV smoke test: PASS; 96 predictions, 8 metric rows,
  192 assignment rows, and 40 coefficient rows.
- Full suite after all generated artifacts: **40 passed, 0 failed**.
- R: 4.5.1 / Bioconductor 3.21; `R/verify_environment.R` passed.
- `renv::status()` passed: `No issues found -- the project is in a
  consistent state.`

## Full-data execution

- CEL archive audit: PASS (246 GSE44076; 390 GSE41258)
- GSE44076 raw RMA/QC: PASS (49,386 probe sets)
- GSE41258 raw RMA/QC: PASS (22,283 probe sets)
- U0/U1/U2 and paired limma artifacts: PASS; current code/results consistency
  validated without repeating unchanged raw RMA or limma fits.
- Enrichment: PASS; the ranked-only route was re-executed after adding explicit
  core leading-edge membership (2,517 ranked terms with nonempty gene IDs).
- grouped repeated nested CV: PASS (135 fold-metric rows; 5,310 held-out predictions)
- compact panels and model reload: PASS
- 1,000-iteration permutation tests: PASS for all three tasks. The parent shell
  timed out after the workstation slept, but the preserved workers completed,
  exited, and atomically wrote the fresh 1,000-iteration result rows at
  2026-08-01 12:49 local time; the final contract tests independently validated
  iteration counts and null-hypothesis names.
- external validation and 1,000-iteration grouped bootstrap: PASS, with
  undefined predictive-value bootstrap summaries explicitly NA when no valid
  draw was estimable.
- publication figure build: PASS (12 generated figure stems, 36 PNG/PDF/SVG
  renderings in the final builder pass; the full results directory also retains
  required upstream QC and modeling figures).
- publication table build: PASS (10 principal journal tables).
- manuscript/artifact build: PASS (18 manuscript sections; 568 nonempty files
  indexed in the final artifact manifest).

## Final verification commands

```powershell
.\\.venv-publication\\Scripts\\python.exe -m pip check
.\\.venv-publication\\Scripts\\python.exe -m compileall -q src scripts tests
.\\.venv-publication\\Scripts\\ruff.exe check src scripts tests
.\\.venv-publication\\Scripts\\python.exe scripts\\run_synthetic_smoke_test.py
.\\.venv-publication\\Scripts\\python.exe -m pytest -q
& 'D:\\R\\R-4.5.1\\bin\\Rscript.exe' R\\verify_environment.R
& 'D:\\R\\R-4.5.1\\bin\\Rscript.exe' -e 'renv::status()'
```

## GitHub Actions

- Baseline branch CI: **PASS**. Workflow run #9 completed successfully for
  baseline commit `70998d6ffb537ddf720adb46ff217c373edce9fb`:
  <https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/30693185908>
- That successful job executed dependency installation, `pip check`,
  `compileall`, Ruff, the no-data test selection, synthetic smoke, and
  configuration/documentation parsing.
- CI for the current unpushed revision is **BLOCKED** because GitHub CLI is not
  authenticated in this host session. No current-revision remote PASS is claimed.
- The manual `full-data.yml` workflow was **NOT RUN** remotely. It requires a
  self-hosted runner labeled `colorectal-field-data`; required runner and data
  availability were not verified. The equivalent full-data suite was executed
  locally and passed as recorded above.

Docker build/runtime verification is BLOCKED by the disabled host WSL service;
see `reports/docker_verification.md`.
