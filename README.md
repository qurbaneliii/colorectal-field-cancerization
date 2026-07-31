# Transcriptomic field cancerization in colorectal cancer

Reproducible hybrid R/Python analysis for **GSE44076** (healthy colon,
tumor-adjacent normal colon, and paired colorectal tumor) with **GSE41258**
used only for the scientifically valid tumor-versus-normal-colon external
validation.

The project deliberately separates:

- raw-CEL preprocessing and biological inference (R/Bioconductor);
- leakage-safe predictive modeling and reporting (Python);
- primary field-effect inference from external tumor/normal validation.

No random sample-level train/test split, SMOTE, or deep learning is used.
Every predictive split is grouped by patient/donor, and every learned
preprocessing operation is fitted inside its training fold.

## Quick start

```powershell
python -m venv --system-site-packages .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
.\.venv\Scripts\python scripts/run_data_audit.py
.\.venv\Scripts\python scripts/run_processed_matrix_pipeline.py
.\.venv\Scripts\python scripts/run_modeling.py
.\.venv\Scripts\python scripts/run_external_validation.py
.\.venv\Scripts\python scripts/build_manuscript_outputs.py
.\.venv\Scripts\python -m pytest -q
```

For publication-grade raw-CEL processing, install R and run:

```powershell
Rscript R/00_install_packages.R
Rscript R/01_extract_and_audit_cel.R
Rscript R/02_preprocess_gse44076.R
Rscript R/03_preprocess_gse41258.R
Rscript R/04_differential_expression.R
Rscript R/05_functional_enrichment.R
```

The committed `renv.lock` is a bootstrap lock for R 4.5.1/Bioconductor 3.21.
Because R was unavailable during this execution, its package section is not
represented as resolved. `R/00_install_packages.R` installs the declared
packages and replaces it with the fully resolved snapshot; commit that resolved
lock before a final manuscript freeze.

`make all` runs the complete raw-CEL R/Python route on systems with R and GNU
Make. `make sensitivity-all` reproduces the deposited-series-matrix route
executed in the current R-free environment.
Configuration lives in `config/analysis.yaml` and `config/paths.yaml`.

## Scientific boundary

The Python processed-matrix route consumes the normalized expression deposited
in GEO series matrices. It supports metadata audit, exploratory analysis,
leakage-safe model development, and external validation. It is **not** a
substitute for executing raw-CEL RMA/QC in Bioconductor. Reports label the
provenance of each result explicitly.

Raw data and generated matrices are ignored by Git. See
`reports/final_analysis_report.md` for the executed-state audit and limitations.
