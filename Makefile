PYTHON ?= python
RSCRIPT ?= Rscript

.PHONY: setup-r audit preprocess processed-matrix qc differential-expression train permutation validate enrichment figures manuscript test sensitivity-all all

setup-r:
	$(RSCRIPT) R/00_install_packages.R

audit:
	$(PYTHON) scripts/run_data_audit.py

preprocess: setup-r
	$(RSCRIPT) R/01_extract_and_audit_cel.R
	$(RSCRIPT) R/02_preprocess_gse44076.R
	$(RSCRIPT) R/03_preprocess_gse41258.R
	$(PYTHON) scripts/convert_r_outputs.py

processed-matrix:
	$(PYTHON) scripts/run_processed_matrix_pipeline.py

qc:
	$(PYTHON) scripts/run_processed_matrix_pipeline.py --figures-only

differential-expression:
	$(RSCRIPT) R/04_differential_expression.R

train:
	$(PYTHON) scripts/run_modeling.py

permutation:
	$(PYTHON) scripts/run_permutation_tests.py

validate:
	$(PYTHON) scripts/run_external_validation.py

enrichment:
	$(RSCRIPT) R/05_functional_enrichment.R

figures:
	$(PYTHON) scripts/build_manuscript_outputs.py --figures

manuscript:
	$(PYTHON) scripts/build_manuscript_outputs.py

test:
	$(PYTHON) -m pytest -q

sensitivity-all: audit processed-matrix train permutation validate manuscript test

all: audit preprocess qc differential-expression train permutation validate enrichment manuscript test
