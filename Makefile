PYTHON ?= python
RSCRIPT ?= Rscript

.PHONY: help restore-r verify-env audit extract preprocess convert processed-matrix qc differential-expression sensitivity train compact-panel permutation validate enrichment figures manuscript test smoke full-data-test sensitivity-all all

help:
	@echo "make all              Full raw-CEL publication pipeline"
	@echo "make sensitivity-all  Deposited-matrix sensitivity pipeline"
	@echo "make test             Data-independent unit tests"
	@echo "make full-data-test   Unit plus local full-data tests"

restore-r:
	$(RSCRIPT) -e "renv::restore(prompt = FALSE)"

verify-env:
	$(RSCRIPT) R/verify_environment.R
	$(PYTHON) -m pip check

audit:
	$(PYTHON) scripts/run_data_audit.py

extract:
	$(RSCRIPT) R/01_extract_and_audit_cel.R

preprocess:
	$(RSCRIPT) R/02_preprocess_gse44076.R
	$(RSCRIPT) R/03_preprocess_gse41258.R

convert:
	$(PYTHON) scripts/convert_r_outputs.py

processed-matrix:
	$(PYTHON) scripts/run_processed_matrix_pipeline.py

qc: preprocess

differential-expression:
	$(RSCRIPT) R/04_differential_expression.R

sensitivity:
	$(PYTHON) scripts/run_raw_vs_processed_sensitivity.py

train:
	$(PYTHON) scripts/run_modeling.py --provenance raw_cel_rma

compact-panel:
	$(PYTHON) scripts/run_compact_panel_analysis.py --provenance raw_cel_rma

permutation:
	$(PYTHON) scripts/run_permutation_tests.py --provenance raw_cel_rma

validate:
	$(PYTHON) scripts/run_external_validation.py --provenance raw_cel_rma

enrichment:
	$(RSCRIPT) R/05_functional_enrichment.R

figures:
	$(PYTHON) scripts/build_manuscript_outputs.py --figures

manuscript:
	$(PYTHON) scripts/build_manuscript_outputs.py

test:
	$(PYTHON) -m pytest -q -m "not full_data"

full-data-test:
	$(PYTHON) -m pytest -q

smoke:
	$(PYTHON) scripts/run_synthetic_smoke_test.py

sensitivity-all: audit processed-matrix
	$(PYTHON) scripts/run_modeling.py --provenance geo_deposited_series_matrix
	$(PYTHON) scripts/run_compact_panel_analysis.py --provenance geo_deposited_series_matrix
	$(PYTHON) scripts/run_permutation_tests.py --provenance geo_deposited_series_matrix
	$(PYTHON) scripts/run_external_validation.py --provenance geo_deposited_series_matrix
	$(PYTHON) scripts/build_manuscript_outputs.py
	$(PYTHON) -m pytest -q

all: audit extract preprocess convert differential-expression sensitivity enrichment train compact-panel permutation validate manuscript full-data-test
