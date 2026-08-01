# Task C model card

- Objective: secondary tumor versus adjacent-normal model
- Model family: Elastic Net logistic regression
- Training accession: GSE44076 only
- Expression provenance: `raw_cel_rma`
- Locked primary signature: FOXQ1, CEMIP, ETV4
- Per-gene stability status: `results/tables/final_compact_signature.csv`
- Threshold: selected later from GSE44076 repeated OOF rank-transport probabilities
- Exact primary artifact: `models/task_c_primary_full_signature_model.joblib`
- Cross-platform transport artifact: created separately during external validation
- Status: retrospective biomarker-discovery model; not clinically ready
