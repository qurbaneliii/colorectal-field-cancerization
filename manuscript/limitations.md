# Limitations

This study is retrospective and uses two public microarray cohorts. Healthy and
cancer-bearing participants were not randomized or prospectively matched, so
unmeasured demographic, clinical, anatomic, medication, and processing factors
may remain. Recorded covariates were incomplete, particularly in GSE41258.

Bulk tissue conflates epithelial, stromal, endothelial, and immune signals.
MCP-counter is a marker-based sensitivity method, not a cell-count gold standard,
and composition-PC residualization may remove biologically meaningful field
signals as well as confounding. The immediate-early/stress audit is limited by a
curated gene list and incomplete preanalytical metadata.

Task B has internal grouped nested validation only. No independent cohort with
both genuinely healthy donor mucosa and tumor-adjacent mucosa was available, so
there is no external field-effect validation. Very high internal discrimination
may still reflect cohort-specific structure.

Task C uses a different normal comparator in GSE41258 and one of three primary
genes is absent on GPL96. Its two-gene transport refit is not the exact serialized
primary model. Strong ROC-AUC coexists with imperfect locked-threshold specificity
and calibration. Patient-cluster bootstrap accounts for repeated tissues but
does not correct dataset shift.

The models are retrospective biomarker-discovery tools. They have not undergone
assay validation, decision-curve analysis, prospective calibration, clinical
impact testing, regulatory review, or evaluation for individual patient care.
No causal biological mechanism or longitudinal progression is inferred from
the three cross-sectional tissue groups.
