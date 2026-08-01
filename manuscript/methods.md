# Methods

## Study design and datasets

This retrospective secondary analysis used public Affymetrix data. GSE44076 on
GPL13667 contained 246 audited arrays: 50 healthy-donor mucosa, 98
adjacent-normal mucosa, and 98 tumors from the same 98 cancer patients [2].
Healthy donors received unique donor IDs; paired tumor and adjacent samples
shared patient IDs. GSE41258 on GPL96 contained 390 deposited arrays [3]. The
external endpoint retained one canonical sample per patient and tissue after
excluding non-primary lesions and designated technical replicates, yielding
233 arrays from 190 patients (52
normal-colon and 181 primary-tumor
arrays); 43 patients contributed both tissues. No cross-platform
normalization or pooled fitting was performed.

## Raw-data audit, preprocessing, and quality control

CEL membership was reconciled to GEO accessions and platform CDFs before
analysis. GSE44076 was processed with the `oligo` route appropriate to HG-U219;
GSE41258 used the `affy` route appropriate to HG-U133A. Each cohort was
background corrected, quantile normalized, and summarized independently using
robust multi-array average [5]. Probe sets were mapped to gene symbols using
platform-specific annotation; ambiguous mappings and collapse rules are retained
in supplementary mapping tables.

Quality assessment comprised pre/post-RMA all-array density distributions,
boxplots, relative log expression, MA diagnostics, PCA, sample correlations,
hierarchical clustering, robust sample median/IQR diagnostics, and executed
arrayQualityMetrics reports [14]. For GPL96, probe-level models supplied NUSE;
no surrogate NUSE was generated for the oligo platform. A sample could be
excluded only after at least two independent severe custom technical failures;
PCA separation alone was never sufficient. Borderline arrays remained in the
primary analysis and were removed only in a labeled sensitivity analysis.

## Differential expression and field-gene evidence

Gene-wise linear models were fit with limma empirical Bayes moderation [6]. For
adjacent-normal versus healthy tissue, U0 included tissue only, U1 included
tissue, continuous age, and sex, and U2 additionally included left/right
location. Stage was not added because it is structurally undefined for healthy
donors. The paired tumor-versus-adjacent model used patient fixed effects and a
tissue coefficient. All design matrices and coefficient names were audited for
full rank. P values were adjusted by the Benjamini-Hochberg procedure [7]; the
main reporting threshold was FDR <0.05 and absolute log2 fold change >=0.5.

Tier 1 high-confidence field genes required strong U1 evidence (FDR <0.01 and
absolute adjusted log2 fold change >=1), consistent U0/U1 direction, U2 support,
at least 80% sample-direction consistency, deposited-matrix support, QC-exclusion
robustness, consistent direction in adequately represented age/sex/location
strata, and a field-compatible three-group trajectory. Tier 2 genes were labeled
provisional. Model coefficients did not define biological evidence tiers.

## Composition, stress, and enrichment sensitivity

MCP-counter 1.2.0 estimated eight immune and two stromal population scores from
human gene symbols [15]. The U3 sensitivity model added two composition principal
components learned from these scores. Genes retaining significance and direction
were called composition robust; attenuation was reported rather than interpreted
as proof of cell-intrinsic or cell-extrinsic causality. A 24-gene literature-based
immediate-early/stress set [18] was audited against measured symbols and removed
in a sensitivity enrichment analysis.

GO Biological Process and Reactome analyses used `clusterProfiler` [16] and the
Reactome knowledgebase [17]. Over-representation analyses used the tested mapped
gene universe, direction-specific gene sets, and BH correction. Ranked analyses
used signed limma statistics. Redundant GO terms were reduced by semantic
similarity; every output records the set, direction, database, universe, and
analysis type.

## Predictive modeling and internal validation

Task B (healthy versus adjacent-normal) was the primary predictive endpoint;
Task C (tumor versus adjacent-normal) was secondary; Task A was supporting.
Elastic Net logistic regression [8] was prespecified as the interpretable model,
with calibrated linear SVM and random forest comparators. All scaling, variance
filtering, univariate selection, residualization, hyperparameter tuning, feature
ranking, and compact-panel selection occurred within training partitions.

Internal performance used three deterministic repetitions of five-fold outer
cross-validation with donor/patient groups kept intact. Four-fold grouped inner
cross-validation selected hyperparameters and panel size from candidates 3, 5,
10, 15, and 20 plus the fold-specific nonzero model. The smallest panel within
the configured macro-F1 tolerance, balanced-accuracy floor, and Brier-loss limit
was selected. Outer-fold means and variability describe model-building-policy
performance; patient/donor bootstrap intervals of sample-aggregated repeated OOF
predictions describe a separate estimand. Calibration was assessed with curves,
Brier score, log loss, intercept, and slope where estimable [12].

Stability used all repeat-fold opportunities in its denominator. A gene was
strictly stable only with selection frequency >=0.65 and sign consistency >=0.80;
compact-consensus and exploratory labels were reported separately. Null tests
were task-specific: Task A held healthy labels fixed while permuting the
tumor/adjacent component, Task B permuted group-level healthy/adjacent labels,
and Task C exchanged tumor/adjacent labels within patients. Each used 1,000
permutations.

## Locked Task C threshold and external evaluation

The Task C panel and probability threshold were derived without GSE41258 labels.
Within-sample percentile ranks were computed over the label-independent common
gene universe. Repeated grouped GSE44076 OOF probabilities were averaged by
patient and the threshold maximizing balanced accuracy was selected, with ties
resolved toward 0.5 and then the lower threshold. Because FOXQ1 was
absent from GPL96, the exact three-gene serialized primary model was not directly
transportable. A distinct two-gene refit using CEMIP, ETV4 and
GSE44076 outcomes was serialized and labeled as a transport refit. External
labels were used only for final evaluation. The primary external estimand used
all canonical patient-tissue rows and patient-cluster bootstrap intervals (1,000
replicates); one-array-per-patient and threshold 0.5 were sensitivities.

## Reproducibility and reporting

Python and R dependencies were locked in `requirements-lock.txt`/`pyproject.toml`
and `renv.lock`, respectively. Seeds, grids, thresholds, panel sizes, bootstrap
counts, and permutation counts are centralized in `config/analysis.yaml`.
Reporting follows TRIPOD principles where applicable [13], while recognizing
that this is retrospective molecular biomarker discovery rather than a clinical
prediction-model validation study.
