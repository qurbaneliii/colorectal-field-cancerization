# Limits of this experiment

The frozen method is a minimal quantile/median-PM comparator with no learned
probe-affinity or between-batch variance model. Changing summarization as well as
normalization prevents isolation of a pure normalization effect. Grid interpolation
is approximate. U219-derived target transfer to U133A is an experimental policy,
not validated platform calibration. Within-sample ranks reduce scale dependence
without establishing assay equivalence.

Panels are historical selections. Fixed-panel CV and conditional donor bootstrap
do not estimate the uncertainty of complete feature discovery, hyperparameter
selection or clinical deployment. The experiment does not rerun the stringent
Task B whole-gene sensitivity workflow, whose feature selection and residualization
require a separately frozen probe-level preprocessing graph at every nested fold.
That experiment remains NOT RUN. No empirical result will remove GSE44076's
case/acquisition aliasing.

Primary joint-RMA results, 0.41 threshold provenance, the historical 0.843 Task B
sensitivity and all historical TAR pins are preserved. GSE41258 is already
exposed retrospective data. There is no external-outcome optimization, new
independent cohort, clinical diagnosis claim or independent human approval.
Issue #7 remains open for established-method comparison, full-development nested
validation and a reviewed acceptance decision.
