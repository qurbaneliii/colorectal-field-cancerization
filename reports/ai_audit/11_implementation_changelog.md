# Scientific implementation changelog

- Reconstructed GEO metadata and the clinical join independently. Preserved the
  original GSE41258 scan-date values and documented the field. Corrected the
  covariate screen's missing-sentinel and single-tissue handling.
- Preserved the original external JSON and its metrics table. Generated the
  active JSON and Markdown from the archived current prediction evidence.
- Added an active external artifact verifier and invoked it before manuscript
  generation. Future external runs record input, model and result hashes.
  The general manifest now records inventory information.
- Bound the historical Task C threshold workflow to its development inputs and
  platform intersection. The rounded threshold grid now applies the recorded
  cutoff exactly. Labeled that threshold as conditional on the full-development
  panel and parameters.
- Mapped tumor probabilities by fitted class identity before saving or scoring
  external predictions. Log loss aligns columns with its label ordering.
- Implemented separate real-data experiments for Task B marker and stress
  ablations, a nested Task C rank-threshold policy, Task C archived score and
  platform diagnostics, and a locked GSE4107 transport sensitivity.
- Updated manuscript methods and limitations without replacing published
  numerical results. New experiments and their distinct estimands live in the
  audit package.
- Added regression tests for held-out patients at both threshold levels,
  exact cutoff behavior, probability order, stale metrics and threshold
  mismatch, scan-date preservation and portable text hashes.
