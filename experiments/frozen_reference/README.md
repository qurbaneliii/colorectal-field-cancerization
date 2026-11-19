# Reproduce the experimental comparison

Use the preserved Python requirements lock and R 4.5.1 scientific library from the
primary raw-CEL workflow. Keep caches and generated intermediate molecular data
outside tracked paths. Set `R_LIBS_USER` to the verified library when needed.

From the repository root, with verified extracted CELs under `<cel_root>/<accession>`:

```text
Rscript --vanilla experiments/frozen_reference/extract_single_array.R <cel_root> <pm_cache>
Rscript --vanilla experiments/frozen_reference/check_raw_invariance.R <cel_root> <raw_invariance.csv>
python experiments/frozen_reference/run_comparison.py <pm_cache> <baseline_processed_dir> <hosted_bundle_root> <archive_cache> <new_output_dir>
python experiments/frozen_reference/verify_results.py <new_output_dir>
python experiments/frozen_reference/summarize_results.py <new_output_dir> <report_dir>
```

The baseline is the archived current-GEO joint-RMA matrix. The hosted bundle is
from run 43846818826 and supplies all original patient-grouped assignments.
The archive cache must contain both TARs matching the preserved current-GEO pins.
Per-CEL compressed digests must match the audited member inventories. Extraction
resumes only when its source-code and CEL digests match. Use a new output directory
for each model run; do not overwrite previous evidence.

Local execution used `phase-iii-pm-cache-v2` and generated
`phase-iii-measured-final` in the existing external audit cache. The earlier
`phase-iii-measured-b00e552` attempt stopped on external eligibility alignment
and is not a completed run. Raw input files,
processed matrices and numeric reference caches are not fetched automatically.
The measured execution records its analysis commit, input hashes and output hashes.
The committed result files preserve their measured bytes through `.gitattributes`.

This is fixed-panel conditional research. See the protocol and limitations before
interpreting a successful software run as evidence about clinical transport.
