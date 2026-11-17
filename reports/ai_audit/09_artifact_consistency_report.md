# Scientific artifact consistency

## F04 | Historical external result reconciled

The 19-gene JSON, its matching 14-common-gene Table 5, and the old external
performance figure predate the current two-gene refit. Git history traces the
JSON and CSV to `62cce06` and the figure to `a5a53d7`. Their original Git
contents are preserved in `reports/archive/`. The active JSON and
Markdown were reconciled from the checked-in predictions and cards at
`c2a36a1`. The JSON records `historical_prediction_replay` so it can be
distinguished from a new model run.

`python scripts/verify_external_results.py` replays the stored score and class
mapping, checks patient/tissue identity and array counts against GEO metadata,
recomputes the metrics, checks the selected threshold and overlapping genes,
compares Tables 6 and 9, verifies evidence hashes and checks the rendered report.
CI also checks the primary figure's confusion counts and the manuscript's
reported external metrics against the same evidence. A fresh external run produces both report
formats from the same summary. `build_manuscript_outputs.py` checks the
contract before reading the external results.

## F07 | Inventory and generation provenance

The external run now records separate checksums for the expression matrices,
training and refitted model binaries, configuration, cohort metadata and
threshold evidence. It records the analysis source digest, run HEAD,
timestamp, Python version, platform and dependency-lock hash. Repository text
hashes normalize CRLF to LF, while matrices and models are hashed byte-for-byte.
The archival replay can verify the published score records and model-card
contents. It cannot recover the original serialized estimator or timestamp
from those records.

The general result manifest generator records the time and HEAD at inventory,
with no inferred generating command. The tracked pre-audit manifest remains a
historical snapshot, including rows for archived paths. Live generation
manifests for raw QC, DE, enrichment and every model family are outstanding.
Current checks cover the external prediction-to-report chain. Manuscript and
figure text still depend on the generation scripts and their input versions.
