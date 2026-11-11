# Raw-CEL reproduction continuation

This audit continues the analysis at `0ede337ffc360bc96dea6724f974bd5ef7f184d4`
on PR #4. The historical manifest remains an immutable record of the earlier
analysis. The fresh run uses a separate worktree and records actual inputs.

## Completion status

| Check | Status | Evidence |
| --- | --- | --- |
| Current GEO acquisition, metadata and CEL integrity | PASS | Both current archives and all 636 CELs verified. GSE41258 resumed from a retained partial file. |
| Historical TAR SHA-256 equality | FAIL | Both current digests differ from the historical manifest. |
| Cause of the historical differences | BLOCKED | Requires the historical TARs or historical per-member fingerprints. |
| Pipeline on the original historical bytes | NOT RUN | The original pinned TARs are unavailable. |
| Windows R 4.5.1 / Bioconductor 3.21 | PASS | 219 exact package-version matches. The pipeline used gzip CSV with Python Parquet conversion. |
| Current GEO raw RMA, core QC and GPL96 NUSE | PASS | Both cohorts and all 636 arrays. |
| limma, covariate, composition and stress sensitivities | PASS | Fresh fits and archived comparisons. |
| Nested grouped CV, biomarker panels, permutations and bootstraps | PASS | Fresh model fits, three 1,000-permutation tests and 1,000-draw cluster intervals. |
| Task C threshold and external validation | PASS | Fresh nested selection, GSE44076 lock and GSE41258 evaluation reproduced the archived primary result. |
| Acquisition-matched Task B validation | BLOCKED | JGAD000785 requires approved access and cohort eligibility review. |

[Machine-readable stage status](../../results/ai_audit/current_geo_stage_status.csv)

## Archive identity

| Accession | Historical SHA-256 | Current GEO SHA-256 | Verified members |
|---|---|---|---:|
| GSE44076 | `8db70acb37d0fc85e85d090459312c8e97960f3cdcd19d1998ee03fb6be5f522` | `059abdaa76f731c9d1af116875bb60663fb875c296bbf01e353de48eeed0a030` | 246 |
| GSE41258 | `965cc38ccf8cf5f803de2bf3e15da0232fd629196cd5c9e75f415a106902576d` | `740aaf0c5ec4c0ba72107b1d4ef50b5a309cb5cab88aa29ea73fc48594f778ee` | 390 |

Both downloads match the historical byte counts and ordered CEL inventories.
Every gzip stream was read through its CRC and size trailer, and every CEL
header matched the expected platform. The member manifests record ordered
tar header, gzip and decompressed CEL SHA-256 fingerprints. Eight evenly
spaced members per archive also matched the corresponding files served by
the official GEO sample endpoints, byte-for-byte in both representations.
This corroborates the current source. [GSE44076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44076),
[GSE41258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41258)

The four committed manifest revisions `62cce06`, `294d5ce`, `20d344b` and
`8452768` retain the same historical digests. They provide member names and
order, with no historical per-member digests. The original tar files were
unavailable in the checked local repositories and drives. The difference is
therefore confirmed at the tar level. Repacking, gzip header changes and
changes to CEL payloads cannot be distinguished for the historical inputs
until an original tar or its member fingerprints are recovered. The
`--historical-archive` comparison implements that distinction and requires
the original archive to match its historical pin.

## Acquisition and staging

GSE41258 resumed at byte 569,329,071 and reached 1,378,918,400 bytes. The
server returned a matching HTTP Content-Range and the completed archive
passed all CEL checks. The audited archives were promoted to the verified
local cache under explicit current-GEO pins. The supporting matrices and
clinical supplement matched their historical digests. The staging manifest
records all five observed hashes alongside their historical pins.
An independent Linux runner also downloaded and verified the same five
pinned inputs. This confirms the current archive digests independently of
the local resumed transfers. [Hosted reproduction run](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/43008659539)

For a current-GEO run in a separate checkout, use:

```text
python scripts/stage_geo_inputs.py --download --cache <cache> --worktree <worktree> --alternate-pins results/ai_audit/current_geo_archive_pins.json --evidence <input-manifest.json>
```

The downloader checks the promised range and available remote validators
before appending. It checks the entire file digest before promotion and
checks the cache before opening a network connection. Omitting the alternate
pin file requests the historical input identities.

The metadata audit ran on the verified current archives and retained the
published cohorts and eligibility. It reported zero count discrepancies.
All six raw-input tests in `tests/test_metadata.py` passed.

The R extraction ran in 41.0 s, and all 636 extracted gzip files matched the
archive member hashes. Native GSE44076 RMA and the intensity, RLE, PCA and
clustering diagnostics ran in 525.2 s. The native call deferred the separate
`arrayQualityMetrics` report while its dependency closure was being restored.
GSE41258 RMA, core array diagnostics and the original GPL96 affyPLM NUSE
stage completed in 3,049.9 s. The 390 NUSE medians and IQRs agreed with the
archive within 1.0 × 10⁻¹⁴. Core QC retained all 636 arrays and recovered
GSM1077736 as the sole one-metric review case. The fresh probe-set counts
were 49,386 and 22,283 and the mapped gene counts were 18,490 and 12,650.

The initial Linux restore built 176 packages and failed at `fs` (missing
`libuv1-dev`), `oligo` (the Bioconductor source endpoint did not serve the
locked 1.71.7), and `MCPcounter` (a duplicate staging path). The hosted job
now installs the system library and builds the two source tarballs from the
revisions recorded in the lock. The Windows host also installed Rtools45
and compiled and linked an R test DLL with GCC 14.3.0. The Windows scientific
run used the locked core packages and the supported gzip CSV/Parquet route.
[oligo 1.71.7 source](https://github.com/bioc/oligo/commit/e073acb9f29886d4cc073f634c109d0bb76b2f81),
[MCPcounter source](https://github.com/ebecht/MCPcounter/tree/b6eac73e91c246fcff0bb1a5c68a816cd588fc48)

The second hosted run installed 201 locked packages, including `oligo`
1.71.7 and both QC and enrichment dependencies. The duplicate renv staging
path persisted for `MCPcounter`. The hosted workflow now installs that
package directly from its pinned source tarball into the locked project
library. The verified GEO input cache was saved before this restore.
[Second hosted run](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/43013057519)

The next hosted run reused the verified GEO cache and passed locked R and
Python dependency verification. Its raw stage stopped at
`pthread_create` error 22 from the threaded `preprocessCore` build.
The hosted recipe now rebuilds the same locked `preprocessCore` version
with the package's serial configure option and checks RMA background
correction and quantile normalization before the scientific stages.
The first serial attempt still linked the threaded cache entry. The
current recipe installs the locked Bioconductor source with the configure
option explicitly on `R CMD INSTALL`.
[Hosted restore and RMA attempt](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/43017039168),
[Bioconductor threading diagnosis](https://github.com/Bioconductor/bioconductor_docker/issues/22)

For the nondefault Windows toolchain install, `RTOOLS45_HOME` pointed to
the Rtools45 directory. The PATH began with its `usr/bin` and
`x86_64-w64-mingw32.static.posix/bin` directories. R 4.5.1 ran with an
isolated user library. The retry used the 2026-08-07 Posit CRAN snapshot and
the locked Bioconductor 3.21 repositories. The build retained the package
versions recorded in `renv.lock`. R Arrow was excluded from the retry because
the pipeline supports gzip CSV with Parquet conversion in Python.
[Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html)

The August 2026 repository lacked the locked Windows `stringi` binary and
fell back to compiling bundled ICU. An earlier Posit snapshot supplied the
exact `stringi` 1.8.7 Windows binary. The ZIP passed CRC and package MD5
checks, then loaded on R 4.5.1 and reported ICU 74.1. Its SHA-256 was
`13d3463d867eed30040aab6473cca1420d1a0fc84bac3b7309ee1a031386a5c0`.
The native retry installed that verified binary and reused the packages
already restored. [Posit snapshot package](https://packagemanager.posit.co/cran/2025-09-08/bin/windows/contrib/4.5/stringi_1.8.7.zip)

Native verification found 219 exact package-version matches. Optional
`arrow` 25.0.0 was the only absent lock entry. The installed
`MCPcounter` DESCRIPTION reported the exact Git SHA and subdirectory from
the lock. Some Windows binaries were built under a different R 4.5 patch
release. The analysis runtime was R 4.5.1.

## Fresh raw-CEL limma results

The fresh raw-CEL limma run completed in 52.9 s after restoring the pinned
`R.utils` package for the gzip CSV reader. All eight model/contrast rows
matched the archived tested universe of 18,490 genes. The largest absolute
log2 fold-change difference was 2.1 × 10⁻¹⁴ and the largest BH FDR difference
was 3.9 × 10⁻¹³. Every significant-gene count and corresponding direction
agreed. The adjusted field contrast retained 1,580 hits and the paired
tumor-versus-distant-mucosa contrast retained 4,502. This reproduces the
archived statistical outputs to numerical precision on the current archive.
Historical TAR and CEL byte identity remains unresolved.

Covariate sensitivity reproduced 1,662, 1,580 and 1,521 significant genes
across U0, U1 and U2 and the 101/1,370/7 evidence-tier counts. The fresh
composition model retained effect size, direction and FDR for 74 of 101
high-confidence genes. Ten of ten flagged stromal/ECM genes attenuated by
at least 25%. The stress-gene sensitivity flagged three high-confidence
genes, leaving 98. These are within-cohort sensitivities.

The full raw-versus-deposited pipeline compared 15,879 shared genes using
corresponding unadjusted and paired designs. The log-fold-change Pearson
correlations were 0.9977–0.9979 and sign agreement was 99.4–99.6%. This
comparison uses the pipeline's separate deposited mapping. The first audit's
same-annotation limma refits held more inputs fixed.

## Fresh grouped nested CV

The original three-task, three-model run completed in 1,617.4 s. It emitted
135 fold metrics and 5,310 held-out predictions. All 135 fold-level macro-F1
values and selected hyperparameter sets matched the archive. Mean outer-fold
elastic-net macro-F1 was 0.9413, 0.9926 and 0.9847 for Tasks A, B and C.
Patient/donor bootstrap intervals were generated from fresh aggregated
out-of-fold predictions using 1,000 draws per task and model.
Fresh compact-panel analysis reproduced the archived summary and selected
CLC/DYNC1H1/FOS/VIP/SNORA12 for Task B and FOXQ1/CEMIP/ETV4 for Task C.
The original Task B confounding stage also ran. Its archived numeric rows
reproduced to within 6.7 × 10⁻¹⁴, including macro-F1 0.9772 after the
prespecified age, sex and location residualization of the selected panel.

All three original group-preserving permutation tests completed in
2,324.6 s, using 1,000 permutations for each prespecified null. The observed
macro-F1, null means and Monte Carlo p-values reproduced the archive.
Each p-value was 1/1,001. The Task B null concerns the observed group
association, which includes the acquisition protocol.

## Task B raw versus deposited expression

The continuation reran the same grouped, fold-local Task B ablations on the
fresh RMA gene matrix from the current GSE44076 archive. It reused the same
metadata, annotation coverage, seeds, panel rule and scoring code as the
deposited-expression analysis. All six macro-F1 values agreed exactly. The
adjusted sensitivity retained macro-F1 0.8431 (conditional donor-bootstrap
interval 0.7769–0.9008). ROC-AUC was 0.9318 on fresh raw RMA and 0.9320 on
deposited expression in that sensitivity. Both analyses use the same cohort
and inherit the acquisition alias.

The prespecified four-gene GSE4107 transfer sensitivity used the current
GSE44076 raw-CEL gene matrix for training and the same deposited GSE4107
evaluation matrix. At the frozen 0.50 threshold, it classified all 22 donors
as healthy. Balanced accuracy was 0.500, macro-F1 0.313 and ROC-AUC 0.983
(deposited training ROC-AUC 0.992). The direction and operating-point failure
persisted. The score remains sensitive to the shared acquisition alias.

## Task C raw-CEL reproduction

The original nested and compact-panel stages recovered the archived
FOXQ1/CEMIP/ETV4 panel and hyperparameters. Its fresh threshold and external
validation stages completed and emitted the same predictions as the
separate frozen-panel refit described below. The saved external result
contract and all six targeted model/prediction tests passed.

A separate run held the archived FOXQ1/CEMIP/ETV4 panel and hyperparameters
fixed, refit on 196 fresh GSE44076 RMA arrays, and selected the transport
threshold on repeated GSE44076 grouped OOF predictions. The label-independent
platform intersection retained 12,039 genes, including CEMIP and ETV4 from
the panel. The threshold relocked at 0.41 before GSE41258 outcomes were read.
This refit assesses reproducibility conditional on the archived selection.

On the canonical GSE41258 set (233 arrays, 190 patients), fresh rank-transport
ROC-AUC was 0.9829, macro-F1 0.8282, sensitivity 1.0000, specificity 0.5577,
Brier 0.0533 and calibration intercept/slope −2.0661/2.6198. All 13 primary
patient-cluster bootstrap interval endpoints reproduced exactly. The
prespecified 0.50 sensitivity had macro-F1 0.8898 and specificity 0.7115.
These retrospective operating points quantify threshold and calibration
transport on this cohort. Neither external outcome set entered fitting,
signature selection, or threshold selection.

## Task C sensitivity to deposited expression

The deposited-expression refit held the panel, hyperparameters, mapping,
12,039-gene intersection, development labels and 233 external evaluation
arrays fixed. The GEO-deposited matrices supplied expression. GSE44076
alone selected its threshold of 0.42. The contrast also persists at the
same 0.41 cutoff.

| Expression route and operating point | ROC-AUC | Macro-F1 | Sensitivity | Specificity | Brier | Calibration intercept / slope |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Fresh raw RMA, locked 0.41 | 0.983 | 0.828 | 1.000 | 0.558 | 0.053 | −2.066 / 2.620 |
| GEO deposited, locked 0.42 | 0.986 | 0.951 | 0.972 | 0.942 | 0.032 | 1.372 / 1.420 |
| GEO deposited, raw reference cutoff 0.41 | 0.986 | 0.951 | 0.972 | 0.942 | 0.032 | 1.372 / 1.420 |

Discrimination was similar across expression routes. The specificity and
calibration differences show sensitivity to the expression inputs even
after within-sample ranking. Raw RMA remains the primary analysis. The
deposited result uses the same retrospective external cohort.

## Reproduction boundaries

A run using current GEO pins is a fresh raw-CEL run on a separately identified
input set. Similar resulting estimates would support scientific robustness
to that input version. Exact reproduction requires the historical bytes.
GSE41258 and GSE4107 remain evaluation cohorts. The JGAD000785 study lists
endoscopic biopsies from 55 colorectal cancer cases and 17 healthy
participants in its RNA-seq cohort, with multiple mucosa sites per donor.
An acquisition-matched Task B test requires an approved Type I data-use
request and eligible analysis environment, followed by a patient/site
eligibility audit. GSE44076 sampled mucosa 10 cm from the tumor. The reported
JGAD000785 1–2 cm specimens and other anatomical sites require a prespecified
distance and site rule before any evaluation. RNA-seq transfer introduces a
platform change. [NBDC dataset](https://humandbs-production.ddbj.nig.ac.jp/dataset/JGAD000785),
[NBDC study](https://humandbs-production.ddbj.nig.ac.jp/research/hum0430)
