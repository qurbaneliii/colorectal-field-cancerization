# Raw-CEL reproduction continuation

This audit continues the analysis at `0ede337ffc360bc96dea6724f974bd5ef7f184d4`
on PR #4. The historical manifest remains an immutable record of the earlier
analysis. The fresh run uses a separate worktree and records actual inputs.

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

The initial Linux restore built 176 packages and failed at `fs` (missing
`libuv1-dev`), `oligo` (the Bioconductor source endpoint did not serve the
locked 1.71.7), and `MCPcounter` (a duplicate staging path). The hosted job
now installs the system library and builds the two source tarballs from the
revisions recorded in the lock. The Windows host also installed Rtools45
and compiled and linked an R test DLL with GCC 14.3.0. The Windows scientific
run used the locked core packages and the supported gzip CSV/Parquet route.
[oligo 1.71.7 source](https://github.com/bioc/oligo/commit/e073acb9f29886d4cc073f634c109d0bb76b2f81),
[MCPcounter source](https://github.com/ebecht/MCPcounter/tree/b6eac73e91c246fcff0bb1a5c68a816cd588fc48)

## Fresh raw-CEL limma results

The fresh raw-CEL limma run completed in 52.9 s after restoring the pinned
`R.utils` package for the gzip CSV reader. All eight model/contrast rows
matched the archived tested universe of 18,490 genes. The largest absolute
log2 fold-change difference was 2.1 × 10⁻¹⁴ and the largest BH FDR difference
was 3.9 × 10⁻¹³. Every significant-gene count and corresponding direction
agreed. The adjusted field contrast retained 1,580 hits and the paired
tumor-versus-distant-mucosa contrast retained 4,502. This reproduces the
archived statistical outputs to numerical precision on the current archive.
It does not resolve historical TAR or CEL byte identity.

Covariate sensitivity reproduced 1,662, 1,580 and 1,521 significant genes
across U0, U1 and U2 and the 101/1,370/7 evidence-tier counts. The fresh
composition model retained effect size, direction and FDR for 74 of 101
high-confidence genes. Ten of ten flagged stromal/ECM genes attenuated by
at least 25%. The stress-gene sensitivity flagged three high-confidence
genes, leaving 98. These are within-cohort sensitivities.

## Task B raw versus deposited expression

The continuation reran the same grouped, fold-local Task B ablations on the
fresh RMA gene matrix from the current GSE44076 archive. It reused the same
metadata, annotation coverage, seeds, panel rule and scoring code as the
deposited-expression analysis. All six macro-F1 values agreed exactly. The
adjusted sensitivity retained macro-F1 0.8431 (conditional donor-bootstrap
interval 0.7769–0.9008). ROC-AUC was 0.9318 on fresh raw RMA and 0.9320 on
deposited expression in that sensitivity. The agreement concerns the same
cohort and does not resolve the acquisition alias.

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
