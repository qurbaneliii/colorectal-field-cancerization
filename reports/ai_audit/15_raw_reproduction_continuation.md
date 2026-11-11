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

## Reproduction boundaries

A run using current GEO pins is a fresh raw-CEL run on a separately identified
input set. Similar resulting estimates would support scientific robustness
to that input version. Exact reproduction requires the historical bytes.
GSE41258 and GSE4107 remain evaluation cohorts. JGAD000785 remains controlled
access and requires an approved data-use request and eligible analysis
environment. [NBDC dataset](https://humandbs-production.ddbj.nig.ac.jp/dataset/JGAD000785)
