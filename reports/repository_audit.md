# Repository audit

## Audited state

The audit preserves the existing hybrid R/Python implementation and reads the
five deposited GEO inputs without modifying `data/raw/`.

## Data discovered

| accession   | filename                      |   size_bytes |
|:------------|:------------------------------|-------------:|
| GSE41258    | GSE41258_RAW.tar              |   1378918400 |
| GSE41258    | GSE41258_clinical_data.txt.gz |         3135 |
| GSE41258    | GSE41258_series_matrix.txt.gz |     17537584 |
| GSE44076    | GSE44076_RAW.tar              |    532213760 |
| GSE44076    | GSE44076_series_matrix.txt.gz |     32596137 |

- GSE44076 raw TAR members: 246
- GSE41258 raw TAR members: 390
- Raw files were read but not modified.

## Runtime

- Python: 3.12.10
- Rscript: Rscript (R) version 4.5.1 (2025-06-13)
- Git: 9b43cc12c879c751f40576d33ae5719e42702eb4

## Missing or blocked requirements

- R was discovered and the raw-CEL, Bioconductor QC, limma, composition, and enrichment stages have executable repository routes.
- Docker verification is a separate host-level gate and is not inferred from
  the presence of a Dockerfile.

## Reusable components and modifications

The existing hybrid R/Python pipeline is reused and extended with deterministic
metadata parsing, provenance capture, independent platform annotation,
group-aware modeling, external validation, tests, and manuscript/report
generation.
