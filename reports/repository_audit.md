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
- Rscript: unavailable: [WinError 2] The system cannot find the file specified
- Git: 62cce06e252cfc01582ab8bc1664d8c1cc037834

## Missing or blocked requirements

- R is not installed in the execution environment, so raw-CEL RMA,
  Bioconductor QC, limma differential expression, and enrichment cannot be
  executed here. Reproducible R scripts and a Docker environment are provided.
- The deposited normalized GEO series matrices remain usable for the explicitly
  labeled processed-matrix Python route.

## Reusable components and modifications

No prior components existed. The repository now contains deterministic metadata
parsers, provenance capture, independent platform annotation, group-aware
modeling, external validation, tests, and manuscript/report generators.
