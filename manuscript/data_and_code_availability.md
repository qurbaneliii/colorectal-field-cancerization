# Data and code availability

Expression data are available from NCBI GEO under GSE44076 and GSE41258. Raw
archives and large generated matrices are intentionally excluded from Git.
Checksums, deterministic metadata, compact result tables, environment locks,
commands, and model documentation are included in the repository. Users must
place the five named GEO files under `data/raw/<accession>/` and run `make all`.
