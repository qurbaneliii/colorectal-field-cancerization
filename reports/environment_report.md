# Environment report

## Python

- Interpreter: CPython 3.12.10
- Isolated environment: `.venv-publication`
- Lock: `requirements-lock.txt`
- Editable package installed with `--no-deps`
- `pip check`: no broken requirements
- Runtime/test imports: passed

The pre-existing `.venv` is not the publication environment and contains an
unrelated inherited OpenTelemetry/protobuf conflict. No result in this release
uses that environment.

## R

- R: 4.5.1
- Bioconductor: 3.21
- Lock: populated `renv.lock`
- `R/verify_environment.R`: passed

Some binary/source packages report that they were built under a later R 4.5
patch release. They loaded and executed successfully under R 4.5.1; this is a
warning, not a package-load failure.

## Host tools

- Git: available
- Docker client/Desktop: installed
- Docker Linux engine: unavailable because the host WSL service is disabled
- GNU Make: unavailable on this Windows host; the equivalent commands were
  executed directly
