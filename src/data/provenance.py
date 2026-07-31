from __future__ import annotations

from pathlib import Path


VALID_PROVENANCE = {"raw_cel_rma", "geo_deposited_series_matrix"}


def expression_path(processed_root: Path, accession: str, level: str, provenance: str) -> Path:
    if provenance not in VALID_PROVENANCE:
        raise ValueError(f"Unsupported expression provenance: {provenance}")
    if level not in {"probe", "gene"}:
        raise ValueError(f"Unsupported expression level: {level}")
    return processed_root / f"{accession}_{level}_expression_{provenance}.parquet"


def sample_metadata_path(processed_root: Path, accession: str, provenance: str) -> Path:
    if provenance not in VALID_PROVENANCE:
        raise ValueError(f"Unsupported expression provenance: {provenance}")
    return processed_root / f"{accession}_sample_metadata_{provenance}.csv"


def result_root(base: Path, provenance: str) -> Path:
    """Raw-CEL results use publication paths; sensitivity outputs are namespaced."""
    if provenance not in VALID_PROVENANCE:
        raise ValueError(f"Unsupported expression provenance: {provenance}")
    return base if provenance == "raw_cel_rma" else base / provenance
