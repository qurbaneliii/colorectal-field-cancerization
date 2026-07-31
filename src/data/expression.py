from __future__ import annotations

import csv
import gzip
import io
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


def read_geo_series_matrix(path: Path) -> pd.DataFrame:
    """Read only the expression table from a GEO series matrix."""
    begin_line = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("!series_matrix_table_begin"):
                begin_line = line_number + 1
                break
    if begin_line is None:
        raise ValueError(f"No expression table found in {path}")
    frame = pd.read_csv(
        path,
        sep="\t",
        skiprows=begin_line,
        comment="!",
        quotechar='"',
        low_memory=False,
    )
    if frame.columns[0] != "ID_REF":
        raise ValueError(f"Unexpected first series-matrix column: {frame.columns[0]}")
    frame = frame.rename(columns={"ID_REF": "probe_id"}).set_index("probe_id")
    frame = frame.apply(pd.to_numeric, errors="raise")
    if frame.index.duplicated().any():
        raise ValueError("Duplicate probe IDs in series matrix")
    values = frame.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Series matrix contains missing or infinite expression values")
    return frame


def _open_annotation_stream(url: str):
    response = urllib.request.urlopen(url, timeout=180)
    if url.lower().endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(fileobj=response), encoding="utf-8", errors="replace")
    return io.TextIOWrapper(response, encoding="utf-8", errors="replace")


def fetch_platform_mapping(
    url: str,
    id_column: str,
    symbol_column: str,
    entrez_column: str,
    output_path: Path,
) -> pd.DataFrame:
    """Fetch the official GEO annotation table and retain an unambiguous probe mapping."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with _open_annotation_stream(url) as handle:
        header: list[str] | None = None
        rows: list[dict[str, str]] = []
        for line in handle:
            stripped = line.rstrip("\r\n")
            if not stripped or stripped.startswith("#") or stripped.startswith("!platform_table_begin"):
                continue
            if stripped.startswith("!platform_table_end"):
                break
            fields = next(csv.reader([stripped], delimiter="\t"))
            if header is None:
                if id_column in fields and symbol_column in fields:
                    header = fields
                continue
            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))
            row = dict(zip(header, fields, strict=False))
            rows.append(
                {
                    "probe_id": row.get(id_column, "").strip(),
                    "gene_symbol_raw": row.get(symbol_column, "").strip(),
                    "entrez_id_raw": row.get(entrez_column, "").strip(),
                }
            )
    if header is None or not rows:
        raise ValueError(f"Could not parse GEO platform annotation from {url}")
    raw = pd.DataFrame(rows)

    def one_value(value: str) -> str | None:
        value = value.strip()
        if not value or value == "---":
            return None
        for separator in (" /// ", " // ", ";"):
            parts = [part.strip() for part in value.split(separator) if part.strip() and part != "---"]
            if len(parts) > 1:
                return parts[0] if len(set(parts)) == 1 else None
            if parts:
                value = parts[0]
        return value if value and value != "---" else None

    raw["gene_symbol"] = raw["gene_symbol_raw"].map(one_value)
    raw["entrez_id"] = raw["entrez_id_raw"].map(one_value)
    raw = raw.dropna(subset=["probe_id", "gene_symbol"])
    probe_symbol_counts = raw.groupby("probe_id")["gene_symbol"].nunique()
    valid_probes = probe_symbol_counts[probe_symbol_counts.eq(1)].index
    mapping = (
        raw[raw["probe_id"].isin(valid_probes)]
        .drop_duplicates(["probe_id", "gene_symbol"])
        .sort_values(["gene_symbol", "probe_id"])
        .reset_index(drop=True)
    )
    mapping.to_csv(output_path, index=False)
    return mapping


def aggregate_probe_expression(
    probe_expression: pd.DataFrame, mapping: pd.DataFrame
) -> pd.DataFrame:
    shared = probe_expression.index.intersection(mapping["probe_id"])
    if shared.empty:
        raise ValueError("No annotated probes overlap the expression matrix")
    ordered_mapping = mapping.drop_duplicates("probe_id").set_index("probe_id").loc[shared]
    annotated = probe_expression.loc[shared].copy()
    annotated.insert(0, "gene_symbol", ordered_mapping["gene_symbol"])
    gene = annotated.groupby("gene_symbol", sort=True).median(numeric_only=True)
    if gene.index.duplicated().any():
        raise AssertionError("Gene aggregation did not create unique symbols")
    return gene


def expression_for_model(
    gene_expression: pd.DataFrame, metadata: pd.DataFrame, included_only: bool = True
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    selected = metadata.copy()
    if included_only:
        selected = selected[selected["inclusion_status"].eq("included")]
    samples = selected["geo_accession"].tolist()
    missing = set(samples).difference(gene_expression.columns)
    if missing:
        raise ValueError(f"Expression missing samples: {sorted(missing)[:5]}")
    x = gene_expression.loc[:, samples].T.to_numpy(dtype=np.float64)
    y = selected["tissue_class"].to_numpy()
    groups = selected["donor_or_patient_group"].astype(str).to_numpy()
    sample_ids = selected["geo_accession"].to_numpy()
    return x, y, groups, sample_ids, gene_expression.index.astype(str).tolist()
