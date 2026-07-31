from __future__ import annotations

import numpy as np
import pandas as pd

from src.data.validation import validate_expression


def test_processed_expression_orientation_and_integrity(root):
    for accession in ("GSE44076", "GSE41258"):
        expression = pd.read_parquet(
            root / f"data/processed/{accession}_gene_expression.parquet"
        )
        metadata = pd.read_csv(
            root / f"data/metadata/{accession.lower()}_samples.csv",
            dtype={"patient_id": str},
        )
        validate_expression(expression, metadata)
        assert expression["gene_symbol"].is_unique
        assert len(expression) > len(expression.columns) - 1
        values = expression.drop(columns="gene_symbol").to_numpy()
        assert np.isfinite(values).all()
        assert not np.isnan(values).any()


def test_unique_sample_identifiers(root):
    for accession in ("gse44076", "gse41258"):
        metadata = pd.read_csv(root / f"data/metadata/{accession}_samples.csv")
        assert metadata["geo_accession"].is_unique
        assert metadata.loc[metadata["inclusion_status"].eq("included"), "tissue_class"].notna().all()
