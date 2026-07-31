from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import yaml

from src.data.provenance import expression_path
from src.data.validation import validate_expression

pytestmark = pytest.mark.full_data


def test_processed_expression_orientation_and_integrity(root):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    provenance = config["modeling"]["expression_provenance"]
    for accession in ("GSE44076", "GSE41258"):
        expression = pd.read_parquet(expression_path(root / "data/processed", accession, "gene", provenance))
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
        if provenance == "raw_cel_rma":
            assert np.nanmedian(values) < 20


def test_unique_sample_identifiers(root):
    for accession in ("gse44076", "gse41258"):
        metadata = pd.read_csv(root / f"data/metadata/{accession}_samples.csv")
        assert metadata["geo_accession"].is_unique
        assert metadata.loc[metadata["inclusion_status"].eq("included"), "tissue_class"].notna().all()
