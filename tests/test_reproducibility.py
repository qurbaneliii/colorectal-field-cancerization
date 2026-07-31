from __future__ import annotations

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from src.data.cross_platform import common_gene_symbols


def test_external_gene_intersection(root):
    primary = pd.read_parquet(root / "data/processed/GSE44076_gene_expression.parquet")
    external = pd.read_parquet(root / "data/processed/GSE41258_gene_expression.parquet")
    common = common_gene_symbols(primary, external)
    saved = (
        root / "data/processed/common_genes_GSE44076_GSE41258.txt"
    ).read_text(encoding="utf-8").splitlines()
    assert common == saved
    assert len(common) > 10000


def test_model_serialization_and_reload(root):
    path = root / "models/task_c_locked_elastic_net.joblib"
    if not path.exists():
        import pytest

        pytest.skip("Locked model is created by scripts/run_modeling.py")
    artifact = joblib.load(path)
    assert artifact["signature_genes"]
    assert set(artifact["classes"]) == {"adjacent_normal", "tumor"}


def test_external_metrics_match_saved_predictions(root):
    predictions = pd.read_csv(root / "results/metrics/external_validation_predictions.csv")
    metrics = pd.read_csv(root / "results/metrics/external_validation_metrics.csv").set_index(
        "representation"
    )
    for representation, frame in predictions.groupby("representation"):
        y_true = frame["y_true"].eq("tumor").astype(int)
        observed = roc_auc_score(y_true, frame["probability_tumor"])
        assert np.isclose(observed, metrics.loc[representation, "roc_auc"])


def test_locked_signature_is_primary_derived_and_cross_platform_subset(root):
    artifact = joblib.load(root / "models/task_c_locked_elastic_net.joblib")
    primary_samples = set(
        pd.read_csv(root / "data/metadata/gse44076_samples.csv")["geo_accession"]
    )
    external_samples = set(
        pd.read_csv(root / "data/metadata/gse41258_samples.csv")["geo_accession"]
    )
    assert set(artifact["training_samples"]).issubset(primary_samples)
    assert set(artifact["training_samples"]).isdisjoint(external_samples)
    rank_artifact = joblib.load(root / "models/task_c_locked_rank_elastic_net.joblib")
    assert set(rank_artifact["signature_genes"]).issubset(artifact["signature_genes"])
