from __future__ import annotations

import joblib
import numpy as np
import pandas as pd
import pytest
import yaml
from sklearn.metrics import roc_auc_score

from src.data.cross_platform import common_gene_symbols
from src.data.provenance import expression_path, result_root
from src.modeling.stability import feature_stability

pytestmark = pytest.mark.full_data


def test_external_gene_intersection(root):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    provenance = config["modeling"]["expression_provenance"]
    primary = pd.read_parquet(expression_path(root / "data/processed", "GSE44076", "gene", provenance))
    external = pd.read_parquet(expression_path(root / "data/processed", "GSE41258", "gene", provenance))
    common = common_gene_symbols(primary, external)
    saved = (root / f"data/processed/common_genes_GSE44076_GSE41258_{provenance}.txt").read_text(encoding="utf-8").splitlines()
    assert common == saved
    assert len(common) > 10000


def test_model_serialization_and_reload(root):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    provenance = config["modeling"]["expression_provenance"]
    path = result_root(root / "models", provenance) / "task_c_primary_full_signature_model.joblib"
    if not path.exists():
        import pytest

        pytest.skip("Locked model is created by scripts/run_modeling.py")
    artifact = joblib.load(path)
    assert artifact["feature_genes"]
    assert set(artifact["model"].named_steps["model"].classes_) == {"adjacent_normal", "tumor"}


def test_external_metrics_match_saved_predictions(root):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    metrics_root = result_root(root / "results/metrics", config["modeling"]["expression_provenance"])
    predictions = pd.read_csv(metrics_root / "external_validation_predictions.csv")
    metrics = pd.read_csv(metrics_root / "external_validation_metrics.csv").set_index(
        ["representation", "evaluation_set", "threshold_policy"]
    )
    for keys, frame in predictions.groupby(
        ["representation", "evaluation_set", "threshold_policy"]
    ):
        # External labels are deliberately preserved instead of being silently
        # rewritten to the primary cohort's shorter class names.
        assert set(frame["y_true"]) <= {"normal_colon", "primary_tumor"}
        y_true = frame["y_true"].eq("primary_tumor").astype(int)
        observed = roc_auc_score(y_true, frame["probability_primary_tumor"])
        assert np.isclose(observed, metrics.loc[keys, "roc_auc"])


def test_locked_signature_is_primary_derived_and_cross_platform_subset(root):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    models_root = result_root(root / "models", config["modeling"]["expression_provenance"])
    artifact = joblib.load(models_root / "task_c_primary_full_signature_model.joblib")
    primary_samples = set(
        pd.read_csv(root / "data/metadata/gse44076_samples.csv")["geo_accession"]
    )
    external_samples = set(
        pd.read_csv(root / "data/metadata/gse41258_samples.csv")["geo_accession"]
    )
    assert set(artifact["training_sample_ids"]).issubset(primary_samples)
    assert set(artifact["training_sample_ids"]).isdisjoint(external_samples)
    rank_artifact = joblib.load(models_root / "task_c_cross_platform_transport_model.joblib")
    assert set(rank_artifact["feature_genes"]).issubset(artifact["feature_genes"])
    assert rank_artifact["refitting_occurred"] is True
    assert artifact["refitting_occurred"] is False


def test_feature_stability_counts_repeat_fold_keys_once():
    coefficients = pd.DataFrame(
        {
            "task": ["task_c"] * 4,
            "repeat": [0, 0, 1, 1],
            "outer_fold": [0, 0, 0, 1],
            "class": ["tumor"] * 4,
            "gene_symbol": ["GENE1"] * 4,
            "coefficient": [1.0, 1.1, 0.9, 1.2],
            "selected_nonzero": [True] * 4,
        }
    )
    summary = feature_stability(coefficients, total_outer_folds=4)
    assert int(summary.loc[0, "selected_folds"]) == 3
    assert summary.loc[0, "selection_frequency"] == 0.75
