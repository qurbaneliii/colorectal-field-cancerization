from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts import run_development_audit as development
from scripts.run_data_audit import write_covariate_audit
from scripts.select_task_c_threshold import threshold_grid
from scripts.verify_external_results import read_components, verify
from src.data import metadata
from src.data.expression import read_geo_series_matrix
from src.data.provenance import artifact_sha256
from src.modeling.external_result import external_summary
from src.modeling.external_validation import locked_binary_predictions


def test_active_external_snapshot_replays_saved_predictions(root):
    assert verify(root)["locked_threshold"] == pytest.approx(0.41)


def test_manuscript_and_primary_figure_use_current_external_evidence(root):
    result = verify(root)["primary_evaluation"]
    manuscript = " ".join((root / "manuscript/results.md").read_text(encoding="utf-8").split())
    for label, key in (
        ("ROC-AUC was", "roc_auc"), ("macro-F1", "f1_macro"),
        ("balanced accuracy", "balanced_accuracy"), ("specificity", "specificity"),
    ):
        assert f"{label} {result[key]:.3f}" in manuscript
    figure = (root / "results/figures/external_validation_primary.svg").read_text(encoding="utf-8")
    assert f"n={result['arrays']} arrays; {result['unique_patients']} patients" in figure
    normal, tumor = result["normal_colon_arrays"], result["primary_tumor_arrays"]
    tn = round(normal * result["specificity"])
    tp = round(tumor * result["sensitivity"])
    assert all(f"<!-- {count} -->" in figure for count in (tn, normal - tn, tumor - tp, tp))


def test_stale_external_metrics_and_threshold_fail(root):
    sources = read_components(root)
    sources["metrics"].loc[sources["metrics"]["is_primary"], "specificity"] = 0.98
    with pytest.raises(ValueError, match="specificity disagrees"):
        external_summary(**sources)
    sources = read_components(root)
    sources["full_card"]["transport_threshold"] = 0.50
    with pytest.raises(ValueError, match="threshold provenance"):
        external_summary(**sources)


def test_external_probability_columns_follow_fitted_class_order():
    model = SimpleNamespace(
        named_steps={"model": SimpleNamespace(classes_=["tumor", "adjacent_normal"])},
        predict_proba=lambda _: np.array([[0.8, 0.2], [0.1, 0.9]]),
    )
    frame, metrics = locked_binary_predictions(
        model,
        np.empty((2, 0)),
        np.array(["A", "B"]),
        np.array(["1", "2"]),
        np.array(["primary_tumor", "normal_colon"]),
        "test",
    )
    assert frame["probability_primary_tumor"].tolist() == [0.8, 0.1]
    assert metrics["roc_auc"] == metrics["specificity"] == 1
    assert metrics["log_loss"] == pytest.approx(-np.log([0.8, 0.9]).mean())


def test_threshold_at_exact_grid_value_uses_inclusive_cutoff():
    config = {
        "threshold_grid_start": 0.10,
        "threshold_grid_stop": 0.60,
        "threshold_grid_step": 0.10,
        "default_threshold": 0.50,
    }
    grid = threshold_grid(np.array([0.49, 0.50]), np.array(["adjacent_normal", "tumor"]), config)
    assert grid.loc[grid["selected"], "threshold"].item() == 0.50
    assert grid.loc[grid["selected"], "balanced_accuracy"].item() == 1.0


def test_nested_threshold_predictions_exclude_every_scored_patient(monkeypatch):
    groups = np.repeat(np.arange(16), 2)
    labels = np.tile(["adjacent_normal", "tumor"], 16)
    values = np.column_stack([groups, labels == "tumor"])
    cohorts = pd.DataFrame(
        {"tissue_class": labels, "donor_or_patient_group": groups, "geo_accession": range(32)}
    )
    captured_fit_sizes = []

    def fit(x, y, fit_groups, genes, seed, config):
        captured_fit_sizes.append(len(set(fit_groups)))
        return set(fit_groups), [0, 1], genes, {}

    def score(fit_groups, x, label):
        assert not fit_groups.intersection(x[:, 0])
        return np.where(x[:, 1] == 1, 0.8, 0.2)

    monkeypatch.setattr(development, "fit_rank_policy", fit)
    monkeypatch.setattr(development, "positive_scores", score)
    predictions, decisions = development.task_c_nested_thresholds(
        values,
        cohorts,
        ["group", "simulated_signal"],
        {
            "modeling": {"seeds": [44076], "outer_splits": 4, "inner_splits": 3},
            "external_validation": {
                "threshold_grid_start": 0.1,
                "threshold_grid_stop": 0.9,
                "threshold_grid_step": 0.1,
                "default_threshold": 0.5,
            },
        },
    )
    assert predictions["sample_id"].nunique() == 32
    assert len(predictions) == 32
    assert len(decisions) == 4
    assert set(captured_fit_sizes) == {8, 12}


def test_external_metadata_preserves_scan_date_source(monkeypatch, tmp_path):
    records = {
        "!Sample_geo_accession": [["GSM1"]],
        "!Sample_title": [["Primary Tumor test"]],
        "!Sample_platform_id": [["GPL96"]],
        "!Sample_supplementary_file": [["ftp://example.com/GSM1.CEL.gz"]],
        "!Sample_characteristics_ch1": [
            ["tissue: Primary Tumor"], ["patient id: 1"], ["included in analysis: Yes"],
            ["scan date: 37659"],
        ],
    }
    monkeypatch.setattr(metadata, "parse_series_metadata", lambda _: records)
    monkeypatch.setattr(metadata, "tar_members", lambda _: ["GSM1.CEL.gz"])
    frame = metadata.build_gse41258_metadata(tmp_path, tmp_path, r"_ez|rehyb")
    assert frame.loc[0, "scan_date"] == "37659"
    assert frame.loc[0, "inclusion_status"] == "included"


def test_versioned_text_hash_is_checkout_independent(tmp_path):
    path = tmp_path / "result.csv"
    path.write_bytes(b"x,y\r\n1,2\r\n")
    expected = artifact_sha256(path)
    path.write_bytes(b"x,y\n1,2\n")
    assert artifact_sha256(path) == expected


def test_missing_expression_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_geo_series_matrix(tmp_path / "unavailable_series_matrix.txt.gz")


def test_missing_stage_and_raw_scan_dates_are_excluded_from_screening(tmp_path):
    cohort = pd.DataFrame({
        "inclusion_status": ["included"] * 24,
        "tissue_class": ["normal_colon"] * 12 + ["primary_tumor"] * 12,
        "stage": ["--"] * 12 + ["II", "III"] * 6,
        "scan_date": ["25/02/2004", "37659"] * 12,
    })
    write_covariate_audit({"metadata": tmp_path, "reports": tmp_path}, cohort, cohort)
    audit = pd.read_csv(tmp_path / "covariate_missingness_usability.csv")
    assert audit.loc[audit["covariate"].eq("stage"), "available_values"].tolist() == [12, 12]
    dates = audit.loc[audit["covariate"].eq("scan_date")]
    assert dates["available_values"].tolist() == [24, 24]
    assert not dates["usable_for_association"].any()
