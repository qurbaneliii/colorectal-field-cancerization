from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest
import yaml

from scripts.run_compact_panel_analysis import choose_panel
from scripts.select_task_c_threshold import threshold_grid
from src.modeling.residualization import TrainingCompositionPCs, TrainingCovariateResidualizer
from src.modeling.stability import feature_stability


def test_training_only_residualizer_does_not_refit_on_test_rows():
    train_covariates = np.array([[0.0], [1.0], [2.0], [3.0]])
    train_expression = np.column_stack(
        [2.0 + 4.0 * train_covariates[:, 0], 7.0 - train_covariates[:, 0]]
    )
    residualizer = TrainingCovariateResidualizer().fit(
        train_expression, train_covariates
    )
    coefficients = residualizer.coefficients_.copy()
    mean = residualizer.covariate_mean_.copy()
    test_covariates = np.array([[1000.0], [-1000.0]])
    test_expression = np.zeros((2, 2))
    residualizer.transform(test_expression, test_covariates)
    np.testing.assert_array_equal(residualizer.coefficients_, coefficients)
    np.testing.assert_array_equal(residualizer.covariate_mean_, mean)


def test_composition_pcs_are_fit_on_training_rows_only():
    train = np.array(
        [[0.0, 1.0, 2.0], [1.0, 2.0, 4.0], [2.0, 4.0, 7.0], [3.0, 8.0, 11.0]]
    )
    transformer = TrainingCompositionPCs(n_components=2).fit(train)
    training_mean = transformer.scaler_.mean_.copy()
    components = transformer.pca_.components_.copy()
    transformed = transformer.transform(np.array([[1000.0, -500.0, 250.0]]))
    assert transformed.shape == (1, 2)
    np.testing.assert_array_equal(transformer.scaler_.mean_, training_mean)
    np.testing.assert_array_equal(transformer.pca_.components_, components)


def test_panel_rule_uses_only_inner_validation_summary():
    inner = pd.DataFrame(
        {
            "inner_fold": [0, 1, 0, 1],
            "panel_size": ["3", "3", "5", "5"],
            "effective_gene_count": [3, 3, 5, 5],
            "f1_macro": [0.98, 0.98, 0.99, 0.99],
            "balanced_accuracy": [0.98, 0.98, 0.99, 0.99],
            "brier_score": [0.02, 0.02, 0.01, 0.01],
            "log_loss": [0.10, 0.10, 0.08, 0.08],
        }
    )
    modeling = {
        "compact_panel_practical_tolerance": 0.02,
        "compact_panel_min_balanced_accuracy": 0.80,
        "compact_panel_max_brier_degradation": 0.02,
    }
    selected, decision = choose_panel(inner, modeling)
    assert selected == "3"
    assert decision.loc[decision["selected"], "panel_size"].tolist() == ["3"]


def test_threshold_grid_uses_configured_rule_and_grid(root: Path):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    external = config["external_validation"]
    probability = np.array([0.1, 0.3, 0.7, 0.9])
    labels = np.array(["adjacent_normal", "adjacent_normal", "tumor", "tumor"])
    curve = threshold_grid(probability, labels, external)
    expected = int(
        round(
            (external["threshold_grid_stop"] - external["threshold_grid_start"])
            / external["threshold_grid_step"]
        )
        + 1
    )
    assert len(curve) == expected
    assert curve["selection_rule"].nunique() == 1


def test_strict_and_consensus_stability_labels_are_disjoint():
    coefficients = pd.DataFrame(
        {
            "task": ["task_b"] * 5,
            "repeat": [0, 0, 0, 0, 0],
            "outer_fold": [0, 1, 2, 3, 4],
            "class": ["adjacent_normal"] * 5,
            "gene_symbol": ["GENE"] * 5,
            "coefficient": [1.0, 1.0, 1.0, -1.0, -1.0],
            "selected_nonzero": [True, True, True, True, True],
        }
    )
    summary = feature_stability(
        coefficients,
        total_outer_folds=5,
        strict_frequency=0.65,
        strict_sign_consistency=0.80,
    )
    assert not bool(summary.loc[0, "strictly_stable_gene"])
    assert not bool(summary.loc[0, "compact_consensus_gene"])
    assert summary.loc[0, "stability_class"] == "exploratory_selected_gene"


def test_external_labels_are_not_referenced_by_threshold_selection(root: Path):
    source = (root / "scripts/select_task_c_threshold.py").read_text(encoding="utf-8").lower()
    assert "gse41258_samples.csv" not in source
    assert "external_validation_predictions" not in source
    assert 'columns=["gene_symbol"]' in source


@pytest.mark.full_data
def test_task_c_model_cards_separate_exact_and_transport_artifacts(root: Path):
    primary = json.loads(
        (root / "models/task_c_primary_full_signature_model_card.json").read_text()
    )
    transport = json.loads(
        (root / "models/task_c_cross_platform_transport_model_card.json").read_text()
    )
    assert primary["refitting_occurred"] is False
    assert transport["refitting_occurred"] is True
    assert primary["external_labels_used"] is False
    assert transport["external_labels_used"] is False
    assert set(transport["feature_genes"]) < set(primary["feature_genes"])
    assert set(transport["excluded_genes"]) == {"FOXQ1"}
    canonical = json.loads((root / "models/task_c_model_card.json").read_text())
    assert canonical == json.loads(
        (root / "models/task_c_primary_full_signature_model_card.json").read_text()
    )
    assert (root / "models/task_c_primary_model.joblib").stat().st_size > 0


@pytest.mark.full_data
def test_external_original_labels_and_patient_tissue_rows_are_preserved(root: Path):
    predictions = pd.read_csv(root / "results/metrics/external_validation_predictions.csv")
    primary = predictions.loc[
        predictions["representation"].eq("within_sample_percentile_rank")
        & predictions["evaluation_set"].eq("canonical_patient_tissue")
        & predictions["threshold_policy"].eq("primary_gse44076_locked")
    ]
    metadata = pd.read_csv(root / "data/metadata/gse41258_samples.csv").set_index(
        "geo_accession"
    )
    expected = metadata.loc[primary["sample_id"], "tissue_class"].to_numpy()
    np.testing.assert_array_equal(primary["y_true"].to_numpy(), expected)
    assert not primary.duplicated(["patient_id", "y_true"]).any()
    assert primary["patient_id"].nunique() == 190
    assert len(primary) == 233


@pytest.mark.full_data
def test_external_bootstrap_is_patient_clustered(root: Path):
    intervals = pd.read_csv(
        root / "results/metrics/external_patient_cluster_bootstrap_ci.csv"
    )
    primary = intervals.loc[
        intervals["representation"].eq("within_sample_percentile_rank")
        & intervals["evaluation_set"].eq("canonical_patient_tissue")
        & intervals["threshold_policy"].eq("primary_gse44076_locked")
    ]
    assert set(primary["cluster_unit"]) == {"patient"}
    assert set(primary["bootstrap_iterations"]) == {1000}


@pytest.mark.full_data
def test_task_specific_permutation_nulls_are_named_correctly(root: Path):
    frame = pd.read_csv(root / "results/metrics/group_preserving_permutation_tests.csv")
    nulls = frame.set_index("task")["null_hypothesis"].to_dict()
    assert "conditional null" in nulls["task_a_three_class"]
    assert "group-level null" in nulls["task_b_field_effect"]
    assert "paired null" in nulls["task_c_tumor_vs_adjacent"]
    assert set(frame["iterations"]) == {1000}


@pytest.mark.full_data
def test_manuscript_claims_and_model_genes_match_cards(root: Path):
    manuscript = (root / "manuscript/manuscript.md").read_text(encoding="utf-8")
    task_b = json.loads((root / "models/task_b_final_model_card.json").read_text())
    task_c = json.loads(
        (root / "models/task_c_primary_full_signature_model_card.json").read_text()
    )
    assert "[REF]" not in manuscript
    assert "GSE41258 cannot\nindependently validate healthy-versus-adjacent" in manuscript
    assert "Task C was secondary" in manuscript
    assert ", ".join(task_b["feature_genes"]) in manuscript
    assert ", ".join(task_c["feature_genes"]) in manuscript
    assert "arrayQualityMetrics completed for both cohorts" in manuscript
    for accession in ("GSE44076", "GSE41258"):
        report = root / "results/qc" / f"{accession}_array_quality_metrics" / "index.html"
        assert report.exists() and report.stat().st_size > 0


@pytest.mark.full_data
def test_serialized_threshold_provenance_has_no_external_tuning(root: Path):
    artifact = joblib.load(root / "models/task_c_cross_platform_transport_model.joblib")
    assert artifact["external_labels_used"] is False
    assert artifact["threshold"] == pytest.approx(0.41)
    assert "GSE44076" in artifact["threshold_origin"]


@pytest.mark.full_data
def test_required_publication_figure_renderings_are_nonempty(root: Path):
    stems = [
        "study_objectives_and_dataset_roles",
        "sample_count_flowchart",
        "raw_array_qc_summary",
        "GSE44076_pre_rma_density_all_arrays",
        "GSE44076_post_rma_density_all_arrays",
        "GSE41258_pre_rma_density_all_arrays",
        "GSE41258_post_rma_density_all_arrays",
        "unadjusted_adjacent_vs_healthy_volcano",
        "adjusted_adjacent_vs_healthy_volcano",
        "unadjusted_vs_adjusted_field_logfc",
        "paired_tumor_vs_adjacent_volcano",
        "final_field_signature_heatmap",
        "top_field_gene_group_distributions",
        "tissue_composition_by_group",
        "field_gene_evidence_tier_summary",
        "high_confidence_field_enrichment",
        "task_a_out_of_fold_confusion_matrix",
        "task_b_roc_pr_calibration",
        "task_b_panel_size_performance",
        "task_c_roc_pr_calibration",
        "task_c_panel_size_performance",
        "stable_feature_coefficient_plot_task_c",
        "task_c_oof_threshold_curve",
        "external_validation_roc_pr",
        "external_calibration",
        "external_confusion_matrix",
        "external_patient_structure_sensitivity",
        "signature_intersection",
        "biological_predictive_gene_intersection",
    ]
    for stem in stems:
        for suffix in (".png", ".pdf", ".svg"):
            path = root / "results/figures" / f"{stem}{suffix}"
            assert path.exists() and path.stat().st_size > 0, path


@pytest.mark.full_data
def test_runtime_artifacts_match_configured_counts_and_thresholds(root: Path):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    modeling = config["modeling"]
    decisions = pd.read_csv(root / "results/metrics/nested_panel_selection_decisions.csv")
    numeric_sizes = set(
        pd.to_numeric(decisions["panel_size"], errors="coerce").dropna().astype(int)
    )
    assert numeric_sizes == set(modeling["compact_panels"])
    permutations = pd.read_csv(
        root / "results/metrics/group_preserving_permutation_tests.csv"
    )
    assert set(permutations["iterations"]) == {modeling["permutation_iterations"]}
    nested_bootstrap = pd.read_csv(
        root / "results/metrics/nested_panel_oof_bootstrap_ci.csv"
    )
    assert set(nested_bootstrap["bootstrap_iterations"]) == {
        modeling["bootstrap_iterations"]
    }
    external_bootstrap = pd.read_csv(
        root / "results/metrics/external_patient_cluster_bootstrap_ci.csv"
    )
    completed = external_bootstrap["bootstrap_iterations"].gt(0)
    assert set(
        external_bootstrap.loc[completed, "bootstrap_iterations"].astype(int)
    ) == {
        config["external_validation"]["bootstrap_iterations"]
    }
    # Undefined predictive values (for example NPV when every prediction is
    # positive) have no estimable bootstrap draws and must remain explicitly NA.
    unestimable = external_bootstrap.loc[~completed]
    assert not unestimable.empty
    assert unestimable[
        ["mean", "median", "standard_deviation", "ci_lower", "ci_upper"]
    ].isna().all(axis=None)


@pytest.mark.full_data
def test_canonical_audit_artifacts_are_complete(root: Path):
    per_gene = pd.read_csv(
        root / "results/tables/covariate_adjusted_field_effect_concordance.csv"
    )
    assert len(per_gene) == 2 * 18490
    assert {
        "unadjusted_log2fc",
        "adjusted_log2fc",
        "effect_difference",
        "unadjusted_fdr",
        "adjusted_fdr",
        "effect_direction_agreement",
        "significance_retained",
    }.issubset(per_gene.columns)
    assert set(per_gene["adjusted_model"]) == {
        "U1_age_sex_adjusted",
        "U2_age_sex_location_adjusted",
    }

    composition = pd.read_csv(root / "results/tables/tissue_composition_sensitivity.csv")
    stress = pd.read_csv(root / "results/tables/preanalytical_stress_gene_audit.csv")
    threshold = pd.read_csv(root / "results/tables/task_c_threshold_selection.csv")
    task_c = pd.read_csv(root / "results/tables/task_c_final_signature.csv")
    assert len(composition) == 18490
    assert not stress.empty
    assert threshold["selected"].sum() == 1
    assert set(task_c["gene_symbol"]) == {"FOXQ1", "CEMIP", "ETV4"}


@pytest.mark.full_data
def test_biological_predictive_intersection_and_final_report(root: Path):
    intersection = pd.read_csv(
        root / "results/tables/biological_predictive_gene_intersection.csv"
    )
    assert not intersection.empty
    assert intersection["leading_edge_enrichment_gene"].any()
    assert intersection["task_b_predictive_gene"].sum() == 5
    assert intersection["task_c_predictive_gene"].sum() == 3
    report = (root / "reports/final_analysis_report.md").read_text(encoding="utf-8")
    required = [
        "Executive summary",
        "Adjusted field-effect analysis",
        "Task B confounding analysis",
        "Threshold locking",
        "External cohort structure",
        "Raw-vs-deposited sensitivity",
        "Final scientific conclusions",
    ]
    assert all(f"## {heading}" in report for heading in required)


@pytest.mark.full_data
def test_no_below_threshold_gene_is_called_strictly_stable(root: Path):
    config = yaml.safe_load((root / "config/analysis.yaml").read_text(encoding="utf-8"))
    for filename in ("task_b_final_signature.csv", "final_compact_signature.csv"):
        signature = pd.read_csv(root / "results/tables" / filename)
        called = signature["panel_label"].eq("strictly_stable_gene")
        assert signature.loc[called, "selection_frequency"].ge(
            config["modeling"]["stable_selection_frequency"]
        ).all()
        assert signature.loc[called, "sign_consistency"].ge(
            config["modeling"]["stable_sign_consistency"]
        ).all()
