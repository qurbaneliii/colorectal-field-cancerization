from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import f1_score
from sklearn.model_selection import GridSearchCV

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_compact_panel_analysis import fit_panel_model, rank_training_features
from scripts.run_task_b_confounding_sensitivity import demographic_matrix
from scripts.select_task_c_threshold import threshold_grid
from src.data.cross_platform import within_sample_percentile_rank
from src.data.expression import aggregate_probe_expression, read_geo_series_matrix
from src.data.provenance import file_sha256
from src.modeling.evaluation import classification_metrics, prediction_probabilities
from src.modeling.pipelines import build_pipeline, parameter_grid
from src.modeling.residualization import TrainingCompositionPCs, TrainingCovariateResidualizer
from src.modeling.splitters import stratified_group_splits


TASK_B_PARAMETERS = {
    "model__C": 0.2,
    "model__l1_ratio": 0.5,
    "model__class_weight": "balanced",
    "univariate__k": 20,
    "variance_quantile__quantile": 0.5,
}


def fit_compact(
    train_x: np.ndarray,
    train_y: np.ndarray,
    genes: list[str],
    parameters: dict,
    size: int,
    seed: int,
    modeling: dict,
) -> tuple:
    ranked, _ = rank_training_features(train_x, train_y, genes, parameters, seed, modeling)
    panel = ranked[:size]
    indices = [genes.index(gene) for gene in panel]
    fitted = fit_panel_model(train_x[:, indices], train_y, parameters, seed, modeling)
    return fitted, indices, panel


def fit_rank_policy(
    train_x: np.ndarray,
    train_y: np.ndarray,
    train_groups: np.ndarray,
    genes: list[str],
    seed: int,
    modeling: dict,
) -> tuple:
    cv = stratified_group_splits(train_y, train_groups, modeling["inner_splits"], seed)
    search = GridSearchCV(
        build_pipeline("elastic_net", seed, modeling).set_params(memory=None),
        parameter_grid("elastic_net", modeling),
        cv=cv,
        scoring="f1_macro",
        refit=False,
        error_score="raise",
        n_jobs=1,
    )
    search.fit(train_x, train_y)
    fitted, indices, panel = fit_compact(
        train_x, train_y, genes, search.best_params_, 3, seed, modeling
    )
    return fitted, indices, panel, search.best_params_


def positive_scores(fitted, matrix: np.ndarray, label: str) -> np.ndarray:
    classes, probabilities = prediction_probabilities(fitted, matrix)
    return probabilities[:, list(classes).index(label)]


def task_c_nested_thresholds(
    matrix: np.ndarray, metadata: pd.DataFrame, genes: list[str], config: dict
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit and score the complete rank/selection/threshold policy inside grouped splits."""
    modeling = config["modeling"]
    threshold_config = config["external_validation"]
    y = metadata["tissue_class"].to_numpy()
    groups = metadata["donor_or_patient_group"].to_numpy()
    sample_ids = metadata["geo_accession"].to_numpy()
    outcomes, decisions = [], []
    for repeat, seed in enumerate(modeling["seeds"]):
        for fold, (train, test) in enumerate(
            stratified_group_splits(y, groups, modeling["outer_splits"], seed)
        ):
            x_train, y_train, train_groups = matrix[train], y[train], groups[train]
            oof = np.full(len(train), np.nan)
            for inner, (fit, calibration) in enumerate(
                stratified_group_splits(
                    y_train, train_groups, modeling["inner_splits"], seed + 12000 + fold
                )
            ):
                fitted, columns, _, _ = fit_rank_policy(
                    x_train[fit], y_train[fit], train_groups[fit], genes,
                    seed + 13000 + inner + fold * 10, modeling,
                )
                oof[calibration] = positive_scores(fitted, x_train[calibration][:, columns], "tumor")
            if not np.isfinite(oof).all():
                raise ValueError("Incomplete inner threshold predictions")
            grid = threshold_grid(oof, y_train, threshold_config)
            threshold = float(grid.loc[grid["selected"], "threshold"].item())
            fitted, columns, panel, parameters = fit_rank_policy(
                x_train, y_train, train_groups, genes, seed + fold, modeling
            )
            scores = positive_scores(fitted, matrix[test][:, columns], "tumor")
            decisions.append(
                dict(
                    repeat=repeat, fold=fold, locked_threshold=threshold,
                    train_patients=len(set(train_groups)), test_patients=len(set(groups[test])),
                    genes=";".join(panel), hyperparameters=json.dumps(parameters, sort_keys=True),
                )
            )
            for index, score in zip(test, scores, strict=True):
                outcomes.append(
                    dict(
                        repeat=repeat, fold=fold, sample_id=sample_ids[index],
                        patient_id=groups[index], y_true=y[index], probability_tumor=score,
                        locked_threshold=threshold,
                    )
                )
    return pd.DataFrame(outcomes), pd.DataFrame(decisions)


def task_b_ablations(
    expression: pd.DataFrame, metadata: pd.DataFrame, config: dict
) -> tuple[pd.DataFrame, pd.DataFrame]:
    modeling = config["modeling"]
    y = metadata["tissue_class"].to_numpy()
    groups = metadata["donor_or_patient_group"].to_numpy()
    sample_ids = metadata["geo_accession"].to_numpy()
    genes = expression.index.tolist()
    values = expression[sample_ids].T.to_numpy(dtype=np.float32)
    demographics = demographic_matrix(metadata)
    markers = pd.read_csv(ROOT / "data/reference/mcpcounter_genes_b6eac73e.tsv", sep="\t")
    marker_genes = set(markers["HUGO symbols"])
    stress = set(pd.read_csv(ROOT / config["biology"]["stress_gene_set"], sep="\t")["gene_symbol"])
    score_values = np.column_stack([
        expression.loc[expression.index.intersection(group["HUGO symbols"]), sample_ids]
        .mean(axis=0).to_numpy()
        for _, group in markers.groupby("Cell population")
    ])
    if not np.isfinite(score_values).all():
        raise ValueError("Composition marker panel coverage is incomplete")
    excluded = [i for i, gene in enumerate(genes) if gene not in marker_genes | stress]
    eligible = [genes[i] for i in excluded]
    stress_indices = [i for i, gene in enumerate(genes) if gene in stress]
    outcomes, decisions = [], []
    for repeat, seed in enumerate(modeling["seeds"]):
        for fold, (train, test) in enumerate(
            stratified_group_splits(y, groups, modeling["outer_splits"], seed)
        ):
            pca = TrainingCompositionPCs(config["biology"]["composition_principal_components"])
            train_cov = np.column_stack([demographics[train], pca.fit_transform(score_values[train])])
            test_cov = np.column_stack([demographics[test], pca.transform(score_values[test])])
            stripped = values[:, excluded]
            adjustment = TrainingCovariateResidualizer().fit(stripped[train], train_cov)
            baseline_genes = [genes[i] for i in stress_indices]
            scenarios = {
                "fold_selected_expression": (values[train], values[test], genes),
                "fold_selected_without_markers_or_stress": (stripped[train], stripped[test], eligible),
                "fold_selected_residualized": (
                    adjustment.transform(stripped[train], train_cov),
                    adjustment.transform(stripped[test], test_cov),
                    eligible,
                ),
                "demographics_only": (demographics[train], demographics[test], None),
                "composition_proxy_only": (score_values[train], score_values[test], None),
                "stress_panel_only": (
                    values[train][:, stress_indices], values[test][:, stress_indices], None
                ),
            }
            for scenario, (train_x, test_x, feature_names) in scenarios.items():
                if feature_names is not None:
                    model, columns, panel = fit_compact(
                        train_x, y[train], feature_names, TASK_B_PARAMETERS, 5, seed + fold, modeling
                    )
                    test_x = test_x[:, columns]
                    decisions.append(dict(scenario=scenario, repeat=repeat, fold=fold, genes=panel))
                else:
                    model = fit_panel_model(train_x, y[train], TASK_B_PARAMETERS, seed + fold, modeling)
                    if scenario == "stress_panel_only":
                        decisions.append(
                            dict(scenario=scenario, repeat=repeat, fold=fold, genes=baseline_genes)
                        )
                scores = positive_scores(model, test_x, "healthy")
                for index, score in zip(test, scores, strict=True):
                    outcomes.append(
                        dict(
                            scenario=scenario, repeat=repeat, fold=fold, sample_id=sample_ids[index],
                            patient_id=groups[index], y_true=y[index], probability_healthy=score,
                        )
                    )
    return pd.DataFrame(outcomes), pd.DataFrame(decisions)


def summarize(outcomes: pd.DataFrame, positive: str, negative: str, threshold: float) -> dict:
    scores = outcomes[f"probability_{positive}"].to_numpy()
    labels = np.array([negative, positive])
    observed = outcomes["y_true"].to_numpy()
    predicted = np.where(scores >= threshold, positive, negative)
    return classification_metrics(observed, predicted, np.column_stack([1 - scores, scores]), labels)


def task_b_summary(outcomes: pd.DataFrame, seed: int) -> pd.DataFrame:
    scores = outcomes.pivot_table(
        index=["sample_id", "y_true"], columns="scenario", values="probability_healthy"
    ).reset_index()
    scenarios = sorted(outcomes["scenario"].unique())
    labels = scores["y_true"].to_numpy()
    calls = {
        scenario: np.where(scores[scenario].ge(0.5), "healthy", "adjacent_normal")
        for scenario in scenarios
    }
    rng = np.random.default_rng(seed)
    resamples = rng.integers(len(scores), size=(1000, len(scores)))
    rows = []
    for scenario in scenarios:
        metrics = summarize(
            scores.rename(columns={scenario: "probability_healthy"}), "healthy", "adjacent_normal", 0.5
        )
        baseline = calls["fold_selected_expression"]
        pairs = np.array([
            (
                f1_score(labels[idx], calls[scenario][idx], average="macro"),
                f1_score(labels[idx], baseline[idx], average="macro"),
            )
            for idx in resamples
        ])
        low, high = np.quantile(pairs[:, 0], [0.025, 0.975])
        delta_low, delta_high = np.quantile(pairs[:, 0] - pairs[:, 1], [0.025, 0.975])
        rows.append(
            dict(
                scenario=scenario, patients=len(scores), **metrics,
                macro_f1_ci_lower=low, macro_f1_ci_upper=high,
                macro_f1_delta_ci_lower=delta_low, macro_f1_delta_ci_upper=delta_high,
            )
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--series-matrix", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "results/ai_audit")
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    config_path = ROOT / "config/analysis.yaml"
    metadata_path = ROOT / "data/metadata/gse44076_samples.csv"
    primary_map = ROOT / "data/metadata/GSE44076_probe_gene_mapping_raw_cel_rma.csv"
    external_map = ROOT / "data/metadata/GSE41258_probe_gene_mapping_raw_cel_rma.csv"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    metadata = pd.read_csv(metadata_path).query("inclusion_status == 'included'")
    probe_expression = read_geo_series_matrix(args.series_matrix)
    mapping = pd.read_csv(primary_map).query("mapping_status == 'mapped_unique'")
    expression = aggregate_probe_expression(probe_expression, mapping)
    field = metadata[metadata["tissue_class"].isin(["healthy", "adjacent_normal"])]
    b_predictions, b_panels = task_b_ablations(expression, field, config)
    seed = config["project"]["random_seed"]
    task_b_summary(b_predictions, seed).to_csv(output / "task_b_development_ablations.csv", index=False)
    frequency = []
    for scenario, rows in b_panels.groupby("scenario"):
        counts = Counter(gene for panel in rows["genes"] for gene in panel)
        frequency.extend(
            dict(scenario=scenario, gene=gene, folds=count, total_folds=len(rows))
            for gene, count in sorted(counts.items())
        )
    pd.DataFrame(frequency).to_csv(output / "task_b_development_stability.csv", index=False)
    cancer = metadata[metadata["tissue_class"].isin(["adjacent_normal", "tumor"])]
    external_genes = pd.read_csv(external_map).query("mapping_status == 'mapped_unique'")["gene_symbol"]
    common = expression.index.intersection(external_genes).tolist()
    ranks = within_sample_percentile_rank(expression.loc[common, cancer["geo_accession"]].T.to_numpy())
    c_predictions, c_thresholds = task_c_nested_thresholds(ranks, cancer, common, config)
    c_thresholds.to_csv(output / "task_c_nested_threshold_decisions.csv", index=False)
    c_metrics = []
    for (repeat, fold), predictions in c_predictions.groupby(["repeat", "fold"]):
        for policy, threshold in (
            ("inner_locked", float(predictions["locked_threshold"].iloc[0])),
            ("default", config["external_validation"]["default_threshold"]),
        ):
            c_metrics.append(
                dict(repeat=repeat, fold=fold, policy=policy, threshold=threshold,
                     **summarize(predictions, "tumor", "adjacent_normal", threshold))
            )
    pd.DataFrame(c_metrics).to_csv(output / "task_c_nested_threshold_metrics.csv", index=False)
    provenance = {
        "expression_provenance": "geo_deposited_series_matrix",
        "accession": "GSE44076",
        "series_matrix_sha256": file_sha256(args.series_matrix),
        "probe_rows": len(probe_expression),
        "gene_rows": len(expression),
        "common_universe_genes": len(common),
        "external_expression_or_labels_loaded": False,
        "sha256": {
            path.relative_to(ROOT).as_posix(): file_sha256(path)
            for path in (
                config_path, metadata_path, primary_map, external_map, Path(__file__).resolve(),
                ROOT / "data/reference/mcpcounter_genes_b6eac73e.tsv",
                ROOT / config["biology"]["stress_gene_set"],
            )
        },
    }
    (output / "development_experiment_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
