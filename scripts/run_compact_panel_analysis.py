from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.provenance import expression_path, result_root, source_files_sha256
from src.modeling.evaluation import (
    classification_metrics,
    grouped_bootstrap_metrics,
    prediction_probabilities,
)
from src.modeling.nested_cv import modal_best_parameters
from src.modeling.pipelines import build_pipeline, selected_feature_names
from src.modeling.splitters import stratified_group_splits
from src.visualization.publication_figures import save_figure


TASKS = {
    "task_a_three_class": ["healthy", "adjacent_normal", "tumor"],
    "task_b_field_effect": ["healthy", "adjacent_normal"],
    "task_c_tumor_vs_adjacent": ["adjacent_normal", "tumor"],
}

MODEL_SOURCE_FILES = [
    "config/analysis.yaml",
    "scripts/run_compact_panel_analysis.py",
    "src/modeling/evaluation.py",
    "src/modeling/nested_cv.py",
    "src/modeling/pipelines.py",
    "src/modeling/splitters.py",
    "src/modeling/stability.py",
]


def fit_panel_model(
    x_train: np.ndarray,
    y_train: np.ndarray,
    parameters: dict[str, object],
    seed: int,
    modeling: dict,
):
    model = build_pipeline("elastic_net", seed, modeling)
    model.set_params(
        **{
            **parameters,
            "variance_quantile__quantile": 0.0,
            "univariate__k": "all",
            "memory": None,
        }
    )
    model.fit(x_train, y_train)
    return model


def rank_training_features(
    x_train: np.ndarray,
    y_train: np.ndarray,
    feature_names: list[str],
    parameters: dict[str, object],
    seed: int,
    modeling: dict,
) -> tuple[list[str], dict[str, float]]:
    full = build_pipeline("elastic_net", seed, modeling).set_params(
        **parameters, memory=None
    )
    full.fit(x_train, y_train)
    selected = selected_feature_names(full, feature_names)
    coefficient = np.asarray(full.named_steps["model"].coef_)
    importance = np.max(np.abs(coefficient), axis=0)
    order = np.argsort(-importance, kind="stable")
    ranked = [selected[index] for index in order]
    importance_by_gene = {
        selected[index]: float(importance[index]) for index in order
    }
    nonzero = [gene for gene in ranked if importance_by_gene[gene] > 1e-12]
    return nonzero or ranked, importance_by_gene


def panel_specs(ranked: list[str], sizes: list[int]) -> dict[str, list[str]]:
    specs = {str(size): ranked[: min(size, len(ranked))] for size in sizes}
    specs["full_nonzero"] = ranked
    if any(not genes for genes in specs.values()):
        raise RuntimeError("Training-only ranking produced an empty panel")
    return specs


def inner_panel_scores(
    x: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    feature_names: list[str],
    parameters: dict[str, object],
    seed: int,
    modeling: dict,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    splits = stratified_group_splits(
        y, groups, int(modeling["inner_splits"]), seed + 7000
    )
    positions = {gene: index for index, gene in enumerate(feature_names)}
    for inner_fold, (train, validation) in enumerate(splits):
        ranked, _ = rank_training_features(
            x[train], y[train], feature_names, parameters, seed + inner_fold, modeling
        )
        for label, genes in panel_specs(ranked, modeling["compact_panels"]).items():
            columns = [positions[gene] for gene in genes]
            model = fit_panel_model(
                x[train][:, columns],
                y[train],
                parameters,
                seed + inner_fold,
                modeling,
            )
            predicted = model.predict(x[validation][:, columns])
            classes, probability = prediction_probabilities(
                model, x[validation][:, columns]
            )
            measured = classification_metrics(
                y[validation], predicted, probability, classes
            )
            rows.append(
                {
                    "inner_fold": inner_fold,
                    "panel_size": label,
                    "effective_gene_count": len(genes),
                    "brier_score": measured.get(
                        "brier_score", measured.get("brier_score_multiclass", np.nan)
                    ),
                    **measured,
                }
            )
    return pd.DataFrame(rows)


def choose_panel(inner: pd.DataFrame, modeling: dict) -> tuple[str, pd.DataFrame]:
    summary = (
        inner.groupby("panel_size", as_index=False)
        .agg(
            inner_folds=("inner_fold", "nunique"),
            effective_gene_count=("effective_gene_count", "median"),
            inner_macro_f1=("f1_macro", "mean"),
            inner_balanced_accuracy=("balanced_accuracy", "mean"),
            inner_brier_score=("brier_score", "mean"),
            inner_log_loss=("log_loss", "mean"),
        )
    )
    best_f1 = float(summary["inner_macro_f1"].max())
    best_brier = float(summary["inner_brier_score"].min())
    eligible = summary[
        summary["inner_macro_f1"].ge(
            best_f1 - float(modeling["compact_panel_practical_tolerance"])
        )
        & summary["inner_balanced_accuracy"].ge(
            float(modeling["compact_panel_min_balanced_accuracy"])
        )
        & summary["inner_brier_score"].le(
            best_brier + float(modeling["compact_panel_max_brier_degradation"])
        )
    ].copy()
    if eligible.empty:
        eligible = summary.nlargest(1, "inner_macro_f1").copy()
    eligible["selection_order"] = np.where(
        eligible["panel_size"].eq("full_nonzero"),
        np.inf,
        pd.to_numeric(eligible["panel_size"], errors="coerce"),
    )
    chosen = str(eligible.sort_values("selection_order").iloc[0]["panel_size"])
    summary["selected"] = summary["panel_size"].eq(chosen)
    summary["selection_rule"] = (
        "smallest panel within configured inner-CV macro-F1 tolerance, "
        "balanced-accuracy floor, and Brier degradation limit"
    )
    return chosen, summary


def average_repeated_predictions(frame: pd.DataFrame) -> pd.DataFrame:
    probability_columns = sorted(
        column
        for column in frame
        if column.startswith("probability_") and frame[column].notna().any()
    )
    averaged = frame.groupby("sample_id", as_index=False).agg(
        {
            "patient_id": "first",
            "y_true": "first",
            **{column: "mean" for column in probability_columns},
        }
    )
    classes = np.asarray(
        [column.removeprefix("probability_") for column in probability_columns]
    )
    averaged["y_pred"] = classes[
        np.argmax(averaged[probability_columns].to_numpy(), axis=1)
    ]
    return averaged


def write_model_artifact(
    task: str,
    genes: list[str],
    final_size: str,
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    nested_metrics: pd.DataFrame,
    config: dict,
    config_text: str,
    provenance: str,
    models_root: Path,
) -> dict[str, object]:
    modeling = config["modeling"]
    selected = metadata[
        metadata["inclusion_status"].eq("included")
        & metadata["tissue_class"].isin(TASKS[task])
    ].copy()
    parameters = modal_best_parameters(nested_metrics, task, "elastic_net")
    x = expression.loc[genes, selected["geo_accession"]].T.to_numpy(dtype=np.float32)
    model = fit_panel_model(
        x,
        selected["tissue_class"].to_numpy(),
        parameters,
        int(config["project"]["random_seed"]),
        modeling,
    )
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    is_task_b = task == "task_b_field_effect"
    card: dict[str, object] = {
        "model_version": "2.0.0",
        "task": task,
        "repository_commit": commit,
        "training_accession": "GSE44076",
        "training_representation": provenance,
        "training_sample_ids": selected["geo_accession"].tolist(),
        "training_patient_group_ids": selected["donor_or_patient_group"].astype(str).tolist(),
        "feature_genes": genes,
        "excluded_genes": [],
        "exclusion_reason": "not_applicable",
        "refitting_occurred": False,
        "external_labels_used": False,
        "final_panel_size": final_size,
        "panel_size_origin": "GSE44076 grouped inner-CV only",
        "hyperparameters": parameters,
        "threshold": None if is_task_b else float(config["external_validation"]["default_threshold"]),
        "threshold_origin": (
            "not_applicable_default_classifier_decision"
            if is_task_b
            else "temporary_default_pending_GSE44076_OOF_threshold_lock"
        ),
        "random_seed": int(config["project"]["random_seed"]),
        "analysis_provenance": provenance,
        "config_sha256": hashlib.sha256(config_text.encode()).hexdigest(),
        "analysis_source_files": MODEL_SOURCE_FILES,
        "analysis_source_sha256": source_files_sha256(ROOT, MODEL_SOURCE_FILES),
        "software_versions": {
            "python": sys.version.split()[0],
            "scikit_learn": sklearn.__version__,
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "intended_use": (
            "Internally validated retrospective field-effect biomarker discovery."
            if is_task_b
            else "Secondary retrospective tumor-normal cross-platform biomarker discovery."
        ),
        "validation_boundary": (
            "Internally validated only; no independent healthy-versus-adjacent cohort."
            if is_task_b
            else "External evaluation is performed only after a label-independent common-gene transport refit."
        ),
        "clinical_readiness": "not clinically ready",
    }
    artifact = {"model": model, **card}
    if is_task_b:
        path = models_root / "task_b_final_elastic_net.joblib"
        card_path = models_root / "task_b_final_model_card.json"
    else:
        path = models_root / "task_c_primary_full_signature_model.joblib"
        card_path = models_root / "task_c_primary_full_signature_model_card.json"
    joblib.dump(artifact, path)
    reloaded = joblib.load(path)
    np.testing.assert_allclose(model.predict_proba(x[:8]), reloaded["model"].predict_proba(x[:8]))
    card_path.write_text(json.dumps(card, indent=2), encoding="utf-8")
    if is_task_b:
        (models_root / "task_b_model_card.json").write_text(
            json.dumps(card, indent=2), encoding="utf-8"
        )
    if not is_task_b:
        joblib.dump(artifact, models_root / "task_c_primary_model.joblib")
        (models_root / "task_c_model_card.json").write_text(
            json.dumps(card, indent=2), encoding="utf-8"
        )
        joblib.dump(artifact, models_root / "task_c_final_elastic_net.joblib")
        (models_root / "task_c_final_model_card.json").write_text(
            json.dumps(card, indent=2), encoding="utf-8"
        )
    return card


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--provenance", choices=["raw_cel_rma", "geo_deposited_series_matrix"]
    )
    args = parser.parse_args()
    os.chdir(ROOT)
    config_text = (ROOT / "config/analysis.yaml").read_text(encoding="utf-8")
    config = yaml.safe_load(config_text)
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    modeling = config["modeling"]
    provenance = args.provenance or modeling["expression_provenance"]
    expression = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE44076", "gene", provenance)
    ).set_index("gene_symbol")
    metadata = pd.read_csv(
        ROOT / paths["metadata"] / "gse44076_samples.csv", dtype={"patient_id": str}
    )
    metrics_root = result_root(ROOT / paths["metrics"], provenance)
    tables_root = result_root(ROOT / paths["tables"], provenance)
    figures_root = result_root(ROOT / paths["figures"], provenance)
    reports_root = result_root(ROOT / paths["reports"], provenance)
    models_root = result_root(ROOT / paths["models"], provenance)
    for directory in (metrics_root, tables_root, figures_root, reports_root, models_root):
        directory.mkdir(parents=True, exist_ok=True)
    nested_metrics = pd.read_csv(metrics_root / "nested_cv_metrics.csv")
    feature_names = expression.index.astype(str).tolist()
    feature_position = {gene: index for index, gene in enumerate(feature_names)}

    decision_frames: list[pd.DataFrame] = []
    outer_metric_rows: list[dict[str, object]] = []
    outer_prediction_rows: list[dict[str, object]] = []
    selected_prediction_rows: list[dict[str, object]] = []
    selected_gene_rows: list[dict[str, object]] = []
    final_panels: dict[str, tuple[str, list[str], dict[str, object]]] = {}

    for task, labels in TASKS.items():
        selected = metadata[
            metadata["inclusion_status"].eq("included")
            & metadata["tissue_class"].isin(labels)
        ].copy()
        samples = selected["geo_accession"].to_numpy()
        x = expression.loc[:, samples].T.to_numpy(dtype=np.float32)
        y = selected["tissue_class"].to_numpy()
        groups = selected["donor_or_patient_group"].astype(str).to_numpy()
        for repeat, seed in enumerate(modeling["seeds"]):
            outer_splits = stratified_group_splits(
                y, groups, int(modeling["outer_splits"]), int(seed)
            )
            for outer_fold, (train, test) in enumerate(outer_splits):
                metric_row = nested_metrics[
                    nested_metrics["task"].eq(task)
                    & nested_metrics["model"].eq("elastic_net")
                    & nested_metrics["repeat"].eq(repeat)
                    & nested_metrics["outer_fold"].eq(outer_fold)
                ]
                if len(metric_row) != 1:
                    raise AssertionError(
                        f"Expected one nested parameter row for {task}/{repeat}/{outer_fold}"
                    )
                parameters = json.loads(metric_row.iloc[0]["best_parameters"])
                inner = inner_panel_scores(
                    x[train],
                    y[train],
                    groups[train],
                    feature_names,
                    parameters,
                    int(seed) + outer_fold,
                    modeling,
                )
                chosen, decision = choose_panel(inner, modeling)
                decision.insert(0, "outer_fold", outer_fold)
                decision.insert(0, "repeat", repeat)
                decision.insert(0, "task", task)
                decision["analysis_provenance"] = provenance
                decision_frames.append(decision)

                ranked, importance = rank_training_features(
                    x[train], y[train], feature_names, parameters, int(seed), modeling
                )
                specs = panel_specs(ranked, modeling["compact_panels"])
                selected_genes = specs[chosen]
                for rank, gene in enumerate(selected_genes, start=1):
                    selected_gene_rows.append(
                        {
                            "task": task,
                            "repeat": repeat,
                            "outer_fold": outer_fold,
                            "selected_panel_size": chosen,
                            "rank": rank,
                            "gene_symbol": gene,
                            "training_only_absolute_coefficient": importance[gene],
                            "analysis_provenance": provenance,
                        }
                    )
                for label, genes in {**specs, "nested_selected_policy": selected_genes}.items():
                    columns = [feature_position[gene] for gene in genes]
                    model = fit_panel_model(
                        x[train][:, columns], y[train], parameters, int(seed), modeling
                    )
                    predicted = model.predict(x[test][:, columns])
                    classes, probability = prediction_probabilities(model, x[test][:, columns])
                    measured = classification_metrics(y[test], predicted, probability, classes)
                    outer_metric_rows.append(
                        {
                            "task": task,
                            "repeat": repeat,
                            "outer_fold": outer_fold,
                            "panel_size": label,
                            "effective_gene_count": len(genes),
                            "n_train": len(train),
                            "n_test": len(test),
                            "analysis_provenance": provenance,
                            **measured,
                        }
                    )
                    for row_position, sample_index in enumerate(test):
                        row = {
                            "task": task,
                            "repeat": repeat,
                            "outer_fold": outer_fold,
                            "panel_size": label,
                            "sample_id": samples[sample_index],
                            "patient_id": groups[sample_index],
                            "y_true": y[sample_index],
                            "y_pred": predicted[row_position],
                            "analysis_provenance": provenance,
                        }
                        row.update(
                            {
                                f"probability_{class_label}": float(
                                    probability[row_position, class_index]
                                )
                                for class_index, class_label in enumerate(classes)
                            }
                        )
                        outer_prediction_rows.append(row)
                        if label == "nested_selected_policy":
                            selected_prediction_rows.append(row)

        parameters = modal_best_parameters(nested_metrics, task, "elastic_net")
        full_inner = inner_panel_scores(
            x,
            y,
            groups,
            feature_names,
            parameters,
            int(config["project"]["random_seed"]) + 9000,
            modeling,
        )
        final_size, final_decision = choose_panel(full_inner, modeling)
        ranked, importance = rank_training_features(
            x,
            y,
            feature_names,
            parameters,
            int(config["project"]["random_seed"]),
            modeling,
        )
        final_genes = panel_specs(ranked, modeling["compact_panels"])[final_size]
        final_panels[task] = (final_size, final_genes, importance)
        final_decision.insert(0, "outer_fold", "full_development_selection")
        final_decision.insert(0, "repeat", "full_development_selection")
        final_decision.insert(0, "task", task)
        final_decision["analysis_provenance"] = provenance
        decision_frames.append(final_decision)

    decisions = pd.concat(decision_frames, ignore_index=True)
    outer_metrics = pd.DataFrame(outer_metric_rows)
    outer_predictions = pd.DataFrame(outer_prediction_rows)
    selected_predictions = pd.DataFrame(selected_prediction_rows)
    selections = pd.DataFrame(selected_gene_rows)
    decisions.to_csv(metrics_root / "nested_panel_selection_decisions.csv", index=False)
    outer_metrics.to_csv(metrics_root / "compact_panel_fold_metrics.csv", index=False)
    outer_predictions.to_csv(metrics_root / "compact_panel_predictions.csv", index=False)
    selected_predictions.to_csv(metrics_root / "nested_panel_outer_predictions.csv", index=False)
    selections.to_csv(metrics_root / "nested_cv_feature_selections.csv", index=False)

    summary = (
        outer_metrics.groupby(["task", "panel_size"], as_index=False)
        .agg(
            outer_folds=("f1_macro", "size"),
            effective_gene_count_median=("effective_gene_count", "median"),
            outer_fold_macro_f1_mean=("f1_macro", "mean"),
            outer_fold_macro_f1_sd=("f1_macro", "std"),
            outer_fold_balanced_accuracy_mean=("balanced_accuracy", "mean"),
            outer_fold_roc_auc_mean=("roc_auc", "mean"),
            outer_fold_pr_auc_mean=("pr_auc", "mean"),
            outer_fold_brier_mean=("brier_score", "mean"),
            outer_fold_log_loss_mean=("log_loss", "mean"),
        )
    )
    summary["analysis_provenance"] = provenance
    summary.to_csv(tables_root / "nested_compact_panel_summary.csv", index=False)
    summary.to_csv(tables_root / "compact_panel_performance.csv", index=False)
    summary.to_csv(tables_root / "table_5_compact_panel_performance.csv", index=False)
    summary[summary["task"].eq("task_b_field_effect")].to_csv(
        tables_root / "task_b_compact_panel_performance.csv", index=False
    )

    stability = pd.read_csv(tables_root / "feature_stability.csv")
    external_genes = set(
        pd.read_parquet(
            expression_path(ROOT / paths["processed"], "GSE41258", "gene", provenance),
            columns=["gene_symbol"],
        )["gene_symbol"]
    )
    signature_tables: dict[str, pd.DataFrame] = {}
    for task in ("task_b_field_effect", "task_c_tumor_vs_adjacent"):
        final_size, genes, importance = final_panels[task]
        task_selections = selections[selections["task"].eq(task)].copy()
        total_folds = len(modeling["seeds"]) * int(modeling["outer_splits"])
        frequency = (
            task_selections.assign(
                fold_key=lambda frame: frame["repeat"].astype(str)
                + ":"
                + frame["outer_fold"].astype(str)
            )
            .groupby("gene_symbol")
            .agg(
                compact_selected_fold_count=("fold_key", "nunique"),
                compact_median_training_rank=("rank", "median"),
            )
        )
        task_stability = stability[stability["task"].eq(task)].copy()
        if task_stability["class"].nunique() != 1:
            raise AssertionError(f"Expected one binary coefficient class for {task}")
        task_stability = task_stability.set_index("gene_symbol")
        table = pd.DataFrame({"gene_symbol": genes})
        table["rank"] = np.arange(1, len(table) + 1)
        table["final_full_data_absolute_coefficient"] = table["gene_symbol"].map(importance)
        table["compact_selected_fold_count"] = table["gene_symbol"].map(
            frequency["compact_selected_fold_count"]
        ).fillna(0).astype(int)
        table["compact_selection_frequency"] = table["compact_selected_fold_count"] / total_folds
        for column in [
            "selection_frequency",
            "sign_consistency",
            "median_coefficient",
            "coefficient_standard_deviation",
            "repeat_coverage",
            "strictly_stable_gene",
            "stability_class",
        ]:
            table[column] = table["gene_symbol"].map(task_stability[column])
        table["compact_consensus_gene"] = table["compact_selection_frequency"].ge(0.50)
        table["panel_label"] = np.where(
            table["strictly_stable_gene"].fillna(False),
            "strictly_stable_gene",
            np.where(
                table["compact_consensus_gene"],
                "compact_consensus_gene",
                "exploratory_panel_gene",
            ),
        )
        table["final_panel_size"] = final_size
        table["analysis_provenance"] = provenance
        if task == "task_c_tumor_vs_adjacent":
            table["present_on_gpl96"] = table["gene_symbol"].isin(external_genes)
        signature_tables[task] = table

    task_b_signature = signature_tables["task_b_field_effect"]
    task_c_signature = signature_tables["task_c_tumor_vs_adjacent"]
    task_b_signature.to_csv(tables_root / "task_b_final_signature.csv", index=False)
    task_c_signature.to_csv(tables_root / "final_compact_signature.csv", index=False)
    task_c_signature.to_csv(tables_root / "task_c_final_signature.csv", index=False)
    task_c_signature.to_csv(tables_root / "table_4_final_signature.csv", index=False)

    task_c_internal = summary[summary["task"].eq("task_c_tumor_vs_adjacent")].copy()
    task_c_internal["candidate_panel"] = task_c_internal["panel_size"]
    task_c_internal["stability_requirement"] = (
        "final genes reported individually; no below-threshold gene called strictly stable"
    )
    task_c_internal["cross_platform_coverage"] = task_c_internal["panel_size"].map(
        lambda label: (
            float(task_c_signature["present_on_gpl96"].mean())
            if label == final_panels["task_c_tumor_vs_adjacent"][0]
            else np.nan
        )
    )
    task_c_internal["external_performance_role"] = (
        "reported after locking; never used for panel selection"
    )
    task_c_internal.to_csv(tables_root / "task_c_panel_decision_matrix.csv", index=False)

    uncertainty_rows: list[dict[str, object]] = []
    for task in ("task_b_field_effect", "task_c_tumor_vs_adjacent"):
        frame = selected_predictions[selected_predictions["task"].eq(task)]
        averaged = average_repeated_predictions(frame)
        metrics = [
            "roc_auc",
            "pr_auc",
            "f1_macro",
            "f1_weighted",
            "balanced_accuracy",
            "brier_score",
            "log_loss",
        ]
        for row in grouped_bootstrap_metrics(
            averaged,
            metrics,
            int(modeling["bootstrap_iterations"]),
            int(config["project"]["random_seed"]),
        ):
            uncertainty_rows.append(
                {
                    "task": task,
                    "estimand": "patient_bootstrap_of_aggregated_repeated_oof",
                    **row,
                    "analysis_provenance": provenance,
                }
            )
    uncertainty = pd.DataFrame(uncertainty_rows)
    uncertainty.to_csv(metrics_root / "nested_panel_oof_bootstrap_ci.csv", index=False)

    cards = {}
    for task in ("task_b_field_effect", "task_c_tumor_vs_adjacent"):
        size, genes, _ = final_panels[task]
        cards[task] = write_model_artifact(
            task,
            genes,
            size,
            expression,
            metadata,
            nested_metrics,
            config,
            config_text,
            provenance,
            models_root,
        )

    for task, stem, title in [
        ("task_b_field_effect", "task_b_panel_size_performance", "Task B field-effect panel size"),
        ("task_c_tumor_vs_adjacent", "task_c_panel_size_performance", "Task C tumor-normal panel size"),
    ]:
        plot = summary[
            summary["task"].eq(task)
            & ~summary["panel_size"].isin(["full_nonzero", "nested_selected_policy"])
        ].copy()
        plot["panel_size_numeric"] = plot["panel_size"].astype(int)
        plot = plot.sort_values("panel_size_numeric")
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        ax.errorbar(
            plot["panel_size_numeric"],
            plot["outer_fold_macro_f1_mean"],
            yerr=plot["outer_fold_macro_f1_sd"],
            marker="o",
            capsize=3,
        )
        ax.set(
            xlabel="Training-only panel size",
            ylabel="Outer-fold macro F1 (mean and fold SD)",
            title=f"{title}; final size chosen only within grouped inner CV",
        )
        save_figure(fig, figures_root / stem)

    task_b_rows = outer_metrics[
        outer_metrics["task"].eq("task_b_field_effect")
        & outer_metrics["panel_size"].eq("nested_selected_policy")
    ]
    task_b_ci = uncertainty[uncertainty["task"].eq("task_b_field_effect")]
    task_b_report = f"""# Task B field-effect model report

Task B (healthy versus adjacent-normal) is the primary predictive endpoint.
Every outer split preserved donor/patient groups. Within each outer-training
partition, feature ranking and panel-size evaluation were repeated in grouped
inner folds; the smallest panel satisfying the configured macro-F1 tolerance,
balanced-accuracy floor, and calibration-loss limit was selected before the
outer test fold was evaluated once.

- Deterministic outer repeats: {len(modeling['seeds'])}
- Outer folds per repeat: {modeling['outer_splits']}
- Mean nested-policy outer-fold macro-F1: {task_b_rows['f1_macro'].mean():.4f}
- Final full-development panel: {', '.join(task_b_signature['gene_symbol'])}
- Panel-size origin: GSE44076 grouped inner CV only
- Independent healthy-versus-adjacent validation: unavailable
- Intended use: retrospective biomarker discovery; not clinically ready

Patient/donor-bootstrap intervals of aggregated repeated OOF predictions:

{task_b_ci.to_markdown(index=False)}
"""
    (reports_root / "task_b_field_effect_model_report.md").write_text(
        task_b_report, encoding="utf-8"
    )
    task_b_card_report = f"""# Task B model card

- Objective: healthy versus adjacent-normal field effect
- Model family: Elastic Net logistic regression
- Training accession: GSE44076 only
- Expression provenance: `{provenance}`
- Final panel: {', '.join(task_b_signature['gene_symbol'])}
- Panel size selected with grouped inner CV only
- Internal validation: repeated patient/donor-grouped nested CV
- External validation: none available for the field-effect task
- Status: internally validated only; not clinically ready
"""
    (reports_root / "task_b_model_card.md").write_text(
        task_b_card_report, encoding="utf-8"
    )
    task_c_report = f"""# Task C model card

- Objective: secondary tumor versus adjacent-normal model
- Model family: Elastic Net logistic regression
- Training accession: GSE44076 only
- Expression provenance: `{provenance}`
- Locked primary signature: {', '.join(task_c_signature['gene_symbol'])}
- Per-gene stability status: `results/tables/final_compact_signature.csv`
- Threshold: selected later from GSE44076 repeated OOF rank-transport probabilities
- Exact primary artifact: `models/task_c_primary_full_signature_model.joblib`
- Cross-platform transport artifact: created separately during external validation
- Status: retrospective biomarker-discovery model; not clinically ready
"""
    (reports_root / "task_c_model_card.md").write_text(task_c_report, encoding="utf-8")

    panel_rationale = f"""# Task C panel selection rationale

Task C panel size was treated as an inner-CV hyperparameter. For each outer
fold, rankings were constructed from inner-training data, candidate sizes
{modeling['compact_panels']} plus the fold-specific full nonzero panel were
evaluated on grouped inner validation folds, and the configured multi-criterion
rule selected a size without using the outer test fold. External labels were
not read by this workflow.

The final full-development panel size was
`{final_panels['task_c_tumor_vs_adjacent'][0]}` and the panel is
{', '.join(task_c_signature['gene_symbol'])}. A gene is called strictly stable
only when the general Elastic Net selection frequency is at least
{modeling['stable_selection_frequency']} and sign consistency is at least
{modeling['stable_sign_consistency']}; other retained genes are named compact
consensus or exploratory panel genes in the signature table.
"""
    (reports_root / "task_c_panel_selection_rationale.md").write_text(
        panel_rationale, encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "provenance": provenance,
                "outer_metric_rows": len(outer_metrics),
                "nested_decision_rows": len(decisions),
                "final_panels": {
                    task: {"size": value[0], "genes": value[1]}
                    for task, value in final_panels.items()
                },
                "model_cards": {task: card["feature_genes"] for task, card in cards.items()},
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
