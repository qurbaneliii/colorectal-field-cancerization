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

from src.data.provenance import expression_path, result_root
from src.modeling.evaluation import classification_metrics, prediction_probabilities
from src.modeling.nested_cv import modal_best_parameters
from src.modeling.pipelines import build_pipeline, selected_feature_names
from src.modeling.splitters import stratified_group_splits
from src.visualization.publication_figures import save_figure


TASKS = {
    "task_a_three_class": ["healthy", "adjacent_normal", "tumor"],
    "task_b_field_effect": ["healthy", "adjacent_normal"],
    "task_c_tumor_vs_adjacent": ["adjacent_normal", "tumor"],
}


def fit_panel_model(
    x_train: np.ndarray,
    y_train: np.ndarray,
    parameters: dict[str, object],
    seed: int,
    modeling: dict,
):
    model = build_pipeline("elastic_net", seed, modeling)
    locked_parameters = {
        **parameters,
        "variance_quantile__quantile": 0.0,
        "univariate__k": "all",
        "memory": None,
    }
    model.set_params(**locked_parameters)
    model.fit(x_train, y_train)
    return model


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

    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    for task, labels in TASKS.items():
        selected = metadata[
            metadata["inclusion_status"].eq("included")
            & metadata["tissue_class"].isin(labels)
        ].copy()
        samples = selected["geo_accession"].to_numpy()
        x = expression.loc[:, samples].T.to_numpy(dtype=np.float32)
        y = selected["tissue_class"].to_numpy()
        groups = selected["donor_or_patient_group"].astype(str).to_numpy()
        feature_names = expression.index.astype(str).tolist()
        for repeat, seed in enumerate(modeling["seeds"]):
            splits = stratified_group_splits(y, groups, int(modeling["outer_splits"]), int(seed))
            for outer_fold, (train, test) in enumerate(splits):
                row = nested_metrics[
                    nested_metrics["task"].eq(task)
                    & nested_metrics["model"].eq("elastic_net")
                    & nested_metrics["repeat"].eq(repeat)
                    & nested_metrics["outer_fold"].eq(outer_fold)
                ]
                if len(row) != 1:
                    raise AssertionError(f"Expected one nested parameter row for {task}/{repeat}/{outer_fold}")
                parameters = json.loads(row.iloc[0]["best_parameters"])
                full = build_pipeline("elastic_net", int(seed), modeling).set_params(
                    **parameters, memory=None
                )
                full.fit(x[train], y[train])
                selected_names = selected_feature_names(full, feature_names)
                coefficient = np.asarray(full.named_steps["model"].coef_)
                importance = np.max(np.abs(coefficient), axis=0)
                ranked_index = np.argsort(-importance, kind="stable")
                ranked = [selected_names[index] for index in ranked_index]
                ranked_importance = [float(importance[index]) for index in ranked_index]
                nonzero = [gene for gene, value in zip(ranked, ranked_importance, strict=True) if value > 1e-12]
                full_panel = nonzero or ranked
                panel_specs: list[tuple[str, list[str]]] = [
                    (str(size), ranked[: min(int(size), len(ranked))])
                    for size in modeling["compact_panels"]
                ]
                panel_specs.append(("full_stable_signature", full_panel))
                for panel_label, genes in panel_specs:
                    if not genes:
                        raise RuntimeError(f"Empty training-only panel for {task}/{repeat}/{outer_fold}")
                    positions = [feature_names.index(gene) for gene in genes]
                    panel_model = fit_panel_model(
                        x[train][:, positions], y[train], parameters, int(seed), modeling
                    )
                    predicted = panel_model.predict(x[test][:, positions])
                    classes, probability = prediction_probabilities(panel_model, x[test][:, positions])
                    measured = classification_metrics(y[test], predicted, probability, classes)
                    metric_rows.append(
                        {
                            "task": task,
                            "repeat": repeat,
                            "outer_fold": outer_fold,
                            "panel_size": panel_label,
                            "effective_gene_count": len(genes),
                            "n_train": len(train),
                            "n_test": len(test),
                            "analysis_provenance": provenance,
                            **measured,
                        }
                    )
                    for sample_position, sample_index in enumerate(test):
                        prediction = {
                            "task": task,
                            "repeat": repeat,
                            "outer_fold": outer_fold,
                            "panel_size": panel_label,
                            "sample_id": samples[sample_index],
                            "patient_id": groups[sample_index],
                            "y_true": y[sample_index],
                            "y_pred": predicted[sample_position],
                            "analysis_provenance": provenance,
                        }
                        prediction.update(
                            {
                                f"probability_{label}": float(probability[sample_position, index])
                                for index, label in enumerate(classes)
                            }
                        )
                        prediction_rows.append(prediction)
                    if panel_label != "full_stable_signature":
                        for rank, gene in enumerate(genes, start=1):
                            feature_rows.append(
                                {
                                    "task": task,
                                    "repeat": repeat,
                                    "outer_fold": outer_fold,
                                    "panel_size": panel_label,
                                    "rank": rank,
                                    "gene_symbol": gene,
                                    "training_only_absolute_coefficient": ranked_importance[
                                        ranked.index(gene)
                                    ],
                                    "analysis_provenance": provenance,
                                }
                            )

    fold_metrics = pd.DataFrame(metric_rows)
    predictions = pd.DataFrame(prediction_rows)
    selections = pd.DataFrame(feature_rows)
    fold_metrics.to_csv(metrics_root / "compact_panel_fold_metrics.csv", index=False)
    predictions.to_csv(metrics_root / "compact_panel_predictions.csv", index=False)
    selections.to_csv(metrics_root / "nested_cv_feature_selections.csv", index=False)

    summary = (
        fold_metrics.groupby(["task", "panel_size"], as_index=False)
        .agg(
            folds=("f1_macro", "size"),
            effective_gene_count_median=("effective_gene_count", "median"),
            macro_f1_mean=("f1_macro", "mean"),
            macro_f1_standard_deviation=("f1_macro", "std"),
            balanced_accuracy_mean=("balanced_accuracy", "mean"),
            roc_auc_mean=("roc_auc", "mean"),
            brier_score_mean=("brier_score", "mean"),
            log_loss_mean=("log_loss", "mean"),
        )
    )
    summary["analysis_provenance"] = provenance
    summary.to_csv(tables_root / "compact_panel_performance.csv", index=False)
    summary.to_csv(tables_root / "table_5_compact_panel_performance.csv", index=False)

    tolerance = float(modeling["compact_panel_practical_tolerance"])
    selected_sizes: dict[str, int] = {}
    selection_rows = []
    for task in TASKS:
        task_summary = summary[summary["task"].eq(task)].copy()
        full_mean = float(
            task_summary.loc[
                task_summary["panel_size"].eq("full_stable_signature"), "macro_f1_mean"
            ].iloc[0]
        )
        numeric = task_summary[task_summary["panel_size"].ne("full_stable_signature")].copy()
        numeric["numeric_size"] = numeric["panel_size"].astype(int)
        eligible = numeric[numeric["macro_f1_mean"].ge(full_mean - tolerance)].sort_values(
            "numeric_size"
        )
        if eligible.empty:
            chosen = int(numeric.sort_values("macro_f1_mean", ascending=False).iloc[0]["numeric_size"])
            reason = "No panel was within the practical tolerance; chose the best nested-CV panel."
        else:
            chosen = int(eligible.iloc[0]["numeric_size"])
            reason = f"Smallest panel within {tolerance:.3f} macro F1 of the full fold-specific signature."
        selected_sizes[task] = chosen
        selection_rows.append(
            {
                "task": task,
                "selected_panel_size": chosen,
                "full_signature_macro_f1": full_mean,
                "selected_panel_macro_f1": float(
                    numeric.loc[numeric["numeric_size"].eq(chosen), "macro_f1_mean"].iloc[0]
                ),
                "rule": reason,
                "analysis_provenance": provenance,
            }
        )
    selection_frame = pd.DataFrame(selection_rows)
    selection_frame.to_csv(tables_root / "compact_panel_selection.csv", index=False)

    task = "task_c_tumor_vs_adjacent"
    chosen = selected_sizes[task]
    chosen_rows = selections[
        selections["task"].eq(task) & selections["panel_size"].eq(str(chosen))
    ].copy()
    total_folds = len(modeling["seeds"]) * int(modeling["outer_splits"])
    chosen_rows["fold_key"] = (
        chosen_rows["repeat"].astype(str) + ":" + chosen_rows["outer_fold"].astype(str)
    )
    final_signature = (
        chosen_rows.groupby("gene_symbol", as_index=False)
        .agg(
            selected_fold_count=("fold_key", "nunique"),
            median_training_rank=("rank", "median"),
            median_absolute_coefficient=("training_only_absolute_coefficient", "median"),
        )
        .assign(selection_frequency=lambda frame: frame["selected_fold_count"] / total_folds)
        .sort_values(
            ["selection_frequency", "median_training_rank", "median_absolute_coefficient"],
            ascending=[False, True, False],
        )
        .head(chosen)
        .reset_index(drop=True)
    )
    final_signature.insert(0, "rank", np.arange(1, len(final_signature) + 1))
    final_signature.insert(0, "task", task)
    final_signature["analysis_provenance"] = provenance
    final_signature.to_csv(tables_root / "final_compact_signature.csv", index=False)
    final_signature.to_csv(tables_root / "table_4_final_signature.csv", index=False)

    selected_meta = metadata[
        metadata["inclusion_status"].eq("included")
        & metadata["tissue_class"].isin(TASKS[task])
    ].copy()
    genes = final_signature["gene_symbol"].tolist()
    parameters = modal_best_parameters(nested_metrics, task, "elastic_net")
    final_model = fit_panel_model(
        expression.loc[genes, selected_meta["geo_accession"]].T.to_numpy(dtype=np.float32),
        selected_meta["tissue_class"].to_numpy(),
        parameters,
        int(config["project"]["random_seed"]),
        modeling,
    )
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.strip()
    model_card = {
        "model_version": "1.0.0",
        "repository_commit": commit,
        "training_accession": "GSE44076",
        "training_sample_ids": selected_meta["geo_accession"].tolist(),
        "training_patient_group_ids": selected_meta["donor_or_patient_group"].astype(str).tolist(),
        "feature_genes": genes,
        "preprocessing_strategy": provenance,
        "hyperparameters": parameters,
        "threshold": float(config["external_validation"]["locked_threshold"]),
        "primary_external_representation": config["external_validation"]["primary_representation"],
        "random_seed": int(config["project"]["random_seed"]),
        "software_versions": {
            "python": sys.version.split()[0],
            "scikit_learn": sklearn.__version__,
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "analysis_provenance": provenance,
        "config_sha256": hashlib.sha256(config_text.encode("utf-8")).hexdigest(),
        "intended_use": "Retrospective transcriptomic biomarker discovery research.",
        "prohibited_claim": "Not a clinically validated diagnostic and not clinically ready.",
    }
    artifact = {"model": final_model, **model_card}
    model_path = models_root / "task_c_final_elastic_net.joblib"
    joblib.dump(artifact, model_path)
    reloaded = joblib.load(model_path)
    check_x = expression.loc[genes, selected_meta["geo_accession"].head(8)].T.to_numpy(
        dtype=np.float32
    )
    np.testing.assert_allclose(
        final_model.predict_proba(check_x), reloaded["model"].predict_proba(check_x)
    )
    (models_root / "task_c_final_model_card.json").write_text(
        json.dumps(model_card, indent=2), encoding="utf-8"
    )
    model_card_report = f"""# Task C model card

- Model: Elastic Net logistic regression
- Training accession: GSE44076 only
- Expression provenance: `{provenance}`
- Compact signature: {len(genes)} genes ({', '.join(genes)})
- Locked threshold: {model_card['threshold']}
- Primary external representation: `{model_card['primary_external_representation']}`
- Reload verification: PASS (probabilities reproduced exactly within floating-point tolerance)

This is a retrospective candidate transcriptomic biomarker-discovery model.
It is not a clinically validated diagnostic and is not clinically ready.
"""
    (reports_root / "task_c_model_card.md").write_text(model_card_report, encoding="utf-8")

    plot = summary.copy()
    plot = plot[plot["panel_size"].ne("full_stable_signature")]
    plot["panel_size_numeric"] = plot["panel_size"].astype(int)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for task_name, frame in plot.groupby("task"):
        frame = frame.sort_values("panel_size_numeric")
        ax.errorbar(
            frame["panel_size_numeric"],
            frame["macro_f1_mean"],
            yerr=frame["macro_f1_standard_deviation"],
            marker="o",
            capsize=3,
            label=task_name,
        )
    ax.set_xlabel("Training-only compact panel size")
    ax.set_ylabel("Outer-fold macro F1 (mean ± fold SD)")
    ax.set_title("Leakage-safe compact-panel performance")
    ax.legend(frameon=False)
    save_figure(fig, figures_root / "panel_size_performance_curve")

    rationale = f"""# Compact-panel selection rationale

Candidate panel construction was repeated independently inside every outer
training fold. Genes were ranked from the inner-selected Elastic Net fitted to
outer-training data only; each restricted model was then evaluated on the
untouched outer test fold. No global panel was used to estimate internal
performance.

{selection_frame.to_markdown(index=False)}

The final Task C gene list was derived only after panel-size performance was
estimated, by aggregating training-fold rankings across all repeat/fold keys.
Selected Task C genes: {', '.join(genes)}.
"""
    (reports_root / "compact_panel_selection_rationale.md").write_text(
        rationale, encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "provenance": provenance,
                "fold_metric_rows": len(fold_metrics),
                "selected_panel_sizes": selected_sizes,
                "task_c_genes": genes,
                "model_path": str(model_path.relative_to(ROOT)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
