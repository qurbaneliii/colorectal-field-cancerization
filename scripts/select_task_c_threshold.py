from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import balanced_accuracy_score, f1_score

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.cross_platform import within_sample_percentile_rank
from src.data.provenance import expression_path, result_root
from src.modeling.evaluation import classification_metrics, prediction_probabilities
from src.modeling.pipelines import build_pipeline
from src.modeling.splitters import stratified_group_splits
from src.visualization.publication_figures import save_figure


def threshold_grid(probability: np.ndarray, labels: np.ndarray, config: dict) -> pd.DataFrame:
    rows = []
    thresholds = np.arange(
        float(config["threshold_grid_start"]),
        float(config["threshold_grid_stop"]) + float(config["threshold_grid_step"]) / 2,
        float(config["threshold_grid_step"]),
    )
    for threshold in thresholds:
        predicted = np.where(probability >= threshold, "tumor", "adjacent_normal")
        rows.append(
            {
                "threshold": float(np.round(threshold, 10)),
                "balanced_accuracy": balanced_accuracy_score(labels, predicted),
                "macro_f1": f1_score(labels, predicted, average="macro"),
                "distance_from_default": abs(threshold - float(config["default_threshold"])),
            }
        )
    frame = pd.DataFrame(rows)
    best = frame["balanced_accuracy"].max()
    candidates = frame[frame["balanced_accuracy"].eq(best)].copy()
    selected_index = candidates.sort_values(
        ["distance_from_default", "threshold"], ascending=[True, True]
    ).index[0]
    frame["selected"] = frame.index == selected_index
    frame["selection_rule"] = (
        "maximize balanced accuracy on patient-aggregated repeated grouped GSE44076 OOF "
        "rank-transport probabilities; ties closest to 0.5 then lower threshold"
    )
    return frame


def main() -> None:
    os.chdir(ROOT)
    config_text = (ROOT / "config/analysis.yaml").read_text(encoding="utf-8")
    config = yaml.safe_load(config_text)
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    provenance = config["modeling"]["expression_provenance"]
    tables = result_root(ROOT / paths["tables"], provenance)
    metrics_root = result_root(ROOT / paths["metrics"], provenance)
    figures = result_root(ROOT / paths["figures"], provenance)
    reports = result_root(ROOT / paths["reports"], provenance)
    models = result_root(ROOT / paths["models"], provenance)
    primary = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE44076", "gene", provenance)
    ).set_index("gene_symbol")
    external_gene_index = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE41258", "gene", provenance),
        columns=["gene_symbol"],
    )["gene_symbol"].astype(str)
    external_genes = set(external_gene_index)
    metadata = pd.read_csv(ROOT / paths["metadata"] / "gse44076_samples.csv")
    metadata = metadata[
        metadata["inclusion_status"].eq("included")
        & metadata["tissue_class"].isin(["adjacent_normal", "tumor"])
    ].reset_index(drop=True)
    card_path = models / "task_c_primary_full_signature_model_card.json"
    artifact_path = models / "task_c_primary_full_signature_model.joblib"
    card = json.loads(card_path.read_text(encoding="utf-8"))
    signature = list(card["feature_genes"])
    common_signature = [gene for gene in signature if gene in external_genes]
    common_universe = [gene for gene in primary.index.astype(str) if gene in external_genes]
    if len(common_signature) < 2:
        raise RuntimeError("Task C cross-platform signature contains fewer than two genes")

    samples = metadata["geo_accession"].to_numpy()
    labels = metadata["tissue_class"].to_numpy()
    groups = metadata["donor_or_patient_group"].astype(str).to_numpy()
    universe_values = primary.loc[common_universe, samples].T.to_numpy(dtype=np.float32)
    ranks = within_sample_percentile_rank(universe_values)
    positions = [common_universe.index(gene) for gene in common_signature]
    matrix = ranks[:, positions]
    modeling = config["modeling"]
    parameters = {
        **card["hyperparameters"],
        "variance_quantile__quantile": 0.0,
        "univariate__k": "all",
        "memory": None,
    }
    rows: list[dict[str, object]] = []
    for repeat, seed in enumerate(modeling["seeds"]):
        for fold, (train, test) in enumerate(
            stratified_group_splits(labels, groups, int(modeling["outer_splits"]), int(seed))
        ):
            model = build_pipeline("elastic_net", int(seed), modeling).set_params(**parameters)
            model.fit(matrix[train], labels[train])
            classes, probability = prediction_probabilities(model, matrix[test])
            tumor_index = int(np.flatnonzero(classes == "tumor")[0])
            for position, sample_index in enumerate(test):
                rows.append(
                    {
                        "repeat": repeat,
                        "outer_fold": fold,
                        "sample_id": samples[sample_index],
                        "patient_id": groups[sample_index],
                        "y_true": labels[sample_index],
                        "probability_tumor": float(probability[position, tumor_index]),
                    }
                )
    predictions = pd.DataFrame(rows)
    averaged = predictions.groupby("sample_id", as_index=False).agg(
        patient_id=("patient_id", "first"),
        y_true=("y_true", "first"),
        probability_tumor=("probability_tumor", "mean"),
        oof_repeats=("repeat", "nunique"),
    )
    if not averaged["oof_repeats"].eq(len(modeling["seeds"])).all():
        raise AssertionError("Each GSE44076 sample must have one OOF prediction per repeat")
    curve = threshold_grid(
        averaged["probability_tumor"].to_numpy(),
        averaged["y_true"].to_numpy(),
        config["external_validation"],
    )
    locked = float(curve.loc[curve["selected"], "threshold"].iloc[0])
    averaged["locked_threshold"] = locked
    averaged["y_pred"] = np.where(
        averaged["probability_tumor"].ge(locked), "tumor", "adjacent_normal"
    )
    averaged["probability_adjacent_normal"] = 1 - averaged["probability_tumor"]
    classes = np.asarray(["adjacent_normal", "tumor"])
    probability = averaged[["probability_adjacent_normal", "probability_tumor"]].to_numpy()
    locked_metrics = classification_metrics(
        averaged["y_true"].to_numpy(),
        averaged["y_pred"].to_numpy(),
        probability,
        classes,
    )
    predictions.to_csv(metrics_root / "task_c_rank_transport_repeated_oof.csv", index=False)
    averaged.to_csv(metrics_root / "task_c_rank_transport_aggregated_oof.csv", index=False)
    curve.to_csv(tables / "task_c_primary_threshold_selection.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.plot(curve["threshold"], curve["balanced_accuracy"], label="Balanced accuracy")
    ax.plot(curve["threshold"], curve["macro_f1"], label="Macro F1", alpha=0.8)
    ax.axvline(locked, color="black", linestyle="--", label=f"Locked: {locked:.2f}")
    ax.set(
        xlabel="Tumor probability threshold",
        ylabel="Aggregated grouped OOF metric",
        title="Task C threshold selected on GSE44076 only",
        ylim=(0.45, 1.02),
    )
    ax.legend()
    save_figure(fig, figures / "task_c_oof_threshold_curve")

    threshold_origin = (
        "GSE44076 patient-grouped repeated OOF rank-transport probabilities; "
        f"{config['external_validation']['threshold_rule']}; external labels not accessed"
    )
    card["transport_threshold"] = locked
    card["transport_threshold_origin"] = threshold_origin
    card["transport_signature_genes"] = common_signature
    card["transport_common_universe_gene_count"] = len(common_universe)
    card["external_labels_used"] = False
    card_path.write_text(json.dumps(card, indent=2), encoding="utf-8")
    artifact = joblib.load(artifact_path)
    artifact.update(card)
    joblib.dump(artifact, artifact_path)
    joblib.dump(artifact, models / "task_c_final_elastic_net.joblib")
    (models / "task_c_final_model_card.json").write_text(
        json.dumps(card, indent=2), encoding="utf-8"
    )

    report = f"""# Task C threshold selection

The cross-platform transport threshold was locked before external labels were
read. Repeated patient-grouped GSE44076 OOF probabilities were generated after
the label-independent GPL13667/GPL96 gene intersection and within-sample rank
transformation. Probabilities were averaged across repeats per sample. The
prespecified rule selected the balanced-accuracy maximum, resolving ties by
proximity to 0.5 and then the lower threshold.

- Locked threshold: {locked:.2f}
- Original Task C signature: {', '.join(signature)}
- Cross-platform signature: {', '.join(common_signature)}
- Label-independent common-gene universe: {len(common_universe)} genes
- Aggregated GSE44076 OOF ROC-AUC: {locked_metrics['roc_auc']:.4f}
- Aggregated GSE44076 OOF balanced accuracy: {locked_metrics['balanced_accuracy']:.4f}
- External labels used for threshold selection: no
"""
    (reports / "task_c_threshold_lock_report.md").write_text(report, encoding="utf-8")
    print(
        json.dumps(
            {
                "locked_threshold": locked,
                "common_signature": common_signature,
                "common_universe_gene_count": len(common_universe),
                "metrics": locked_metrics,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
