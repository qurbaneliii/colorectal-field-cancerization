from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from sklearn.calibration import calibration_curve
from sklearn.metrics import ConfusionMatrixDisplay, PrecisionRecallDisplay, RocCurveDisplay

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.cross_platform import within_sample_percentile_rank
from src.data.provenance import expression_path, result_root
from src.modeling.evaluation import grouped_bootstrap_metrics
from src.modeling.external_validation import locked_binary_predictions
from src.modeling.pipelines import build_pipeline
from src.visualization.publication_figures import save_figure


def deterministic_patient_subset(frame: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Select one array per patient without consulting tissue labels or outcomes."""

    selected = []
    for _, group in frame.groupby("donor_or_patient_group", sort=True):
        scores = group["geo_accession"].map(
            lambda sample: hashlib.sha256(f"{seed}|{sample}".encode()).hexdigest()
        )
        selected.append(group.loc[scores.idxmin()])
    return pd.DataFrame(selected).sort_values("geo_accession").reset_index(drop=True)


def fit_transport_model(
    artifact, primary, primary_meta, genes, common_universe, config, representation
):
    modeling = config["modeling"]
    model = build_pipeline("elastic_net", int(config["project"]["random_seed"]), modeling)
    parameters = {
        **artifact["hyperparameters"],
        "variance_quantile__quantile": 0.0,
        "univariate__k": "all",
        "memory": None,
    }
    model.set_params(**parameters)
    if representation == "within_sample_percentile_rank":
        all_ranks = within_sample_percentile_rank(
            primary.loc[common_universe, primary_meta["geo_accession"]].T.to_numpy(
                dtype=np.float32
            )
        )
        positions = [common_universe.index(gene) for gene in genes]
        matrix = all_ranks[:, positions]
    else:
        matrix = primary.loc[genes, primary_meta["geo_accession"]].T.to_numpy(dtype=np.float32)
    model.fit(matrix, primary_meta["tissue_class"].to_numpy())
    return model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--provenance", choices=["raw_cel_rma", "geo_deposited_series_matrix"])
    args = parser.parse_args()
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    provenance = args.provenance or config["modeling"]["expression_provenance"]
    metrics_root = result_root(ROOT / paths["metrics"], provenance)
    tables_root = result_root(ROOT / paths["tables"], provenance)
    figures_root = result_root(ROOT / paths["figures"], provenance)
    reports_root = result_root(ROOT / paths["reports"], provenance)
    models_root = result_root(ROOT / paths["models"], provenance)
    for directory in (metrics_root, tables_root, figures_root, reports_root, models_root):
        directory.mkdir(parents=True, exist_ok=True)
    model_path = models_root / "task_c_final_elastic_net.joblib"
    if not model_path.exists():
        raise FileNotFoundError("Run scripts/run_compact_panel_analysis.py first")
    artifact = joblib.load(model_path)
    if artifact["analysis_provenance"] != provenance:
        raise AssertionError("Model/expression provenance mismatch")

    primary = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE44076", "gene", provenance)
    ).set_index("gene_symbol")
    external = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE41258", "gene", provenance)
    ).set_index("gene_symbol")
    primary_meta = pd.read_csv(
        ROOT / paths["metadata"] / "gse44076_samples.csv", dtype={"patient_id": str}
    )
    external_meta = pd.read_csv(
        ROOT / paths["metadata"] / "gse41258_samples.csv", dtype={"patient_id": str}
    )
    included = external_meta[external_meta["inclusion_status"].eq("included")].copy()
    included["validation_label"] = included["tissue_class"].map(
        {"normal_colon": "adjacent_normal", "primary_tumor": "tumor"}
    )
    if included["validation_label"].isna().any():
        raise AssertionError("External eligible subset contains invalid tissue labels")
    if included.groupby(["patient_id", "tissue_class"]).size().gt(1).any():
        raise AssertionError("Unresolved multiple arrays per patient/tissue remain")
    canonical = deterministic_patient_subset(included, int(config["project"]["random_seed"]))
    if canonical["donor_or_patient_group"].duplicated().any():
        raise AssertionError("Canonical external subset is not one array per patient")

    signature = list(artifact["feature_genes"])
    common_signature = [gene for gene in signature if gene in primary.index and gene in external.index]
    external_gene_set = set(external.index)
    common_universe = [gene for gene in primary.index if gene in external_gene_set]
    excluded_signature = [gene for gene in signature if gene not in common_signature]
    if len(common_signature) < 2:
        raise RuntimeError("Fewer than two locked signature genes are common to both platforms")
    primary_task = primary_meta[
        primary_meta["tissue_class"].isin(["adjacent_normal", "tumor"])
        & primary_meta["inclusion_status"].eq("included")
    ]
    threshold = float(artifact["threshold"])
    primary_representation = config["external_validation"]["primary_representation"]
    representations = [
        primary_representation,
        *config["external_validation"].get("sensitivity_representations", []),
    ]
    if representations[0] != "within_sample_percentile_rank":
        raise AssertionError("Primary external representation must remain prespecified as ranks")

    prediction_frames = []
    metric_rows = []
    fitted_models = {}
    for representation in representations:
        model = fit_transport_model(
            artifact,
            primary,
            primary_task,
            common_signature,
            common_universe,
            config,
            representation,
        )
        fitted_models[representation] = model
        for evaluation_set, frame in (("canonical_patient", canonical), ("all_array_sensitivity", included)):
            matrix = external.loc[common_signature, frame["geo_accession"]].T.to_numpy(
                dtype=np.float32
            )
            if representation == "within_sample_percentile_rank":
                all_ranks = within_sample_percentile_rank(
                    external.loc[common_universe, frame["geo_accession"]].T.to_numpy(
                        dtype=np.float32
                    )
                )
                positions = [common_universe.index(gene) for gene in common_signature]
                matrix = all_ranks[:, positions]
            prediction, measured = locked_binary_predictions(
                model,
                matrix,
                frame["geo_accession"].to_numpy(),
                frame["donor_or_patient_group"].astype(str).to_numpy(),
                frame["validation_label"].to_numpy(),
                representation,
                threshold,
            )
            prediction["evaluation_set"] = evaluation_set
            prediction["analysis_provenance"] = provenance
            prediction_frames.append(prediction)
            metric_rows.append(
                {
                    "representation": representation,
                    "evaluation_set": evaluation_set,
                    "is_primary": representation == primary_representation
                    and evaluation_set == "canonical_patient",
                    "arrays": len(frame),
                    "unique_patients": frame["donor_or_patient_group"].nunique(),
                    "analysis_provenance": provenance,
                    **measured,
                }
            )
    predictions = pd.concat(prediction_frames, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    ci_rows = []
    for (representation, evaluation_set), frame in predictions.groupby(
        ["representation", "evaluation_set"]
    ):
        summaries = grouped_bootstrap_metrics(
            frame,
            [
                "roc_auc",
                "pr_auc",
                "balanced_accuracy",
                "sensitivity",
                "specificity",
                "f1_macro",
                "f1_weighted",
                "brier_score",
                "log_loss",
                "negative_predictive_value",
                "positive_predictive_value",
                "calibration_intercept",
                "calibration_slope",
            ],
            int(config["external_validation"]["bootstrap_iterations"]),
            int(config["project"]["random_seed"]),
        )
        ci_rows.extend(
            {
                "representation": representation,
                "evaluation_set": evaluation_set,
                "analysis_provenance": provenance,
                **summary,
            }
            for summary in summaries
        )
    ci = pd.DataFrame(ci_rows)
    predictions.to_csv(metrics_root / "external_validation_predictions.csv", index=False)
    metrics.to_csv(metrics_root / "external_validation_metrics.csv", index=False)
    ci.to_csv(metrics_root / "external_validation_bootstrap_ci.csv", index=False)
    metrics.to_csv(tables_root / "table_6_external_validation.csv", index=False)
    predictions.to_csv(tables_root / "supplementary_external_predictions.csv", index=False)

    common_table = pd.DataFrame(
        {
            "gene_symbol": signature,
            "present_gse44076": [gene in primary.index for gene in signature],
            "present_gse41258": [gene in external.index for gene in signature],
            "included_external_signature": [gene in common_signature for gene in signature],
            "exclusion_reason": [
                "" if gene in common_signature else "not in label-independent cross-platform intersection"
                for gene in signature
            ],
        }
    )
    common_table.to_csv(tables_root / "external_common_signature_audit.csv", index=False)
    patient_counts = included.groupby("patient_id").agg(
        arrays=("geo_accession", "size"),
        tissues=("tissue_class", lambda values: ";".join(sorted(set(values)))),
    )
    structure = pd.DataFrame(
        [
            ["eligible_arrays", len(included)],
            ["eligible_unique_patients", included["patient_id"].nunique()],
            ["patients_with_both_tissues", patient_counts["tissues"].str.contains(";").sum()],
            ["tumor_only_patients", patient_counts["tissues"].eq("primary_tumor").sum()],
            ["normal_only_patients", patient_counts["tissues"].eq("normal_colon").sum()],
            ["technical_replicates_removed", external_meta["technical_replicate_candidate"].sum()],
            ["canonical_primary_arrays", len(canonical)],
        ],
        columns=["quantity", "value"],
    )
    structure.to_csv(tables_root / "external_patient_structure.csv", index=False)

    rank_model = fitted_models["within_sample_percentile_rank"]
    joblib.dump(
        {
            "model": rank_model,
            "signature_genes": common_signature,
            "common_gene_universe": common_universe,
            "threshold": threshold,
            "training_accession": "GSE44076",
            "analysis_provenance": provenance,
            "representation": "within_sample_percentile_rank",
            "external_labels_used_for_training_or_selection": False,
        },
        models_root / "task_c_final_rank_elastic_net.joblib",
    )

    primary_predictions = predictions[
        predictions["representation"].eq(primary_representation)
        & predictions["evaluation_set"].eq("canonical_patient")
    ]
    positive = primary_predictions["y_true"].eq("tumor").astype(int)
    probability = primary_predictions["probability_tumor"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    RocCurveDisplay.from_predictions(positive, probability, ax=axes[0, 0])
    PrecisionRecallDisplay.from_predictions(positive, probability, ax=axes[0, 1])
    observed, predicted = calibration_curve(positive, probability, n_bins=8, strategy="quantile")
    axes[1, 0].plot(predicted, observed, marker="o")
    axes[1, 0].plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.8)
    axes[1, 0].set(
        xlabel="Mean predicted tumor probability",
        ylabel="Observed tumor fraction",
        title="External calibration",
    )
    ConfusionMatrixDisplay.from_predictions(
        primary_predictions["y_true"],
        primary_predictions["y_pred"],
        labels=["adjacent_normal", "tumor"],
        display_labels=["Normal colon", "Primary tumor"],
        colorbar=False,
        ax=axes[1, 1],
    )
    axes[1, 1].set_title("Locked threshold confusion matrix")
    fig.suptitle(
        f"GSE41258 patient-balanced external validation (n={len(canonical)} patients; "
        f"{len(common_signature)} genes)"
    )
    save_figure(fig, figures_root / "external_validation_primary")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    RocCurveDisplay.from_predictions(positive, probability, ax=axes[0])
    PrecisionRecallDisplay.from_predictions(positive, probability, ax=axes[1])
    fig.suptitle(f"External discrimination: {primary_representation} (n={len(canonical)})")
    save_figure(fig, figures_root / "external_validation_roc_pr")
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.plot(predicted, observed, marker="o")
    ax.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.8)
    ax.set(
        xlabel="Mean predicted tumor probability",
        ylabel="Observed tumor fraction",
        title=f"External calibration (n={len(canonical)})",
    )
    save_figure(fig, figures_root / "external_calibration")
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ConfusionMatrixDisplay.from_predictions(
        primary_predictions["y_true"],
        primary_predictions["y_pred"],
        labels=["adjacent_normal", "tumor"],
        display_labels=["Normal colon", "Primary tumor"],
        colorbar=False,
        ax=ax,
    )
    ax.set_title(f"External confusion matrix (n={len(canonical)})")
    save_figure(fig, figures_root / "external_confusion_matrix")

    primary_metric = metrics[metrics["is_primary"]].iloc[0]
    report = f"""# External validation report

The prespecified primary transfer representation was within-sample percentile
rank because it is computed independently within each sample, does not require
external labels, reduces reliance on platform-specific intensity scale, and
preserves relative gene ordering. Training-derived z-scores are a sensitivity
analysis. GSE41258 labels were not used for feature, panel, model,
transformation, threshold, or hyperparameter selection.

- Eligible arrays: {len(included)} ({int(included['tissue_class'].eq('primary_tumor').sum())} primary tumor; {int(included['tissue_class'].eq('normal_colon').sum())} normal colon)
- Eligible unique patients: {included['patient_id'].nunique()}
- Primary deterministic one-array-per-patient subset: {len(canonical)}
- Original Task C signature: {len(signature)} genes
- Cross-platform common signature: {len(common_signature)} genes
- Excluded signature genes: {', '.join(excluded_signature) if excluded_signature else 'none'}
- Locked threshold: {threshold}

Primary point estimates: ROC-AUC {primary_metric['roc_auc']:.4f}, PR-AUC
{primary_metric['pr_auc']:.4f}, balanced accuracy
{primary_metric['balanced_accuracy']:.4f}, macro F1
{primary_metric['f1_macro']:.4f}, Brier score
{primary_metric['brier_score']:.4f}, and log loss
{primary_metric['log_loss']:.4f}. Patient-cluster bootstrap intervals and the
all-array sensitivity analysis are saved in `results/metrics`.

This validates only the tumor-versus-normal-colon component. GSE41258 is not an
independent validation cohort for cancer-free healthy versus tumor-adjacent
field cancerization.
"""
    (reports_root / "external_validation_report.md").write_text(report, encoding="utf-8")
    print(
        json.dumps(
            {
                "provenance": provenance,
                "eligible_arrays": len(included),
                "unique_patients": included["patient_id"].nunique(),
                "canonical_primary_arrays": len(canonical),
                "signature_genes_primary": len(signature),
                "signature_genes_common": len(common_signature),
                "primary_metrics": primary_metric.to_dict(),
            },
            indent=2,
            default=str,
        )
    )


if __name__ == "__main__":
    main()
