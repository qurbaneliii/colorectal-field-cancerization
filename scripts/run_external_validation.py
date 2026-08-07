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
from src.data.provenance import expression_path, result_root, source_files_sha256
from src.modeling.evaluation import grouped_bootstrap_metrics
from src.modeling.external_validation import locked_binary_predictions
from src.modeling.pipelines import build_pipeline
from src.visualization.publication_figures import save_figure


TRANSPORT_SOURCE_FILES = [
    "config/analysis.yaml",
    "scripts/run_external_validation.py",
    "src/data/cross_platform.py",
    "src/modeling/evaluation.py",
    "src/modeling/external_validation.py",
    "src/modeling/pipelines.py",
]


def deterministic_subset(
    frame: pd.DataFrame, grouping: list[str], seed: int
) -> pd.DataFrame:
    """Select a canonical array per group without consulting outcome performance."""

    selected = []
    for _, group in frame.groupby(grouping, sort=True):
        scores = group["geo_accession"].map(
            lambda sample: hashlib.sha256(f"{seed}|{sample}".encode()).hexdigest()
        )
        selected.append(group.loc[scores.idxmin()])
    return pd.DataFrame(selected).sort_values("geo_accession").reset_index(drop=True)


def fit_transport_model(
    artifact: dict,
    primary: pd.DataFrame,
    primary_meta: pd.DataFrame,
    genes: list[str],
    common_universe: list[str],
    config: dict,
    representation: str,
):
    modeling = config["modeling"]
    model = build_pipeline("elastic_net", int(config["project"]["random_seed"]), modeling)
    model.set_params(
        **{
            **artifact["hyperparameters"],
            "variance_quantile__quantile": 0.0,
            "univariate__k": "all",
            "memory": None,
        }
    )
    if representation == "within_sample_percentile_rank":
        all_ranks = within_sample_percentile_rank(
            primary.loc[common_universe, primary_meta["geo_accession"]].T.to_numpy(
                dtype=np.float32
            )
        )
        positions = [common_universe.index(gene) for gene in genes]
        matrix = all_ranks[:, positions]
    elif representation == "training_zscore":
        matrix = primary.loc[genes, primary_meta["geo_accession"]].T.to_numpy(
            dtype=np.float32
        )
    else:
        raise ValueError(f"Unsupported representation: {representation}")
    model.fit(matrix, primary_meta["tissue_class"].to_numpy())
    return model


def transport_matrix(
    expression: pd.DataFrame,
    frame: pd.DataFrame,
    genes: list[str],
    common_universe: list[str],
    representation: str,
) -> np.ndarray:
    if representation == "within_sample_percentile_rank":
        ranks = within_sample_percentile_rank(
            expression.loc[common_universe, frame["geo_accession"]].T.to_numpy(
                dtype=np.float32
            )
        )
        positions = [common_universe.index(gene) for gene in genes]
        return ranks[:, positions]
    return expression.loc[genes, frame["geo_accession"]].T.to_numpy(dtype=np.float32)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--provenance", choices=["raw_cel_rma", "geo_deposited_series_matrix"]
    )
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

    primary_model_path = models_root / "task_c_primary_full_signature_model.joblib"
    if not primary_model_path.exists():
        raise FileNotFoundError("Run compact-panel analysis and threshold locking first")
    artifact = joblib.load(primary_model_path)
    if artifact["analysis_provenance"] != provenance:
        raise AssertionError("Model/expression provenance mismatch")
    if "transport_threshold" not in artifact:
        raise AssertionError("Task C transport threshold has not been locked on GSE44076")

    primary = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE44076", "gene", provenance)
    ).set_index("gene_symbol")
    external = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE41258", "gene", provenance)
    ).set_index("gene_symbol")
    primary_meta = pd.read_csv(ROOT / paths["metadata"] / "gse44076_samples.csv")
    external_meta = pd.read_csv(
        ROOT / paths["metadata"] / "gse41258_samples.csv", dtype={"patient_id": str}
    )
    eligible = external_meta[
        external_meta["inclusion_status"].eq("included")
        & external_meta["tissue_class"].isin(["normal_colon", "primary_tumor"])
    ].copy()
    if eligible["technical_replicate_candidate"].any():
        raise AssertionError("Eligible external set contains technical replicate candidates")
    canonical = deterministic_subset(
        eligible,
        ["donor_or_patient_group", "tissue_class"],
        int(config["project"]["random_seed"]),
    )
    if canonical.groupby(["donor_or_patient_group", "tissue_class"]).size().gt(1).any():
        raise AssertionError("Canonical set has multiple arrays per patient and tissue")
    one_per_patient = deterministic_subset(
        eligible, ["donor_or_patient_group"], int(config["project"]["random_seed"])
    )

    signature = list(artifact["feature_genes"])
    external_genes = set(external.index.astype(str))
    common_signature = [
        gene for gene in signature if gene in primary.index and gene in external_genes
    ]
    common_universe = [gene for gene in primary.index.astype(str) if gene in external_genes]
    excluded_signature = [gene for gene in signature if gene not in common_signature]
    if len(common_signature) < 2:
        raise RuntimeError("Fewer than two locked genes occur on both platforms")
    primary_task = primary_meta[
        primary_meta["inclusion_status"].eq("included")
        & primary_meta["tissue_class"].isin(["adjacent_normal", "tumor"])
    ].copy()
    locked_threshold = float(artifact["transport_threshold"])
    default_threshold = float(config["external_validation"]["default_threshold"])
    primary_representation = config["external_validation"]["primary_representation"]
    representations = [
        primary_representation,
        *config["external_validation"].get("sensitivity_representations", []),
    ]
    if primary_representation != "within_sample_percentile_rank":
        raise AssertionError("Primary external representation must be within-sample rank")
    evaluation_sets = {
        "canonical_patient_tissue": canonical,
        "all_array_patient_cluster_sensitivity": eligible,
        "one_array_per_patient_sensitivity": one_per_patient,
    }
    threshold_policies = {
        "primary_gse44076_locked": locked_threshold,
        "default_0_5": default_threshold,
    }

    prediction_frames: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []
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
        for evaluation_set, frame in evaluation_sets.items():
            matrix = transport_matrix(
                external, frame, common_signature, common_universe, representation
            )
            for threshold_policy, threshold in threshold_policies.items():
                prediction, measured = locked_binary_predictions(
                    model,
                    matrix,
                    frame["geo_accession"].to_numpy(),
                    frame["donor_or_patient_group"].astype(str).to_numpy(),
                    frame["tissue_class"].to_numpy(),
                    representation,
                    threshold,
                )
                prediction["evaluation_set"] = evaluation_set
                prediction["threshold_policy"] = threshold_policy
                prediction["analysis_provenance"] = provenance
                prediction_frames.append(prediction)
                metric_rows.append(
                    {
                        "representation": representation,
                        "evaluation_set": evaluation_set,
                        "threshold_policy": threshold_policy,
                        "is_primary": (
                            representation == primary_representation
                            and evaluation_set == "canonical_patient_tissue"
                            and threshold_policy == "primary_gse44076_locked"
                        ),
                        "arrays": len(frame),
                        "unique_patients": frame["donor_or_patient_group"].nunique(),
                        "normal_colon_arrays": int(frame["tissue_class"].eq("normal_colon").sum()),
                        "primary_tumor_arrays": int(
                            frame["tissue_class"].eq("primary_tumor").sum()
                        ),
                        "analysis_provenance": provenance,
                        **measured,
                    }
                )
    predictions = pd.concat(prediction_frames, ignore_index=True)
    metrics = pd.DataFrame(metric_rows)
    ci_rows: list[dict[str, object]] = []
    for keys, frame in predictions.groupby(
        ["representation", "evaluation_set", "threshold_policy"]
    ):
        representation, evaluation_set, threshold_policy = keys
        bootstrap_metrics = [
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
        ]
        if (
            representation == primary_representation
            and evaluation_set == "canonical_patient_tissue"
            and threshold_policy == "primary_gse44076_locked"
        ):
            bootstrap_metrics.extend(["calibration_intercept", "calibration_slope"])
        summaries = grouped_bootstrap_metrics(
            frame,
            bootstrap_metrics,
            int(config["external_validation"]["bootstrap_iterations"]),
            int(config["project"]["random_seed"]),
        )
        ci_rows.extend(
            {
                "representation": representation,
                "evaluation_set": evaluation_set,
                "threshold_policy": threshold_policy,
                "cluster_unit": "patient",
                "analysis_provenance": provenance,
                **summary,
            }
            for summary in summaries
        )
    intervals = pd.DataFrame(ci_rows)
    predictions.to_csv(metrics_root / "external_validation_predictions.csv", index=False)
    metrics.to_csv(metrics_root / "external_patient_tissue_metrics.csv", index=False)
    intervals.to_csv(metrics_root / "external_patient_cluster_bootstrap_ci.csv", index=False)
    metrics.to_csv(metrics_root / "external_validation_metrics.csv", index=False)
    intervals.to_csv(metrics_root / "external_validation_bootstrap_ci.csv", index=False)
    metrics.to_csv(tables_root / "table_6_external_validation.csv", index=False)
    predictions.to_csv(tables_root / "supplementary_external_predictions.csv", index=False)

    common_table = pd.DataFrame(
        {
            "gene_symbol": signature,
            "present_gse44076": [gene in primary.index for gene in signature],
            "present_gse41258": [gene in external.index for gene in signature],
            "included_transport_signature": [gene in common_signature for gene in signature],
            "exclusion_reason": [
                ""
                if gene in common_signature
                else "absent from label-independent GPL13667/GPL96 intersection"
                for gene in signature
            ],
        }
    )
    common_table.to_csv(tables_root / "external_common_signature_audit.csv", index=False)
    patient_structure = eligible.groupby("donor_or_patient_group", as_index=False).agg(
        arrays=("geo_accession", "size"),
        normal_colon_arrays=("tissue_class", lambda x: int((x == "normal_colon").sum())),
        primary_tumor_arrays=(
            "tissue_class",
            lambda x: int((x == "primary_tumor").sum()),
        ),
    )
    patient_structure["has_both_tissues"] = (
        patient_structure["normal_colon_arrays"].gt(0)
        & patient_structure["primary_tumor_arrays"].gt(0)
    )
    patient_structure["canonical_arrays"] = patient_structure[
        ["normal_colon_arrays", "primary_tumor_arrays"]
    ].gt(0).sum(axis=1)
    patient_structure.to_csv(
        tables_root / "external_patient_tissue_structure.csv", index=False
    )

    transport_card = {
        "model_version": "2.0.0-transport",
        "task": "task_c_tumor_vs_normal_colon_cross_platform_transport",
        "training_accession": "GSE44076",
        "training_representation": "within_sample_percentile_rank",
        "feature_genes": common_signature,
        "excluded_genes": excluded_signature,
        "exclusion_reason": "absent from label-independent cross-platform common-gene set",
        "common_gene_universe": common_universe,
        "common_gene_universe_count": len(common_universe),
        "refitting_occurred": True,
        "refitting_description": (
            "Locked feature specification, family, hyperparameters, representation, and threshold; "
            "refitted on GSE44076 common signature before GSE41258 evaluation"
        ),
        "external_labels_used_for_training_selection_or_threshold": False,
        "external_labels_used": False,
        "threshold": locked_threshold,
        "threshold_origin": artifact["transport_threshold_origin"],
        "analysis_provenance": provenance,
        "analysis_source_files": TRANSPORT_SOURCE_FILES,
        "analysis_source_sha256": source_files_sha256(ROOT, TRANSPORT_SOURCE_FILES),
        "clinical_readiness": "not clinically ready",
    }
    joblib.dump(
        {"model": fitted_models[primary_representation], **transport_card},
        models_root / "task_c_cross_platform_transport_model.joblib",
    )
    (models_root / "task_c_cross_platform_transport_model_card.json").write_text(
        json.dumps(transport_card, indent=2), encoding="utf-8"
    )

    primary_predictions = predictions[
        predictions["representation"].eq(primary_representation)
        & predictions["evaluation_set"].eq("canonical_patient_tissue")
        & predictions["threshold_policy"].eq("primary_gse44076_locked")
    ].copy()
    positive = primary_predictions["y_true"].eq("primary_tumor").astype(int)
    probability = primary_predictions["probability_primary_tumor"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    RocCurveDisplay.from_predictions(positive, probability, ax=axes[0, 0])
    PrecisionRecallDisplay.from_predictions(positive, probability, ax=axes[0, 1])
    observed, predicted_probability = calibration_curve(
        positive, probability, n_bins=8, strategy="quantile"
    )
    axes[1, 0].plot(predicted_probability, observed, marker="o")
    axes[1, 0].plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.8)
    axes[1, 0].set(
        xlabel="Mean predicted primary-tumor probability",
        ylabel="Observed primary-tumor fraction",
        title="External calibration",
    )
    ConfusionMatrixDisplay.from_predictions(
        primary_predictions["y_true"],
        primary_predictions["y_pred"],
        labels=["normal_colon", "primary_tumor"],
        display_labels=["Normal colon", "Primary tumor"],
        colorbar=False,
        ax=axes[1, 1],
    )
    axes[1, 1].set_title("GSE44076-locked threshold")
    fig.suptitle(
        f"GSE41258 canonical patient-tissue validation "
        f"(n={len(canonical)} arrays; {canonical['donor_or_patient_group'].nunique()} patients)"
    )
    save_figure(fig, figures_root / "external_validation_primary")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    RocCurveDisplay.from_predictions(positive, probability, ax=axes[0])
    PrecisionRecallDisplay.from_predictions(positive, probability, ax=axes[1])
    fig.suptitle("GSE41258 discrimination: canonical patient-tissue set")
    save_figure(fig, figures_root / "external_validation_roc_pr")
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.plot(predicted_probability, observed, marker="o")
    ax.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.8)
    ax.set(
        xlabel="Mean predicted primary-tumor probability",
        ylabel="Observed primary-tumor fraction",
        title="External calibration",
    )
    save_figure(fig, figures_root / "external_calibration")
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ConfusionMatrixDisplay.from_predictions(
        primary_predictions["y_true"],
        primary_predictions["y_pred"],
        labels=["normal_colon", "primary_tumor"],
        display_labels=["Normal colon", "Primary tumor"],
        colorbar=False,
        ax=ax,
    )
    ax.set_title("External confusion matrix: locked threshold")
    save_figure(fig, figures_root / "external_confusion_matrix")

    counts = patient_structure
    both = int(counts["has_both_tissues"].sum())
    tumor_only = int(
        (counts["primary_tumor_arrays"].gt(0) & counts["normal_colon_arrays"].eq(0)).sum()
    )
    normal_only = int(
        (counts["normal_colon_arrays"].gt(0) & counts["primary_tumor_arrays"].eq(0)).sum()
    )
    technical_removed = int(
        external_meta[
            external_meta["tissue_class"].isin(["normal_colon", "primary_tumor"])
        ]["technical_replicate_candidate"].sum()
    )
    primary_metric = metrics[metrics["is_primary"]].iloc[0]
    structure_report = f"""# External patient structure report

The primary GSE41258 evaluation retains one canonical array per patient and
tissue class. A patient with both normal-colon and primary-tumor tissue
therefore contributes both observations. Uncertainty resamples patients and
retains every selected observation for each sampled patient.

- Eligible arrays: {len(eligible)}
- Canonical patient-tissue arrays: {len(canonical)}
- Unique patients: {eligible['donor_or_patient_group'].nunique()}
- Patients with both tissues: {both}
- Tumor-only patients: {tumor_only}
- Normal-only patients: {normal_only}
- Technical replicate candidates removed: {technical_removed}
- Samples per patient: median {counts['arrays'].median():.0f}, range {counts['arrays'].min()}-{counts['arrays'].max()}
- Optional one-array-per-patient sensitivity: {len(one_per_patient)} arrays
"""
    (reports_root / "external_patient_structure_report.md").write_text(
        structure_report, encoding="utf-8"
    )
    validation_report = f"""# External validation report

This is a secondary cross-platform tumor-normal evaluation, not independent
validation of the healthy-versus-adjacent field-effect task. The exact
{len(signature)}-gene Task C model could not be deployed because
{', '.join(excluded_signature)} {'was' if len(excluded_signature) == 1 else 'were'}
absent on GPL96. The feature specification, Elastic Net family,
hyperparameters, within-sample rank representation, and threshold were locked
using GSE44076. A separate transport model was refitted on GSE44076 using the
{len(common_signature)} common signature genes and evaluated once on GSE41258.
External labels were never used for fitting, selection, transformation, or
threshold choice.

- External labels retained verbatim: `normal_colon`, `primary_tumor`
- Locked threshold: {locked_threshold:.2f} (GSE44076 grouped repeated OOF)
- Default threshold sensitivity: {default_threshold:.2f}
- Primary set: {len(canonical)} arrays from {canonical['donor_or_patient_group'].nunique()} patients
- Primary ROC-AUC: {primary_metric['roc_auc']:.4f}
- Primary PR-AUC: {primary_metric['pr_auc']:.4f}
- Primary balanced accuracy: {primary_metric['balanced_accuracy']:.4f}
- Primary macro-F1: {primary_metric['f1_macro']:.4f}
- Primary NPV: {primary_metric['negative_predictive_value']:.4f}
- Primary Brier score: {primary_metric['brier_score']:.4f}
- Primary log loss: {primary_metric['log_loss']:.4f}

Patient-cluster bootstrap intervals, the all-array analysis, the optional
one-array-per-patient analysis, and the failed or degraded training-z-score
transport are reported as distinct estimands. High ranking discrimination does
not by itself establish threshold or probability calibration transport. This
retrospective biomarker-discovery model is not clinically ready.
"""
    (reports_root / "external_validation_report.md").write_text(
        validation_report, encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "eligible_arrays": len(eligible),
                "canonical_arrays": len(canonical),
                "unique_patients": eligible["donor_or_patient_group"].nunique(),
                "patients_with_both_tissues": both,
                "signature_genes_primary": signature,
                "signature_genes_transport": common_signature,
                "locked_threshold": locked_threshold,
                "primary_metrics": primary_metric.to_dict(),
            },
            indent=2,
            default=str,
        )
    )


if __name__ == "__main__":
    main()
