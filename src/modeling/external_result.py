from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.modeling.evaluation import classification_metrics


EVALUATION_KEYS = ["representation", "evaluation_set", "threshold_policy"]
PRIMARY_KEY = (
    "within_sample_percentile_rank",
    "canonical_patient_tissue",
    "primary_gse44076_locked",
)
LABELS = np.array(["normal_colon", "primary_tumor"])


def external_summary(
    metrics: pd.DataFrame,
    predictions: pd.DataFrame,
    threshold_grid: pd.DataFrame,
    full_card: dict,
    transport_card: dict,
    signature_audit: pd.DataFrame,
) -> dict:
    """Cross-check the locked evaluation against its component artifacts."""
    if metrics.duplicated(EVALUATION_KEYS).any() or metrics["is_primary"].sum() != 1:
        raise ValueError("External metrics require unique evaluation keys and one primary row")
    primary = metrics.loc[metrics["is_primary"]].iloc[0]
    if tuple(primary[EVALUATION_KEYS]) != PRIMARY_KEY:
        raise ValueError("Unexpected primary external evaluation")
    for source in (predictions, metrics):
        if set(source["analysis_provenance"]) != {transport_card["analysis_provenance"]}:
            raise ValueError("External expression provenance differs across artifacts")
    if full_card["analysis_provenance"] != transport_card["analysis_provenance"]:
        raise ValueError("Task C model cards use different expression provenance")
    full_genes = full_card["feature_genes"]
    transport_genes = transport_card["feature_genes"]
    if (
        signature_audit["gene_symbol"].tolist() != full_genes
        or signature_audit.loc[
            signature_audit["included_transport_signature"], "gene_symbol"
        ].tolist() != transport_genes
        or full_card["transport_signature_genes"] != transport_genes
        or [gene for gene in full_genes if gene not in transport_genes]
        != transport_card["excluded_genes"]
    ):
        raise ValueError("Task C signatures differ across the model cards and mapping audit")
    selected = threshold_grid.loc[threshold_grid["selected"], "threshold"].tolist()
    threshold = float(transport_card["threshold"])
    if (
        selected != [threshold]
        or full_card["transport_threshold"] != threshold
        or full_card["transport_threshold_origin"] != transport_card["threshold_origin"]
        or full_card["external_labels_used"] is not False
        or transport_card["external_labels_used"] is not False
    ):
        raise ValueError("Task C development threshold provenance differs across artifacts")

    grouped = predictions.groupby(EVALUATION_KEYS, sort=False)
    if set(grouped.groups) != set(metrics[EVALUATION_KEYS].itertuples(index=False, name=None)):
        raise ValueError("External metrics and predictions cover different evaluations")
    for row in metrics.to_dict(orient="records"):
        key = tuple(row[column] for column in EVALUATION_KEYS)
        frame = grouped.get_group(key)
        probabilities = frame[[f"probability_{label}" for label in LABELS]].to_numpy()
        limits = frame["locked_threshold"].unique()
        if (
            frame["sample_id"].duplicated().any()
            or frame["patient_id"].isna().any()
            or not np.isfinite(probabilities).all()
            or ((probabilities < 0) | (probabilities > 1)).any()
            or not np.allclose(probabilities.sum(axis=1), 1)
            or len(limits) != 1
            or not 0 <= limits[0] <= 1
            or (key[2] == PRIMARY_KEY[2] and limits[0] != threshold)
        ):
            raise ValueError(f"Invalid external predictions for {key}")
        observed = frame["y_true"].to_numpy()
        predicted = frame["y_pred"].to_numpy()
        if set(observed) != set(LABELS) or not np.array_equal(
            predicted, np.where(probabilities[:, 1] >= limits[0], LABELS[1], LABELS[0])
        ):
            raise ValueError(f"External class or threshold mapping disagrees for {key}")
        if (
            row["arrays"] != len(frame)
            or row["unique_patients"] != frame["patient_id"].nunique()
            or row["normal_colon_arrays"] != int((observed == LABELS[0]).sum())
            or row["primary_tumor_arrays"] != int((observed == LABELS[1]).sum())
            or (key[1] == PRIMARY_KEY[1] and frame.duplicated(["patient_id", "y_true"]).any())
        ):
            raise ValueError(f"External sample or patient counts differ for {key}")
        measured = classification_metrics(
            observed, predicted, probabilities, LABELS, include_calibration=True
        )
        for metric, value in measured.items():
            tolerance = 1e-4 if metric.startswith("calibration_") else 1e-8
            if not np.isclose(row[metric], value, rtol=tolerance, atol=1e-10, equal_nan=True):
                raise ValueError(f"External metric {metric} disagrees with predictions for {key}")

    default_key = (*PRIMARY_KEY[:2], "default_0_5")
    default_threshold = float(grouped.get_group(default_key)["locked_threshold"].iloc[0])
    return {
        "schema_version": 1,
        "training_accession": "GSE44076",
        "evaluation_accession": "GSE41258",
        "analysis_provenance": transport_card["analysis_provenance"],
        "endpoint": "Task C, primary tumor versus normal colon in the neoplasia cohort",
        "full_signature_genes": full_genes,
        "transport_signature_genes": transport_genes,
        "excluded_signature_genes": transport_card["excluded_genes"],
        "common_universe_gene_count": transport_card["common_gene_universe_count"],
        "locked_threshold": threshold,
        "threshold_origin": transport_card["threshold_origin"],
        "threshold_selection_scope": full_card.get(
            "transport_threshold_selection_scope", "fixed_full_development_panel_and_hyperparameters"
        ),
        "default_threshold_sensitivity": default_threshold,
        "primary_evaluation": primary.astype(object).where(primary.notna(), None).to_dict(),
    }


def external_report(summary: dict) -> str:
    primary = summary["primary_evaluation"]
    return f"""# External validation report

The GSE41258 endpoint was primary tumor versus normal colon in patients with
colorectal neoplasia. Task B requires a separate cohort with cancer-free controls
and cancer-associated normal mucosa collected under comparable protocols.

The full GSE44076 signature contained {', '.join(summary['full_signature_genes'])}.
The cross-platform model was refitted on GSE44076 with
{', '.join(summary['transport_signature_genes'])}. The threshold was selected on
grouped GSE44076 out-of-fold scores with the full-development panel and
hyperparameters held fixed. GSE41258 outcomes were used for evaluation.
The cohort has been examined in earlier project analyses.

- Locked threshold {summary['locked_threshold']:.2f}
- Default threshold sensitivity {summary['default_threshold_sensitivity']:.2f}
- Canonical set of {primary['arrays']} arrays from {primary['unique_patients']} patients
- ROC-AUC {primary['roc_auc']:.4f} and average precision {primary['pr_auc']:.4f}
- Balanced accuracy {primary['balanced_accuracy']:.4f} and macro-F1 {primary['f1_macro']:.4f}
- Sensitivity {primary['sensitivity']:.4f} and specificity {primary['specificity']:.4f}
- Brier score {primary['brier_score']:.4f} and log loss {primary['log_loss']:.4f}
- Calibration intercept {primary['calibration_intercept']:.4f} and slope {primary['calibration_slope']:.4f}

Patient-cluster bootstrap intervals and the sensitivity estimands are in the
external validation CSVs. Calibration coefficients describe score transport in
this cohort. They were not applied to the model or the locked predictions.
"""


def write_external_report(reports: Path, summary: dict) -> None:
    (reports / "external_validation_report.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    (reports / "external_validation_report.md").write_text(
        external_report(summary), encoding="utf-8"
    )
