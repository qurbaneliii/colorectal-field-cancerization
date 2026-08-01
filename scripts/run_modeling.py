from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.provenance import expression_path, result_root
from src.data.validation import validate_fold_assignments
from src.interpretation.coefficients import coefficient_plot
from src.modeling.evaluation import classification_metrics, grouped_bootstrap_metrics
from src.modeling.nested_cv import run_nested_cv
from src.modeling.stability import feature_stability


TASKS = {
    "task_a_three_class": ["healthy", "adjacent_normal", "tumor"],
    "task_b_field_effect": ["healthy", "adjacent_normal"],
    "task_c_tumor_vs_adjacent": ["adjacent_normal", "tumor"],
}


def average_repeated_predictions(frame: pd.DataFrame) -> pd.DataFrame:
    probability_columns = [
        column for column in frame if column.startswith("probability_") and frame[column].notna().any()
    ]
    averaged = frame.groupby("sample_id", as_index=False).agg(
        {
            "patient_id": "first",
            "y_true": "first",
            **{column: "mean" for column in probability_columns},
        }
    )
    classes = [column.removeprefix("probability_") for column in probability_columns]
    averaged["y_pred"] = np.asarray(classes)[
        np.argmax(averaged[probability_columns].to_numpy(), axis=1)
    ]
    return averaged


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--postprocess-only", action="store_true")
    parser.add_argument(
        "--provenance",
        choices=["raw_cel_rma", "geo_deposited_series_matrix"],
        help="Override config modeling.expression_provenance.",
    )
    args = parser.parse_args()
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    modeling = config["modeling"]
    provenance = args.provenance or modeling["expression_provenance"]
    processed_root = ROOT / paths["processed"]
    expression_file = expression_path(processed_root, "GSE44076", "gene", provenance)
    if not expression_file.exists():
        raise FileNotFoundError(f"Missing {provenance} modeling matrix: {expression_file}")

    metrics_root = result_root(ROOT / paths["metrics"], provenance)
    tables_root = result_root(ROOT / paths["tables"], provenance)
    figures_root = result_root(ROOT / paths["figures"], provenance)
    reports_root = result_root(ROOT / paths["reports"], provenance)
    for directory in (metrics_root, tables_root, figures_root, reports_root):
        directory.mkdir(parents=True, exist_ok=True)

    expression = pd.read_parquet(expression_file).set_index("gene_symbol")
    metadata = pd.read_csv(
        ROOT / paths["metadata"] / "gse44076_samples.csv", dtype={"patient_id": str}
    )
    feature_names = expression.index.astype(str).tolist()
    if args.postprocess_only:
        predictions = pd.read_csv(metrics_root / "nested_cv_predictions.csv")
        metrics = pd.read_csv(metrics_root / "nested_cv_metrics.csv")
        assignments = pd.read_csv(metrics_root / "fold_assignments.csv", dtype={"patient_id": str})
        coefficients = pd.read_csv(metrics_root / "elastic_net_fold_coefficients.csv")
    else:
        outputs: dict[str, list[pd.DataFrame]] = {
            "predictions": [],
            "metrics": [],
            "assignments": [],
            "coefficients": [],
        }
        for task, labels in TASKS.items():
            selected = metadata[
                metadata["inclusion_status"].eq("included")
                & metadata["tissue_class"].isin(labels)
            ].copy()
            samples = selected["geo_accession"].tolist()
            result = run_nested_cv(
                expression[samples].T.to_numpy(dtype=np.float32),
                selected["tissue_class"].to_numpy(),
                selected["donor_or_patient_group"].astype(str).to_numpy(),
                selected["geo_accession"].to_numpy(),
                feature_names,
                task,
                ["elastic_net", "linear_svm", "random_forest"],
                modeling["seeds"],
                modeling["outer_splits"],
                modeling["inner_splits"],
                modeling,
            )
            for key, frame in zip(outputs, result, strict=True):
                outputs[key].append(frame)
        predictions = pd.concat(outputs["predictions"], ignore_index=True)
        metrics = pd.concat(outputs["metrics"], ignore_index=True)
        assignments = pd.concat(outputs["assignments"], ignore_index=True)
        coefficients = pd.concat(outputs["coefficients"], ignore_index=True)
        for frame in (predictions, metrics, assignments, coefficients):
            frame["analysis_provenance"] = provenance
        predictions.to_csv(metrics_root / "nested_cv_predictions.csv", index=False)
        metrics.to_csv(metrics_root / "nested_cv_metrics.csv", index=False)
        assignments.to_csv(metrics_root / "fold_assignments.csv", index=False)
        coefficients.to_csv(metrics_root / "elastic_net_fold_coefficients.csv", index=False)
    validate_fold_assignments(assignments)

    total_outer = len(modeling["seeds"]) * int(modeling["outer_splits"])
    stability = feature_stability(
        coefficients,
        total_outer,
        strict_frequency=float(modeling["stable_selection_frequency"]),
        strict_sign_consistency=float(modeling["stable_sign_consistency"]),
    )
    stability["analysis_provenance"] = provenance
    stability.to_csv(tables_root / "feature_stability.csv", index=False)
    coefficient_plot(
        stability,
        "task_c_tumor_vs_adjacent",
        figures_root / "stable_feature_coefficient_plot_task_c",
    )

    summary_rows: list[dict[str, object]] = []
    ci_rows: list[dict[str, object]] = []
    for (task, model), group in metrics.groupby(["task", "model"]):
        candidate_metrics = [
            "f1_macro",
            "balanced_accuracy",
            "log_loss",
            "roc_auc",
            "pr_auc",
            "brier_score",
            "roc_auc_ovr_macro",
            "brier_score_multiclass",
        ]
        available_metrics = [
            metric
            for metric in candidate_metrics
            if metric in group and not group[metric].isna().all()
        ]
        averaged = average_repeated_predictions(
            predictions[predictions["task"].eq(task) & predictions["model"].eq(model)]
        )
        probability_columns = [
            column for column in averaged if column.startswith("probability_")
        ]
        classes = np.asarray(
            [column.removeprefix("probability_") for column in probability_columns]
        )
        aggregated = classification_metrics(
            averaged["y_true"].to_numpy(),
            averaged["y_pred"].to_numpy(),
            averaged[probability_columns].to_numpy(),
            classes,
        )
        bootstrap = grouped_bootstrap_metrics(
            averaged,
            available_metrics,
            int(modeling["bootstrap_iterations"]),
            int(config["project"]["random_seed"]),
        )
        bootstrap_by_metric = {row["metric"]: row for row in bootstrap}
        for row in bootstrap:
            ci_rows.append(
                {"task": task, "model": model, **row, "analysis_provenance": provenance}
            )
        for metric in available_metrics:
            values = group[metric].dropna()
            repeat_means = group.groupby("repeat")[metric].mean().dropna()
            interval = bootstrap_by_metric[metric]
            summary_rows.append(
                {
                    "task": task,
                    "model": model,
                    "metric": metric,
                    "outer_fold_mean": values.mean(),
                    "outer_fold_sd": values.std(ddof=1),
                    "repeat_mean_sd": repeat_means.std(ddof=1),
                    "aggregated_oof_metric": aggregated.get(metric, np.nan),
                    "aggregated_oof_bootstrap_ci_lower": interval["ci_lower"],
                    "aggregated_oof_bootstrap_ci_upper": interval["ci_upper"],
                    "mean": values.mean(),
                    "median": values.median(),
                    "standard_deviation": values.std(ddof=1),
                    "fold_minimum": values.min(),
                    "fold_maximum": values.max(),
                    "repeat_mean_standard_deviation": repeat_means.std(ddof=1),
                    "analysis_provenance": provenance,
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(tables_root / "table_3_model_comparison.csv", index=False)
    pd.DataFrame(ci_rows).to_csv(
        metrics_root / "patient_group_bootstrap_confidence_intervals.csv", index=False
    )
    metrics.to_csv(tables_root / "supplementary_fold_metrics.csv", index=False)
    metrics[["task", "model", "repeat", "outer_fold", "best_parameters", "analysis_provenance"]].to_csv(
        metrics_root / "nested_cv_hyperparameters.csv", index=False
    )

    decisions = []
    full_signatures = []
    for task in TASKS:
        task_summary = summary[
            summary["task"].eq(task) & summary["metric"].eq("f1_macro")
        ].sort_values("mean", ascending=False)
        elastic_mean = float(task_summary.loc[task_summary["model"].eq("elastic_net"), "mean"].iloc[0])
        decisions.append(
            {
                "task": task,
                "selected_model": "elastic_net",
                "elastic_net_mean_macro_f1": elastic_mean,
                "best_comparator_point_estimate": float(task_summary["mean"].iloc[0]),
                "decision_rule": "Elastic Net was prespecified as the primary scientific model; comparators are sensitivity benchmarks.",
                "analysis_provenance": provenance,
            }
        )
        stable = stability[
            stability["task"].eq(task)
            & stability["strictly_stable_gene"]
        ].copy()
        stable["selected_for_full_stable_signature"] = True
        full_signatures.append(stable)
    decision_frame = pd.DataFrame(decisions)
    decision_frame.to_csv(tables_root / "model_selection_decisions.csv", index=False)
    full_signature = pd.concat(full_signatures, ignore_index=True)
    full_signature.to_csv(tables_root / "full_stable_signature.csv", index=False)

    rationale = f"""# Model selection rationale

Elastic Net is the prespecified primary model because the study is a sparse,
interpretable biomarker-discovery analysis. Linear SVM and Random Forest are
benchmarks, not candidates selected by their external-cohort behavior. All
models used the same repeated nested patient/donor-group folds, and all learned
preprocessing and calibration occurred inside training data.

{decision_frame.to_markdown(index=False)}

Expression provenance: `{provenance}`. Compact-panel performance and the final
locked Task C artifact are produced separately by
`scripts/run_compact_panel_analysis.py`.
"""
    (reports_root / "model_selection_rationale.md").write_text(rationale, encoding="utf-8")
    print(
        json.dumps(
            {
                "provenance": provenance,
                "fold_metric_rows": len(metrics),
                "prediction_rows": len(predictions),
                "outer_repeats": len(modeling["seeds"]),
                "models": sorted(metrics["model"].unique()),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
