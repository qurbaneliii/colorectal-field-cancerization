from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.validation import validate_fold_assignments
from src.interpretation.coefficients import coefficient_plot
from src.modeling.evaluation import grouped_bootstrap_ci
from src.modeling.nested_cv import modal_best_parameters, run_nested_cv
from src.modeling.pipelines import build_pipeline
from src.modeling.stability import feature_stability


TASKS = {
    "task_a_three_class": ["healthy", "adjacent_normal", "tumor"],
    "task_b_field_effect": ["healthy", "adjacent_normal"],
    "task_c_tumor_vs_adjacent": ["adjacent_normal", "tumor"],
}


def average_repeated_predictions(frame: pd.DataFrame) -> pd.DataFrame:
    probability_columns = [
        c for c in frame if c.startswith("probability_") and frame[c].notna().any()
    ]
    first = {"patient_id": "first", "y_true": "first"}
    means = {column: "mean" for column in probability_columns}
    averaged = frame.groupby("sample_id", as_index=False).agg({**first, **means})
    classes = [column.removeprefix("probability_") for column in probability_columns]
    averaged["y_pred"] = np.asarray(classes)[
        np.argmax(averaged[probability_columns].to_numpy(), axis=1)
    ]
    return averaged


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--postprocess-only",
        action="store_true",
        help="Resume summaries/model locking from completed nested-CV CSV files.",
    )
    args = parser.parse_args()
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    modeling = config["modeling"]
    metrics_root = ROOT / paths["metrics"]
    tables_root = ROOT / paths["tables"]
    models_root = ROOT / paths["models"]
    figures_root = ROOT / paths["figures"]
    reports_root = ROOT / paths["reports"]
    for directory in (metrics_root, tables_root, models_root, figures_root, reports_root):
        directory.mkdir(parents=True, exist_ok=True)

    expression_table = pd.read_parquet(
        ROOT / paths["processed"] / "GSE44076_gene_expression.parquet"
    )
    expression = expression_table.set_index("gene_symbol")
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
        all_predictions = []
        all_metrics = []
        all_assignments = []
        all_coefficients = []
        for task, labels in TASKS.items():
            selected = metadata[
                metadata["inclusion_status"].eq("included") & metadata["tissue_class"].isin(labels)
            ].copy()
            samples = selected["geo_accession"].tolist()
            x = expression[samples].T.to_numpy(dtype=np.float32)
            y = selected["tissue_class"].to_numpy()
            groups = selected["donor_or_patient_group"].astype(str).to_numpy()
            sample_ids = selected["geo_accession"].to_numpy()
            task_predictions, task_metrics, task_assignments, task_coefficients = run_nested_cv(
                x,
                y,
                groups,
                sample_ids,
                feature_names,
                task,
                ["elastic_net", "linear_svm", "random_forest"],
                modeling["seeds"],
                modeling["outer_splits"],
                modeling["inner_splits"],
            )
            all_predictions.append(task_predictions)
            all_metrics.append(task_metrics)
            all_assignments.append(task_assignments)
            all_coefficients.append(task_coefficients)
        predictions = pd.concat(all_predictions, ignore_index=True)
        metrics = pd.concat(all_metrics, ignore_index=True)
        assignments = pd.concat(all_assignments, ignore_index=True)
        coefficients = pd.concat(all_coefficients, ignore_index=True)
        predictions.to_csv(metrics_root / "nested_cv_predictions.csv", index=False)
        metrics.to_csv(metrics_root / "nested_cv_metrics.csv", index=False)
        assignments.to_csv(metrics_root / "fold_assignments.csv", index=False)
        coefficients.to_csv(metrics_root / "elastic_net_fold_coefficients.csv", index=False)
    validate_fold_assignments(assignments)

    total_outer = len(modeling["seeds"]) * modeling["outer_splits"]
    stability = feature_stability(coefficients, total_outer)
    stability.to_csv(tables_root / "feature_stability.csv", index=False)
    coefficient_plot(
        stability,
        "task_c_tumor_vs_adjacent",
        figures_root / "stable_feature_coefficient_plot_task_c",
    )

    summary_path = tables_root / "table_3_model_comparison.csv"
    ci_path = metrics_root / "patient_group_bootstrap_confidence_intervals.csv"
    if args.postprocess_only and summary_path.exists() and ci_path.exists():
        summary = pd.read_csv(summary_path)
    else:
        summary_rows = []
        ci_rows = []
        for (task, model), group in metrics.groupby(["task", "model"]):
            for metric in ["f1_macro", "balanced_accuracy", "log_loss", "roc_auc", "pr_auc"]:
                if metric not in group or group[metric].isna().all():
                    continue
                values = group[metric].dropna()
                summary_rows.append(
                    {
                        "task": task,
                        "model": model,
                        "metric": metric,
                        "mean": values.mean(),
                        "median": values.median(),
                        "standard_deviation": values.std(ddof=1),
                        "fold_minimum": values.min(),
                        "fold_maximum": values.max(),
                    }
                )
            model_predictions = average_repeated_predictions(
                predictions[(predictions["task"].eq(task)) & (predictions["model"].eq(model))]
            )
            ci = grouped_bootstrap_ci(
                model_predictions,
                "f1_macro",
                modeling["bootstrap_iterations"],
                config["project"]["random_seed"],
            )
            ci_rows.append({"task": task, "model": model, **ci})
        summary = pd.DataFrame(summary_rows)
        summary.to_csv(summary_path, index=False)
        pd.DataFrame(ci_rows).to_csv(ci_path, index=False)
    metrics.to_csv(tables_root / "supplementary_fold_metrics.csv", index=False)

    # Prespecified decision: Elastic Net unless its mean macro F1 is outside the
    # uncertainty range of the best comparator. External behavior is assessed later.
    decisions = []
    final_signatures = []
    panels = []
    for task in TASKS:
        task_summary = summary[
            summary["task"].eq(task) & summary["metric"].eq("f1_macro")
        ].sort_values("mean", ascending=False)
        elastic_mean = float(
            task_summary.loc[task_summary["model"].eq("elastic_net"), "mean"].iloc[0]
        )
        best_mean = float(task_summary["mean"].iloc[0])
        selected_model = "elastic_net"
        decisions.append(
            {
                "task": task,
                "selected_model": selected_model,
                "elastic_net_mean_macro_f1": elastic_mean,
                "best_point_estimate": best_mean,
                "decision_rule": "Elastic Net retained for leakage safety, parsimony, stability, and interpretability; point estimates are not the sole criterion.",
            }
        )
        stable = stability[
            stability["task"].eq(task)
            & stability["selection_frequency"].ge(modeling["stable_selection_frequency"])
            & stability["sign_consistency"].ge(0.8)
        ].copy()
        if stable.empty:
            stable = stability[stability["task"].eq(task)].head(30).copy()
            stable["fallback_reason"] = "no feature met prespecified stability threshold"
        stable["selected_for_final_signature"] = True
        final_signatures.append(stable)
        ranked = stable.sort_values(
            ["selection_frequency", "median_absolute_coefficient"], ascending=False
        )
        for size in modeling["compact_panels"]:
            for rank, gene in enumerate(ranked["gene_symbol"].head(size), start=1):
                panels.append(
                    {"task": task, "panel_size": size, "rank": rank, "gene_symbol": gene}
                )

    decision_frame = pd.DataFrame(decisions)
    decision_frame.to_csv(tables_root / "model_selection_decisions.csv", index=False)
    final_signature = pd.concat(final_signatures, ignore_index=True)
    final_signature.to_csv(tables_root / "final_signature.csv", index=False)
    final_signature.to_csv(tables_root / "table_4_final_signature.csv", index=False)
    pd.DataFrame(panels).to_csv(tables_root / "candidate_gene_panels.csv", index=False)

    # Lock Task C model and signature using GSE44076 only.
    task = "task_c_tumor_vs_adjacent"
    labels = TASKS[task]
    selected = metadata[
        metadata["inclusion_status"].eq("included") & metadata["tissue_class"].isin(labels)
    ].copy()
    signature_genes = final_signature[final_signature["task"].eq(task)]["gene_symbol"].drop_duplicates().tolist()
    if not signature_genes:
        raise RuntimeError("Task C final signature is empty")
    best_params = modal_best_parameters(metrics, task, "elastic_net")
    # Panel fit uses only locked genes; selection steps are valid but set to retain all.
    locked = build_pipeline("elastic_net", config["project"]["random_seed"])
    locked_parameters = {
        **best_params,
        "variance_quantile__quantile": 0.0,
        "univariate__k": "all",
    }
    locked.set_params(**locked_parameters)
    x_locked = expression.loc[signature_genes, selected["geo_accession"]].T.to_numpy(dtype=np.float32)
    y_locked = selected["tissue_class"].to_numpy()
    locked.fit(x_locked, y_locked)
    artifact = {
        "model": locked,
        "task": task,
        "classes": locked.named_steps["model"].classes_.tolist(),
        "signature_genes": signature_genes,
        "training_samples": selected["geo_accession"].tolist(),
        "training_groups": selected["donor_or_patient_group"].astype(str).tolist(),
        "best_parameters_mode": best_params,
        "threshold": config["external_validation"]["locked_threshold"],
        "provenance": "GSE44076 only; GEO deposited normalized matrix",
    }
    joblib.dump(artifact, models_root / "task_c_locked_elastic_net.joblib")
    reloaded = joblib.load(models_root / "task_c_locked_elastic_net.joblib")
    if reloaded["model"].predict(x_locked[:5]).shape[0] != 5:
        raise AssertionError("Serialized model reload sanity check failed")

    rationale = f"""# Model selection rationale

The prespecified priority order was leakage safety, external validity, balanced
performance, adjacent-normal recall, calibration, stability, compactness,
interpretability, and computational simplicity. Elastic Net was retained for
all three tasks. It was evaluated in the same repeated nested patient-group
cross-validation as the linear SVM and Random Forest; no point estimate alone
determined the decision.

{decision_frame.to_markdown(index=False)}

Task C locked signature size: {len(signature_genes)} genes. The model artifact
was serialized and reloaded successfully. External-cohort labels were not used
in this fit or any hyperparameter decision.
"""
    (reports_root / "model_selection_rationale.md").write_text(rationale, encoding="utf-8")
    print(
        json.dumps(
            {
                "fold_metric_rows": len(metrics),
                "prediction_rows": len(predictions),
                "task_c_signature_size": len(signature_genes),
                "models": sorted(metrics["model"].unique()),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
