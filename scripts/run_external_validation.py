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
from sklearn.calibration import calibration_curve
from sklearn.metrics import ConfusionMatrixDisplay, PrecisionRecallDisplay, RocCurveDisplay

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.cross_platform import within_sample_percentile_rank
from src.modeling.evaluation import grouped_bootstrap_metrics
from src.modeling.external_validation import locked_binary_predictions
from src.modeling.pipelines import build_pipeline
from src.visualization.publication_figures import save_figure


def main() -> None:
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    model_path = ROOT / paths["models"] / "task_c_locked_elastic_net.joblib"
    if not model_path.exists():
        raise FileNotFoundError("Run scripts/run_modeling.py before external validation")
    artifact = joblib.load(model_path)
    primary = pd.read_parquet(
        ROOT / paths["processed"] / "GSE44076_gene_expression.parquet"
    ).set_index("gene_symbol")
    external = pd.read_parquet(
        ROOT / paths["processed"] / "GSE41258_gene_expression.parquet"
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
        raise AssertionError("External primary subset contains an invalid tissue")

    signature = artifact["signature_genes"]
    common_signature = [
        gene for gene in signature if gene in primary.index and gene in external.index
    ]
    if len(common_signature) < 2:
        raise RuntimeError("Fewer than two locked signature genes are common to both platforms")
    primary_task = primary_meta[
        primary_meta["tissue_class"].isin(["adjacent_normal", "tumor"])
        & primary_meta["inclusion_status"].eq("included")
    ]
    threshold = float(artifact["threshold"])
    results = []
    predictions = []

    # Representation 1: training-derived z-scaling embedded in the locked model.
    raw_model = artifact["model"]
    if common_signature != signature:
        raw_model = build_pipeline("elastic_net", config["project"]["random_seed"])
        params = {
            **artifact["best_parameters_mode"],
            "variance_quantile__quantile": 0.0,
            "univariate__k": "all",
        }
        raw_model.set_params(**params)
        raw_model.fit(
            primary.loc[common_signature, primary_task["geo_accession"]].T.to_numpy(dtype=np.float32),
            primary_task["tissue_class"].to_numpy(),
        )
    raw_x = external.loc[common_signature, included["geo_accession"]].T.to_numpy(dtype=np.float32)
    raw_pred, raw_metrics = locked_binary_predictions(
        raw_model,
        raw_x,
        included["geo_accession"].to_numpy(),
        included["donor_or_patient_group"].astype(str).to_numpy(),
        included["validation_label"].to_numpy(),
        "training_zscore",
        threshold,
    )
    predictions.append(raw_pred)
    results.append({"representation": "training_zscore", **raw_metrics})

    # Representation 2: ranks are computed within each sample over the locked
    # common-gene universe, with no use of external labels.
    external_gene_set = set(external.index)
    common_universe = [gene for gene in primary.index if gene in external_gene_set]
    signature_positions = [common_universe.index(gene) for gene in common_signature]
    primary_rank_all = within_sample_percentile_rank(
        primary.loc[common_universe, primary_task["geo_accession"]].T.to_numpy(dtype=np.float32)
    )
    external_rank_all = within_sample_percentile_rank(
        external.loc[common_universe, included["geo_accession"]].T.to_numpy(dtype=np.float32)
    )
    rank_model = build_pipeline("elastic_net", config["project"]["random_seed"])
    rank_parameters = {
        **artifact["best_parameters_mode"],
        "variance_quantile__quantile": 0.0,
        "univariate__k": "all",
    }
    rank_model.set_params(**rank_parameters)
    rank_model.fit(
        primary_rank_all[:, signature_positions], primary_task["tissue_class"].to_numpy()
    )
    rank_pred, rank_metrics = locked_binary_predictions(
        rank_model,
        external_rank_all[:, signature_positions],
        included["geo_accession"].to_numpy(),
        included["donor_or_patient_group"].astype(str).to_numpy(),
        included["validation_label"].to_numpy(),
        "within_sample_percentile_rank",
        threshold,
    )
    predictions.append(rank_pred)
    results.append({"representation": "within_sample_percentile_rank", **rank_metrics})
    joblib.dump(
        {
            "model": rank_model,
            "signature_genes": common_signature,
            "common_gene_universe": common_universe,
            "threshold": threshold,
            "provenance": "trained on GSE44076 only; no GSE41258 label tuning",
        },
        ROOT / paths["models"] / "task_c_locked_rank_elastic_net.joblib",
    )

    prediction_frame = pd.concat(predictions, ignore_index=True)
    metric_frame = pd.DataFrame(results)
    ci_rows = []
    for representation, frame in prediction_frame.groupby("representation"):
        summaries = grouped_bootstrap_metrics(
            frame,
            [
                "roc_auc",
                "pr_auc",
                "balanced_accuracy",
                "sensitivity",
                "specificity",
                "f1_macro",
                "brier_score",
            ],
            config["external_validation"]["bootstrap_iterations"],
            config["project"]["random_seed"],
        )
        ci_rows.extend({"representation": representation, **summary} for summary in summaries)
    ci_frame = pd.DataFrame(ci_rows)
    metrics_root = ROOT / paths["metrics"]
    tables_root = ROOT / paths["tables"]
    figures_root = ROOT / paths["figures"]
    metrics_root.mkdir(parents=True, exist_ok=True)
    tables_root.mkdir(parents=True, exist_ok=True)
    prediction_frame.to_csv(metrics_root / "external_validation_predictions.csv", index=False)
    metric_frame.to_csv(metrics_root / "external_validation_metrics.csv", index=False)
    ci_frame.to_csv(metrics_root / "external_validation_bootstrap_ci.csv", index=False)
    metric_frame.assign(
        included_samples=len(included),
        excluded_samples=len(external_meta) - len(included),
        common_signature_genes=len(common_signature),
        threshold=threshold,
    ).to_csv(tables_root / "table_5_external_validation.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for representation, frame in prediction_frame.groupby("representation"):
        positive = frame["y_true"].eq("tumor").astype(int)
        probability = frame["probability_tumor"]
        RocCurveDisplay.from_predictions(
            positive, probability, name=representation, ax=axes[0, 0]
        )
        PrecisionRecallDisplay.from_predictions(
            positive, probability, name=representation, ax=axes[0, 1]
        )
        observed, predicted = calibration_curve(positive, probability, n_bins=8, strategy="quantile")
        axes[1, 0].plot(predicted, observed, marker="o", label=representation)
    axes[1, 0].plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.8)
    axes[1, 0].set(xlabel="Mean predicted tumor probability", ylabel="Observed tumor fraction")
    axes[1, 0].set_title("External calibration")
    axes[1, 0].legend(frameon=False)
    preferred = prediction_frame[
        prediction_frame["representation"].eq("within_sample_percentile_rank")
    ]
    ConfusionMatrixDisplay.from_predictions(
        preferred["y_true"],
        preferred["y_pred"],
        labels=["adjacent_normal", "tumor"],
        display_labels=["Normal colon", "Primary tumor"],
        colorbar=False,
        ax=axes[1, 1],
    )
    axes[1, 1].set_title("Locked rank representation")
    fig.suptitle(
        f"GSE41258 external validation (n={len(included)}; {len(common_signature)} genes)"
    )
    save_figure(fig, figures_root / "GSE41258_external_validation_performance")

    report = {
        "included": int(len(included)),
        "excluded": int(len(external_meta) - len(included)),
        "normal_colon": int(included["tissue_class"].eq("normal_colon").sum()),
        "primary_tumor": int(included["tissue_class"].eq("primary_tumor").sum()),
        "signature_genes_primary": len(signature),
        "signature_genes_common": len(common_signature),
        "technical_replicates_excluded": int(
            external_meta["technical_replicate_candidate"].sum()
        ),
        "metrics": results,
        "claim_boundary": (
            "Validates the tumor-versus-normal-colon component only; it is not independent "
            "validation of healthy-versus-adjacent field cancerization."
        ),
    }
    (ROOT / paths["reports"] / "external_validation_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
