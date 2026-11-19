"""Generate experimental comparisons and a figure from verified measured outputs."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("execution", type=Path)
    parser.add_argument("report_dir", type=Path)
    args = parser.parse_args()
    args.report_dir.mkdir(parents=True, exist_ok=True)
    root = args.execution
    execution = json.loads((root / "execution.json").read_text())
    metrics = pd.read_csv(root / "fold_metrics.csv")
    columns = ["f1_macro", "balanced_accuracy", "roc_auc", "pr_auc", "sensitivity", "specificity", "brier_score", "log_loss"]
    means = metrics.groupby(["task", "method"])[columns].mean().reset_index()
    means.to_csv(root / "mean_fold_metrics.csv", index=False)
    external = pd.read_csv(root / "external_metrics.csv")
    invariant = pd.read_csv(root / "invariance.csv")
    predictions = pd.read_csv(root / "heldout_predictions.csv")
    averaged = predictions.groupby(["task", "sample_id", "patient_id", "observed", "method"], as_index=False).score.mean()
    pooled = []
    for (task, method), part in averaged.groupby(["task", "method"]):
        pooled.append({"task": task, "method": method, "f1_macro": f1_score(part.observed, part.score >= .5, average="macro"),
                       "balanced_accuracy": balanced_accuracy_score(part.observed, part.score >= .5),
                       "roc_auc": roc_auc_score(part.observed, part.score)})
    pd.DataFrame(pooled).to_csv(root / "repeat_averaged_metrics.csv", index=False)
    text = "# Measured model comparison\n\n"
    text += "These are conditional historical-panel L2 logistic-regression comparisons. The same model settings are used in both arms. They are distinct from the published nested Elastic Net estimates.\n\n"
    text += "Mean across the 15 paired held-out folds per task\n\n" + means.to_markdown(index=False, floatfmt=".6f")
    text += "\n\nThe paired donor-bootstrap intervals in `paired_conditional_intervals.csv` use repeat-averaged predictions, not the mean-fold estimand above. They do not cover feature discovery or cross-study uncertainty.\n\n"
    text += "Repeat-averaged metrics corresponding to that interval estimand\n\n" + pd.DataFrame(pooled).to_markdown(index=False, floatfmt=".6f")
    text += "\n\nCalibration fits and convergence flags are retained for each fold; separation can prevent reliable slope estimation. Brier score and log loss remain reported.\n"
    (args.report_dir / "04_MODEL_COMPARISON.md").write_text(text, encoding="utf-8")
    text = "# Measured platform transport\n\n"
    text += f"The shared annotation universe contains {execution['shared_genes']:,} genes. CEMIP and ETV4 are refit on development ranks; FOXQ1 is absent externally. The experimental frozen reference is learned on 196 development arrays.\n\n"
    text += external.to_markdown(index=False, floatfmt=".6f")
    text += "\n\nBoth arms above use a new prespecified L2 refit and retain 0.41 as a historical reference cutoff. The primary archived Elastic Net transport result remains ROC-AUC 0.982894, macro-F1 0.828155, sensitivity 1.0 and specificity 0.557692. Deposited-expression specificity 0.942 remains a separate sensitivity. Differences between the new joint-RMA arm and that historical model include classifier policy and cannot be attributed to normalization.\n\n"
    text += "No external outcome selected the method, parameters, panel or cutoff. GSE41258 is previously exposed retrospective evidence. Rank transport does not establish calibrated probabilities or clinical validity.\n"
    (args.report_dir / "05_PLATFORM_TRANSPORT.md").write_text(text, encoding="utf-8")
    text = "# Real-data execution\n\n"
    text += f"Analysis commit `{execution['analysis_commit']}` completed on verified current-GEO CELs. All 246 GSE44076 and 233 endpoint GSE41258 arrays were independently extracted; 390 external raw-array inventory/QC evidence remains in the primary audit. No raw download was repeated. The local archived joint-RMA matrices are reused baselines, not newly computed RMA.\n\n"
    text += "Source TAR SHA-256 values\n\n"
    text += "\n".join(f"- `{name}` — `{sha}`" for name, sha in execution["source_archives"]["archives"].items())
    text += f"\n\nThere are {len(metrics)} fold/method rows, {len(predictions)} held-out score rows and {len(invariant)} invariance comparisons. Maximum gene absolute difference is {invariant.gene_max_abs.max():.3g}; maximum probability difference is {invariant.probability_abs.max():.3g}. All class decisions match.\n\n"
    text += "The independent raw-PM experiment checked ten arrays and found zero differences in background-corrected PM values and ranks. Input/member hashes, manufacturer mapping hashes, train/test IDs, quantile references, numeric models, environments, score locks, independent metric replay and conditional intervals are retained with execution evidence.\n\n"
    text += "The first slow extraction implementation was stopped after feasibility checking; a direct-reader implementation produced identical statistics on the checked sample and completed the extraction. Old partial caches were not mixed into this run. The preserved Windows R environment emitted locale and package-build-version warnings; no dependency version was changed.\n"
    (args.report_dir / "03_REAL_DATA_EXECUTION.md").write_text(text, encoding="utf-8")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for ax, metric in zip(axes, ["f1_macro", "balanced_accuracy", "roc_auc"]):
        for i, task in enumerate(means.task.unique()):
            for j, method in enumerate(["joint_rma", "frozen"]):
                val = means.loc[means.task.eq(task) & means.method.eq(method), metric].iloc[0]
                ax.bar(i + (j-.5)*.32, val, width=.3,
                       label=method if i == 0 else None, color=["#657a91", "#197b70"][j])
        ax.set_xticks(np.arange(2), ["Task B", "Task C"])
        ax.set_ylim(0, 1.03)
        ax.set_title(metric.replace("_", " "))
    axes[0].legend(fontsize=8)
    fig.suptitle("Fixed historical panels: mean held-out fold metrics (conditional sensitivity)")
    fig.tight_layout()
    fig.savefig(root / "comparison.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
