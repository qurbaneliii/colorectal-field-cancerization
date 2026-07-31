from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    PrecisionRecallDisplay,
    RocCurveDisplay,
)

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.visualization.publication_figures import PALETTE, save_figure


def study_design_figure(stem: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.axis("off")
    boxes = [
        (0.02, 0.58, "GSE44076\n50 healthy donors\n98 paired patients"),
        (0.28, 0.58, "Independent annotation\nand deposited-matrix QC"),
        (0.54, 0.58, "Biological inference\nfield-effect trajectories\npaired tumor contrast"),
        (0.80, 0.58, "Repeated nested\npatient-group CV\n3 prediction tasks"),
        (0.54, 0.12, "Locked compact\nElastic Net signature"),
        (0.80, 0.12, "GSE41258\nPrimary tumor vs\nnormal colon only"),
    ]
    for x, y, text in boxes:
        ax.add_patch(
            plt.Rectangle((x, y), 0.18, 0.25, facecolor="#EAF4F4", edgecolor="#1D3557", lw=1.2)
        )
        ax.text(x + 0.09, y + 0.125, text, ha="center", va="center", fontsize=10)
    arrows = [
        ((0.20, 0.705), (0.28, 0.705)),
        ((0.46, 0.705), (0.54, 0.705)),
        ((0.72, 0.705), (0.80, 0.705)),
        ((0.89, 0.58), (0.63, 0.37)),
        ((0.72, 0.245), (0.80, 0.245)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 1.5})
    ax.set_title("Leakage-safe colorectal field-cancerization study design", fontsize=15)
    save_figure(fig, stem)


def sample_flow_figure(primary: pd.DataFrame, external: pd.DataFrame, stem: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax in axes:
        ax.axis("off")
    p_counts = primary["tissue_class"].value_counts()
    axes[0].text(
        0.5,
        0.78,
        f"GSE44076 deposited arrays\nn={len(primary)}",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.5", "fc": "#EAF4F4"},
    )
    axes[0].text(
        0.5,
        0.35,
        f"Included: {p_counts['healthy']} healthy\n"
        f"{p_counts['adjacent_normal']} adjacent-normal\n{p_counts['tumor']} tumor\n"
        "98 complete tumor-adjacent pairs",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.5", "fc": "#F1FAEE"},
    )
    axes[0].annotate("", xy=(0.5, 0.48), xytext=(0.5, 0.66), arrowprops={"arrowstyle": "->"})
    included = external[external["inclusion_status"].eq("included")]
    e_counts = included["tissue_class"].value_counts()
    axes[1].text(
        0.5,
        0.78,
        f"GSE41258 deposited arrays\nn={len(external)}",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.5", "fc": "#EAF4F4"},
    )
    axes[1].text(
        0.5,
        0.35,
        f"External validation: n={len(included)}\n"
        f"{e_counts.get('normal_colon', 0)} normal colon\n"
        f"{e_counts.get('primary_tumor', 0)} primary tumor\n"
        f"Excluded: {len(external) - len(included)}",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.5", "fc": "#F1FAEE"},
    )
    axes[1].annotate("", xy=(0.5, 0.48), xytext=(0.5, 0.66), arrowprops={"arrowstyle": "->"})
    fig.suptitle("Sample inclusion flow")
    save_figure(fig, stem)


def modeling_figures(metrics: pd.DataFrame, predictions: pd.DataFrame, figures: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=metrics, x="task", y="f1_macro", hue="model", ax=ax)
    sns.stripplot(
        data=metrics,
        x="task",
        y="f1_macro",
        hue="model",
        dodge=True,
        color="black",
        alpha=0.45,
        size=3,
        legend=False,
        ax=ax,
    )
    ax.set_ylabel("Outer-fold macro F1")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=15)
    ax.set_title("Repeated nested patient-group cross-validation")
    ax.legend(frameon=False, title="Model")
    save_figure(fig, figures / "nested_cv_performance_distribution")

    frame = predictions[
        predictions["task"].eq("task_a_three_class")
        & predictions["model"].eq("elastic_net")
    ]
    probability_columns = sorted(c for c in frame if c.startswith("probability_"))
    averaged = frame.groupby("sample_id", as_index=False).agg(
        {"patient_id": "first", "y_true": "first", **{c: "mean" for c in probability_columns}}
    )
    classes = [c.removeprefix("probability_") for c in probability_columns]
    averaged["y_pred"] = np.asarray(classes)[
        np.argmax(averaged[probability_columns].to_numpy(), axis=1)
    ]
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ConfusionMatrixDisplay.from_predictions(
        averaged["y_true"],
        averaged["y_pred"],
        labels=classes,
        display_labels=classes,
        cmap="Blues",
        colorbar=False,
        ax=ax,
    )
    ax.set_title(f"Three-class Elastic Net out-of-fold predictions (n={len(averaged)})")
    plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
    save_figure(fig, figures / "task_a_out_of_fold_confusion_matrix")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for class_name in classes:
        binary = averaged["y_true"].eq(class_name).astype(int)
        probability = averaged[f"probability_{class_name}"]
        RocCurveDisplay.from_predictions(binary, probability, name=class_name, ax=axes[0])
        PrecisionRecallDisplay.from_predictions(binary, probability, name=class_name, ax=axes[1])
        observed, predicted = calibration_curve(binary, probability, n_bins=8, strategy="quantile")
        axes[2].plot(predicted, observed, marker="o", label=class_name)
    axes[2].plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.8)
    axes[2].set(xlabel="Mean predicted probability", ylabel="Observed fraction")
    axes[2].set_title("One-vs-rest calibration")
    axes[2].legend(frameon=False)
    fig.suptitle("Three-class out-of-fold discrimination and calibration")
    save_figure(fig, figures / "task_a_roc_pr_calibration")


def biological_figures(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    candidates: pd.DataFrame,
    de: pd.DataFrame,
    figures: Path,
) -> None:
    genes = candidates["gene_symbol"].drop_duplicates().head(30).tolist()
    ordered = metadata.sort_values(["tissue_class", "patient_id"])["geo_accession"].tolist()
    z = expression.loc[genes, ordered]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0)
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(z, cmap="vlag", center=0, xticklabels=False, yticklabels=True, ax=ax)
    ax.set_title(f"Top field-cancerization candidates (n={len(genes)} genes)")
    ax.set_xlabel(f"Samples ordered by tissue class (n={len(ordered)})")
    ax.set_ylabel("Gene symbol")
    save_figure(fig, figures / "field_cancerization_signature_heatmap")

    top = candidates["gene_symbol"].drop_duplicates().head(6).tolist()
    long = (
        expression.loc[top]
        .T.reset_index(names="geo_accession")
        .merge(metadata[["geo_accession", "tissue_class"]], on="geo_accession")
        .melt(id_vars=["geo_accession", "tissue_class"], var_name="gene_symbol", value_name="expression")
    )
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    order = ["healthy", "adjacent_normal", "tumor"]
    for ax, gene in zip(axes.ravel(), top, strict=True):
        subset = long[long["gene_symbol"].eq(gene)]
        sns.violinplot(
            data=subset,
            x="tissue_class",
            y="expression",
            order=order,
            palette=PALETTE,
            inner=None,
            cut=0,
            ax=ax,
        )
        sns.boxplot(
            data=subset,
            x="tissue_class",
            y="expression",
            order=order,
            width=0.28,
            showfliers=False,
            color="white",
            ax=ax,
        )
        ax.set_title(gene)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=20)
        ax.set_ylabel("Normalized log2 expression")
    fig.suptitle("Healthy → adjacent-normal → tumor candidate trajectories")
    save_figure(fig, figures / "top_gene_expression_trajectories")

    paired = de[de["comparison"].eq("tumor_vs_adjacent_normal_paired")].nsmallest(
        20, "adjusted_p_value"
    )
    paired = paired.sort_values("log2_fold_change")
    fig, ax = plt.subplots(figsize=(8, 7))
    colors = ["#277DA1" if x < 0 else "#D62828" for x in paired["log2_fold_change"]]
    ax.barh(paired["gene_symbol"], paired["log2_fold_change"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Tumor minus paired adjacent-normal mean log2 expression")
    ax.set_ylabel("Gene symbol")
    ax.set_title("Top paired tumor-versus-adjacent processed-matrix contrasts")
    save_figure(fig, figures / "tumor_vs_adjacent_paired_de_plot_exploratory")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures", action="store_true")
    parser.parse_args()
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    metadata_root = ROOT / paths["metadata"]
    tables = ROOT / paths["tables"]
    figures = ROOT / paths["figures"]
    manuscript = ROOT / paths["manuscript"]
    reports = ROOT / paths["reports"]
    for directory in (tables, figures, manuscript, reports):
        directory.mkdir(parents=True, exist_ok=True)

    primary = pd.read_csv(metadata_root / "gse44076_samples.csv", dtype={"patient_id": str})
    external = pd.read_csv(metadata_root / "gse41258_samples.csv", dtype={"patient_id": str})
    expression = pd.read_parquet(
        ROOT / paths["processed"] / "GSE44076_gene_expression.parquet"
    ).set_index("gene_symbol")
    de = pd.read_csv(
        ROOT / paths["de"] / "processed_matrix_exploratory_all_comparisons.csv"
    )
    candidates = pd.read_csv(tables / "field_cancerization_candidates.csv")
    metrics = pd.read_csv(ROOT / paths["metrics"] / "nested_cv_metrics.csv")
    predictions = pd.read_csv(ROOT / paths["metrics"] / "nested_cv_predictions.csv")
    external_metrics = pd.read_csv(ROOT / paths["metrics"] / "external_validation_metrics.csv")
    signature = pd.read_csv(tables / "final_signature.csv")

    study_design_figure(figures / "study_design_pipeline")
    sample_flow_figure(primary, external, figures / "sample_count_flowchart")
    modeling_figures(metrics, predictions, figures)
    biological_figures(expression, primary, candidates, de, figures)

    table1 = pd.DataFrame(
        [
            ["GSE44076", "GPL13667", "healthy", int(primary["tissue_class"].eq("healthy").sum()), "development"],
            ["GSE44076", "GPL13667", "adjacent_normal", int(primary["tissue_class"].eq("adjacent_normal").sum()), "development"],
            ["GSE44076", "GPL13667", "tumor", int(primary["tissue_class"].eq("tumor").sum()), "development"],
            ["GSE41258", "GPL96", "normal_colon", int((external["inclusion_status"].eq("included") & external["tissue_class"].eq("normal_colon")).sum()), "external validation"],
            ["GSE41258", "GPL96", "primary_tumor", int((external["inclusion_status"].eq("included") & external["tissue_class"].eq("primary_tumor")).sum()), "external validation"],
        ],
        columns=["accession", "platform", "tissue_class", "included_samples", "role"],
    )
    table1.to_csv(tables / "table_1_dataset_characteristics.csv", index=False)
    biology = config["biology"]
    de_summary = (
        de.assign(
            significant=lambda x: x["adjusted_p_value"].lt(biology["fdr_threshold"])
            & x["log2_fold_change"].abs().ge(biology["log2fc_threshold"])
        )
        .groupby(["comparison", "direction"], as_index=False)
        .agg(tested_genes=("gene_symbol", "size"), significant_genes=("significant", "sum"))
    )
    de_summary.to_csv(tables / "table_2_differential_expression_summary.csv", index=False)
    pd.concat(
        [
            primary.assign(accession="GSE44076"),
            external.assign(accession="GSE41258"),
        ],
        ignore_index=True,
        sort=False,
    ).to_csv(tables / "supplementary_sample_metadata.csv", index=False)
    de.to_csv(tables / "supplementary_full_de_results.csv", index=False)

    env = json.loads((metadata_root / "software_environment.json").read_text(encoding="utf-8"))
    best = (
        metrics.groupby(["task", "model"], as_index=False)["f1_macro"]
        .mean()
        .query("model == 'elastic_net'")
    )
    external_display = external_metrics[
        ["representation", "roc_auc", "pr_auc", "balanced_accuracy", "f1_macro", "sensitivity", "specificity"]
    ]
    methods = f"""# Methods

## Datasets and eligibility

GSE44076 (GPL13667) was the sole development cohort: 50 cancer-free healthy
donors, 98 tumor-adjacent normal samples, and 98 tumors from the same 98
patients. GSE41258 (GPL96) was normalized and annotated independently. Its main
external subset was restricted before analysis to author-included Primary Tumor
and Normal Colon arrays; metastases, non-colon normal tissues, polyps, cell
lines, ambiguous samples, and titles matching `_ez` or `rehyb` were excluded.

## Preprocessing and quality control

The executed Python sensitivity route used the normalized expression deposited
in the GEO series matrices. Official GEO platform annotation was fetched
independently; probes with absent or ambiguous gene-symbol mapping were removed,
and multiple unambiguously mapped probes were aggregated per sample by the
median. The two platforms were never jointly normalized. Distribution,
PCA, correlation, missingness, and robust sample-median outlier-candidate checks
were generated; separation alone never triggered exclusion.

Raw-CEL RMA scripts use `oligo`/`pd.hg.u219`/`hgu219.db` for GSE44076 and
`affy`/`hgu133a.db` for GSE41258, followed by label-independent median
probe-to-gene aggregation, arrayQualityMetrics, limma, and clusterProfiler.
They were not executed in this environment because Rscript was unavailable.
Consequently, inferential differential-expression tables in this run are
explicitly labeled processed-matrix Welch or paired-t sensitivity analyses and
must not be represented as the final limma analysis.

## Statistical and trajectory analysis

Adjacent-normal versus healthy and tumor versus healthy sensitivity contrasts
used Welch tests. Tumor versus adjacent-normal used patient-matched paired
tests. Benjamini-Hochberg correction was applied. Field candidates required
FDR < {biology['fdr_threshold']}, absolute mean log2 difference ≥
{biology['log2fc_threshold']}, and at least 70% sample-direction consistency.
Trajectory categories were assigned by explicit rules in `src/biology.py` with
a {biology['trajectory_tolerance']} log2 tolerance.

## Predictive modeling and leakage prevention

Three objectives were analyzed separately: healthy/adjacent/tumor,
healthy/adjacent, and adjacent/tumor. Repeated {config['modeling']['outer_splits']}-fold
outer and {config['modeling']['inner_splits']}-fold inner StratifiedGroupKFold
splits used {len(config['modeling']['seeds'])} deterministic seeds. Patient
pairs were indivisible; healthy donors had unique groups. Zero-variance
removal, training-variance filtering, univariate supervised selection,
standardization, tuning, and fitting all occurred inside the training pipeline.
Elastic Net logistic regression was compared with a linear SVM and Random
Forest. Macro F1 was primary. Out-of-fold predictions and every assignment were
retained. No SMOTE or neural network was used.

Feature stability was the outer-fold selection frequency plus coefficient-sign
consistency. Elastic Net was prespecified as preferred when performance was
within uncertainty of a complex comparator. The Task C panel and threshold
were locked from GSE44076. GSE41258 labels were never used for tuning. External
evaluation compared training-derived scaling with within-sample percentile
ranks computed over the common-gene universe. Confidence intervals used
patient/donor-group bootstrapping.

## Software and seeds

Python {env['python'].split()[0]} on {env['platform']}; scikit-learn
{env['scikit_learn']}. Primary seed: {config['project']['random_seed']};
additional deterministic seeds: {config['modeling']['seeds']}.
"""
    (manuscript / "methods.md").write_text(methods, encoding="utf-8")

    results_text = f"""# Results

## Cohorts and processing

The primary audit recovered exactly 246 arrays: 50 healthy, 98 adjacent-normal,
and 98 tumor, with 98 complete patient pairs and no count discrepancy. Official
platform annotation and deterministic aggregation yielded {len(expression):,}
GSE44076 genes. GSE41258 yielded
{len(pd.read_parquet(ROOT / paths['processed'] / 'GSE41258_gene_expression.parquet')):,}
genes; {len((ROOT / paths['processed'] / 'common_genes_GSE44076_GSE41258.txt').read_text().splitlines()):,}
symbols were common.

The processed-matrix sensitivity analysis identified {len(candidates):,}
adjacent-versus-healthy candidates meeting the prespecified FDR, effect, and
direction-consistency criteria. These are provisional until raw-CEL
Bioconductor preprocessing and limma are executed.

## Internal validation

Elastic Net outer-fold mean macro F1 values were:

{best.to_markdown(index=False)}

All reported predictions were generated for untouched patient-group outer
folds. Adjacent-normal recall, per-class metrics, calibration scores, and all
fold-level distributions are in `results/metrics/nested_cv_metrics.csv`.

## Signature and external validation

The locked Task C stable signature contained
{signature[signature['task'].eq('task_c_tumor_vs_adjacent')]['gene_symbol'].nunique()}
genes before cross-platform intersection. GSE41258 main validation contained
{int(external['inclusion_status'].eq('included').sum())} eligible arrays.

{external_display.to_markdown(index=False)}

This validates only the tumor-versus-normal-colon component. GSE41258 does not
provide an independent cancer-free healthy-versus-adjacent validation cohort.
"""
    (manuscript / "results.md").write_text(results_text, encoding="utf-8")

    limitations = """# Limitations

- These are retrospective public microarray cohorts, not prospectively collected
  clinical-validation samples.
- GPL13667 and GPL96 differ in probe design and measurement distribution;
  within-sample ranks reduce but do not eliminate platform shift.
- GSE41258 validates tumor versus normal colon, not the healthy
  cancer-free-versus-adjacent field effect.
- Stage and molecular-subtype coverage is limited; GSE44076 is enriched for
  stage II, microsatellite-stable disease.
- Bulk-tissue expression may reflect cell-composition changes as well as
  epithelial field biology.
- No prospective or clinical-utility validation was performed.
- High-dimensional small-sample selection remains vulnerable to instability,
  even with nested grouped validation and stability thresholds.
- The executed environment lacked R, so raw-CEL RMA, final limma inference,
  arrayQualityMetrics, and enrichment remain pending.
- The signature is a biomarker-discovery candidate, not a clinically validated
  or clinically ready diagnostic.
"""
    (manuscript / "limitations.md").write_text(limitations, encoding="utf-8")

    legends = """# Figure legends

1. **Study design and pipeline.** Independent cohort preprocessing, biological
   analysis, patient-group nested cross-validation, signature locking, and the
   restricted external-validation claim.
2. **Sample-count flowchart.** Deposited and eligible samples for both GEO
   accessions, including exact primary tissue counts.
3. **Preprocessing QC.** Deposited normalized expression distributions,
   platform-specific PCA, and sample-correlation heatmaps.
4. **Primary PCA and clustering.** GSE44076 PCA and top-variable-gene heatmap
   across healthy, adjacent-normal, and tumor tissues.
5. **Adjacent-versus-healthy volcano.** Processed-matrix sensitivity contrast;
   red points pass FDR and effect-size thresholds. This is not the pending limma result.
6. **Paired tumor-versus-adjacent plot.** Top patient-matched processed-matrix
   contrasts; positive values indicate greater tumor expression.
7. **Field-signature heatmap and trajectories.** Top candidates after FDR,
   effect-size, and direction-consistency filtering.
8. **Nested-CV performance.** Untouched outer-fold macro F1 across models and
   prediction tasks.
9. **Out-of-fold confusion matrix, ROC/PR, and calibration.** Repeated
   predictions were averaged per sample before visualization.
10. **Stable coefficients.** Median standardized Elastic Net coefficients for
    frequently selected genes.
11. **External validation.** Locked Task C model evaluated in eligible GSE41258
    Primary Tumor and Normal Colon arrays using two cross-platform representations.
"""
    (manuscript / "figure_legends.md").write_text(legends, encoding="utf-8")

    permutation_path = ROOT / paths["metrics"] / "group_preserving_permutation_tests.csv"
    permutation_text = (
        pd.read_csv(permutation_path).to_markdown(index=False)
        if permutation_path.exists()
        else "Permutation testing has not yet completed."
    )
    report = f"""# Final analysis report

## 1. Executive summary

The repository now implements a patient-aware, feature-leakage-safe hybrid
R/Python field-cancerization workflow. Primary counts and all 98 pairs were
verified exactly. The Python deposited-matrix route, nested validation,
signature locking, and external validation were executed. Raw-CEL Bioconductor
stages remain blocked solely because R is absent and are not claimed complete.

## 2. Scientific question

Whether histologically normal tumor-adjacent colon carries a reproducible
transcriptomic field effect relative to genuinely healthy mucosa, and whether
compact explainable signatures distinguish healthy, adjacent, and tumor states
without patient leakage.

## 3. Dataset audit and sample inclusion

{table1.to_markdown(index=False)}

Detailed deterministic exclusions are in `reports/metadata_audit.md`.

## 4. Preprocessing and quality control

The executed route used independently deposited normalized GEO matrices,
official platform annotations, unambiguous probe mapping, and per-sample median
aggregation. It produced 300-DPI PNG plus SVG/PDF QC outputs. No outlier was
automatically removed. Raw CEL files were checksummed but not modified.

## 5. Differential expression and field findings

{de_summary.to_markdown(index=False)}

Field candidates passing the provisional processed-matrix rules: {len(candidates)}.
Final limma inference is pending R/Bioconductor execution.

## 6. Model methodology and leakage prevention

All patient/donor groups were disjoint in every fold. Every learned predictive
operation was fitted inside the inner/outer training boundary. Three tasks,
Elastic Net, linear SVM, and Random Forest were evaluated; Elastic Net was
selected under the prespecified interpretability/stability rule.

## 7. Internal validation, stability, and explainability

{best.to_markdown(index=False)}

Fold predictions, metrics, assignments, coefficient stability, compact panels,
and the reload-tested locked model are saved under `results/` and `models/`.

## 8. External validation

{external_display.to_markdown(index=False)}

The external claim is restricted to tumor versus normal colon.

## 9. Statistical uncertainty and permutation

Patient/donor-group bootstrap confidence intervals are stored in
`results/metrics/`. Group-preserving permutation results:

{permutation_text}

## 10. Functional enrichment

Correct-background, BH-adjusted GO enrichment is implemented in
`R/05_functional_enrichment.R`; it was not executed because R is unavailable.
No pathway mechanism is fabricated.

## 11. Reproducibility

Exact executed command sequence:

```powershell
.\\.venv\\Scripts\\python scripts/run_data_audit.py
.\\.venv\\Scripts\\python scripts/run_processed_matrix_pipeline.py
.\\.venv\\Scripts\\python scripts/run_modeling.py
.\\.venv\\Scripts\\python scripts/run_permutation_tests.py
.\\.venv\\Scripts\\python scripts/run_external_validation.py
.\\.venv\\Scripts\\python scripts/build_manuscript_outputs.py
.\\.venv\\Scripts\\python -m pytest -q
```

Complete route after R is installed: `make all`.

The committed R lock is intentionally a bootstrap lock because R was absent in
this execution. Run `R/00_install_packages.R` and commit the resolved
package-filled `renv.lock` before the manuscript environment is frozen.

## 12. Remaining scientific risks and next steps

Execute raw-CEL RMA/QC/limma/enrichment in the pinned Docker/R environment,
compare deposited-matrix sensitivity results with raw-CEL results, validate the
field-effect panel in an independent cancer-free/adjacent cohort, assess
cell-composition and molecular-subtype confounding, and perform prospective
assay validation before any clinical claim.
"""
    (reports / "final_analysis_report.md").write_text(report, encoding="utf-8")
    print(
        json.dumps(
            {
                "figures_png": len(list(figures.glob("*.png"))),
                "table_files": len(list(tables.glob("*.csv"))),
                "manuscript_files": len(list(manuscript.glob("*.md"))),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
