from __future__ import annotations

import json
import os
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.visualization.publication_figures import PALETTE, save_figure

TABLES = ROOT / "results/tables"
METRICS = ROOT / "results/metrics"
FIGURES = ROOT / "results/figures"


def study_design() -> None:
    fig, ax = plt.subplots(figsize=(12, 6.5))
    ax.axis("off")
    boxes = [
        (0.03, 0.58, 0.25, 0.27, "GSE44076 development\n246 arrays\n50 healthy | 98 adjacent\n98 tumor"),
        (0.37, 0.58, 0.25, 0.27, "Primary Task B\nhealthy vs adjacent-normal\ngrouped nested CV\ninternal validation only"),
        (0.72, 0.58, 0.25, 0.27, "Secondary Task C\ntumor vs adjacent-normal\nlocked on GSE44076"),
        (0.20, 0.10, 0.27, 0.25, "Biological inference\nU0/U1/U2 limma\nfield evidence tiers\ncomposition + stress sensitivity"),
        (0.59, 0.10, 0.30, 0.25, "GSE41258 external transfer\n233 canonical patient-tissue arrays\nnormal_colon vs primary_tumor\n2-gene rank transport refit"),
    ]
    for x, y, width, height, label in boxes:
        ax.add_patch(
            plt.Rectangle(
                (x, y), width, height, transform=ax.transAxes, facecolor="#F4F7FA",
                edgecolor="#264653", linewidth=1.5
            )
        )
        ax.text(x + width / 2, y + height / 2, label, transform=ax.transAxes,
                ha="center", va="center", fontsize=9.2)
    arrows = [((0.28, 0.715), (0.37, 0.715)), ((0.62, 0.715), (0.72, 0.715)),
              ((0.18, 0.58), (0.30, 0.35)), ((0.845, 0.58), (0.74, 0.35))]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, xycoords="axes fraction",
                    arrowprops={"arrowstyle": "->", "lw": 1.6, "color": "#264653"})
    ax.set_title(
        "Study objectives and dataset roles: field-effect discovery is distinct from tumor-normal transfer",
        fontsize=14,
    )
    save_figure(fig, FIGURES / "study_objectives_and_dataset_roles")


def sample_inclusion_flow() -> None:
    primary = pd.read_csv(ROOT / "data/metadata/gse44076_samples.csv")
    external = pd.read_csv(ROOT / "data/metadata/gse41258_samples.csv")
    included = external[external["inclusion_status"].eq("included")]
    fig, ax = plt.subplots(figsize=(12, 6.5))
    ax.axis("off")
    boxes = [
        (0.08, 0.63, 0.34, 0.17, f"GSE44076 deposited arrays\nn={len(primary)}", "#E8F1F2"),
        (0.58, 0.63, 0.34, 0.17, f"GSE41258 deposited arrays\nn={len(external)}", "#E8F1F2"),
        (
            0.05,
            0.18,
            0.40,
            0.25,
            "Included: 50 healthy | 98 adjacent-normal | 98 tumor\n"
            "98 complete tumor-adjacent pairs",
            "#EFF8EA",
        ),
        (
            0.55,
            0.18,
            0.40,
            0.25,
            f"Canonical external: n={len(included)} arrays\n"
            "52 normal colon | 181 primary tumor\n"
            f"{included['donor_or_patient_group'].nunique()} unique patients; "
            f"excluded={len(external) - len(included)}",
            "#EFF8EA",
        ),
    ]
    for x, y, width, height, label, color in boxes:
        ax.add_patch(
            plt.Rectangle(
                (x, y),
                width,
                height,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor="#264653",
                linewidth=1.5,
            )
        )
        ax.text(
            x + width / 2,
            y + height / 2,
            label,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )
    for x in (0.25, 0.75):
        ax.annotate(
            "",
            xy=(x, 0.43),
            xytext=(x, 0.63),
            xycoords="axes fraction",
            arrowprops={"arrowstyle": "->", "lw": 1.6},
        )
    ax.set_title("Sample inclusion flow and analysis-ready cohorts", fontsize=14)
    save_figure(fig, FIGURES / "sample_count_flowchart")


def raw_qc_summary() -> None:
    qc = pd.read_csv(TABLES / "table_2_raw_cel_qc_summary.csv")
    plot = qc.melt(
        id_vars="accession",
        value_vars=[
            "aqm_flagged_any",
            "aqm_two_or_more_flags",
            "one_custom_failure_metric",
            "two_or_more_custom_failure_metrics",
            "automatic_exclusions",
        ],
        var_name="criterion",
        value_name="array_count",
    )
    labels = {
        "aqm_flagged_any": "AQM >=1 flag",
        "aqm_two_or_more_flags": "AQM >=2 flags",
        "one_custom_failure_metric": "Custom exactly 1",
        "two_or_more_custom_failure_metrics": "Custom >=2",
        "automatic_exclusions": "Excluded",
    }
    plot["criterion"] = plot["criterion"].map(labels)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    sns.barplot(data=plot, x="criterion", y="array_count", hue="accession", ax=ax)
    ax.set(
        xlabel="Review/exclusion criterion",
        ylabel="Arrays",
        title="Raw-CEL QC review flags under the conservative two-failure exclusion rule",
    )
    ax.tick_params(axis="x", rotation=25)
    ax.legend(title="Dataset", frameon=False)
    save_figure(fig, FIGURES / "raw_array_qc_summary")


def volcano(frame: pd.DataFrame, title: str, stem: str) -> None:
    plotted = frame.copy()
    plotted["minus_log10_fdr"] = -np.log10(plotted["adjusted_p_value"].clip(lower=1e-300))
    significant = plotted["adjusted_p_value"].lt(0.05) & plotted[
        "log2_fold_change"
    ].abs().ge(0.5)
    fig, ax = plt.subplots(figsize=(7.5, 6))
    ax.scatter(
        plotted.loc[~significant, "log2_fold_change"],
        plotted.loc[~significant, "minus_log10_fdr"],
        s=7,
        alpha=0.30,
        color="#8D99AE",
    )
    ax.scatter(
        plotted.loc[significant, "log2_fold_change"],
        plotted.loc[significant, "minus_log10_fdr"],
        s=9,
        alpha=0.62,
        color="#C44536",
    )
    ax.axvline(-0.5, color="black", linewidth=0.7, linestyle="--")
    ax.axvline(0.5, color="black", linewidth=0.7, linestyle="--")
    ax.axhline(-np.log10(0.05), color="black", linewidth=0.7, linestyle="--")
    ax.set(xlabel="log2 fold change", ylabel="-log10 BH FDR", title=title)
    ax.text(
        0.02,
        0.98,
        f"FDR < 0.05 and |log2FC| >= 0.5: {int(significant.sum()):,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
    )
    save_figure(fig, FIGURES / stem)


def de_figures() -> None:
    de = pd.read_csv(TABLES / "supplementary_full_raw_cel_de_results.csv")
    unadjusted = de[
        de["model"].eq("U0_unadjusted")
        & de["comparison"].eq("adjacent_vs_healthy")
    ]
    adjusted = de[
        de["model"].eq("U1_age_sex_adjusted")
        & de["comparison"].eq("adjacent_vs_healthy")
    ]
    paired = de[
        de["model"].eq("P_patient_fixed_effect")
        & de["comparison"].eq("tumor_vs_adjacent")
    ]
    volcano(
        unadjusted,
        "Unadjusted adjacent-normal versus healthy mucosa",
        "unadjusted_adjacent_vs_healthy_volcano",
    )
    volcano(
        adjusted,
        "Age/sex-adjusted adjacent-normal versus healthy mucosa",
        "adjusted_adjacent_vs_healthy_volcano",
    )
    volcano(
        paired,
        "Paired patient-fixed-effect tumor versus adjacent-normal mucosa",
        "paired_tumor_vs_adjacent_volcano",
    )


def top_gene_distributions() -> None:
    high = pd.read_csv(TABLES / "final_high_confidence_field_signature.csv")
    genes = high.sort_values("u1_fdr").head(6)["gene_symbol"].tolist()
    expression = pd.read_parquet(
        ROOT / "data/processed/GSE44076_gene_expression_raw_cel_rma.parquet"
    ).set_index("gene_symbol")
    metadata = pd.read_csv(ROOT / "data/metadata/gse44076_samples.csv").set_index(
        "geo_accession"
    )
    long = (
        expression.loc[genes]
        .T.rename_axis("geo_accession")
        .reset_index()
        .melt(id_vars="geo_accession", var_name="gene_symbol", value_name="expression")
    )
    long["tissue_class"] = long["geo_accession"].map(metadata["tissue_class"])
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), sharey=False)
    order = ["healthy", "adjacent_normal", "tumor"]
    for ax, gene in zip(axes.flat, genes, strict=True):
        frame = long[long["gene_symbol"].eq(gene)]
        sns.boxplot(
            data=frame,
            x="tissue_class",
            y="expression",
            order=order,
            palette=PALETTE,
            hue="tissue_class",
            legend=False,
            showfliers=False,
            ax=ax,
        )
        sns.stripplot(
            data=frame,
            x="tissue_class",
            y="expression",
            order=order,
            color="black",
            alpha=0.25,
            size=2,
            ax=ax,
        )
        ax.set_title(gene)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=25)
    fig.suptitle("Top adjusted high-confidence field genes across all three tissue groups")
    fig.tight_layout()
    save_figure(fig, FIGURES / "top_field_gene_group_distributions")


def tier_summary() -> None:
    tiers = pd.read_csv(TABLES / "field_gene_evidence_tiers.csv")
    counts = tiers["evidence_tier"].value_counts().rename_axis("tier").reset_index(name="genes")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    sns.barplot(data=counts, x="genes", y="tier", color="#457B9D", ax=ax)
    ax.set(
        xlabel="Genes",
        ylabel="Evidence tier",
        title="Field-gene evidence tiers after covariate, subgroup, QC, and deposited-data checks",
    )
    for container in ax.containers:
        ax.bar_label(container, fmt="%d")
    save_figure(fig, FIGURES / "field_gene_evidence_tier_summary")


def enrichment_figure() -> None:
    frame = pd.read_csv(TABLES / "high_confidence_field_enrichment.csv").copy()
    frame = frame.nsmallest(15, "adjusted_p_value").sort_values("adjusted_p_value")
    frame["label"] = frame["Description"].map(lambda value: textwrap.shorten(value, 58))
    frame["display_label"] = (
        frame["analysis_set"].str.replace("_", " ") + " | " + frame["label"]
    )
    frame["minus_log10_fdr"] = -np.log10(frame["adjusted_p_value"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.barplot(
        data=frame,
        x="minus_log10_fdr",
        y="display_label",
        hue="analysis_set",
        dodge=False,
        ax=ax,
    )
    ax.set(
        xlabel="-log10 BH FDR",
        ylabel="Analysis set and nonredundant term",
        title="High-confidence and composition-stratified field-gene enrichment",
    )
    ax.legend(title="Analysis set", frameon=False, fontsize=8)
    save_figure(fig, FIGURES / "high_confidence_field_enrichment")


def external_sensitivity() -> None:
    metrics = pd.read_csv(METRICS / "external_patient_tissue_metrics.csv")
    frame = metrics[
        metrics["representation"].eq("within_sample_percentile_rank")
        & metrics["threshold_policy"].eq("primary_gse44076_locked")
    ].copy()
    long = frame.melt(
        id_vars="evaluation_set",
        value_vars=["roc_auc", "f1_macro", "balanced_accuracy", "brier_score"],
        var_name="metric",
        value_name="value",
    )
    fig, ax = plt.subplots(figsize=(10, 5.8))
    sns.barplot(data=long, x="metric", y="value", hue="evaluation_set", ax=ax)
    ax.set(
        xlabel="Metric",
        ylabel="Point estimate",
        title="External patient-structure sensitivity at the GSE44076-locked threshold",
        ylim=(0, 1.05),
    )
    ax.legend(title="Evaluation set", frameon=False, fontsize=8)
    save_figure(fig, FIGURES / "external_patient_structure_sensitivity")


def signature_intersection() -> None:
    audit = pd.read_csv(TABLES / "external_common_signature_audit.csv")
    matrix = audit.set_index("gene_symbol")[["present_gse44076", "present_gse41258"]]
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for row, gene in enumerate(matrix.index):
        for column, dataset in enumerate(matrix.columns):
            present = bool(matrix.loc[gene, dataset])
            ax.scatter(
                column,
                row,
                s=450,
                marker="o" if present else "x",
                color="#2A9D8F" if present else "#C44536",
                linewidth=2,
            )
    ax.set_xticks([0, 1], ["GSE44076 / GPL13667", "GSE41258 / GPL96"])
    ax.set_yticks(range(len(matrix)), matrix.index)
    ax.invert_yaxis()
    ax.set_title("Locked Task C signature and label-independent cross-platform intersection")
    ax.text(0.5, -0.20, "CEMIP and ETV4 form the transport subset; FOXQ1 is absent on GPL96",
            transform=ax.transAxes, ha="center")
    save_figure(fig, FIGURES / "signature_intersection")


def main() -> None:
    os.chdir(ROOT)
    study_design()
    sample_inclusion_flow()
    raw_qc_summary()
    de_figures()
    top_gene_distributions()
    tier_summary()
    enrichment_figure()
    external_sensitivity()
    signature_intersection()
    generated = [
        "study_objectives_and_dataset_roles",
        "sample_count_flowchart",
        "raw_array_qc_summary",
        "unadjusted_adjacent_vs_healthy_volcano",
        "adjusted_adjacent_vs_healthy_volcano",
        "paired_tumor_vs_adjacent_volcano",
        "top_field_gene_group_distributions",
        "field_gene_evidence_tier_summary",
        "high_confidence_field_enrichment",
        "external_patient_structure_sensitivity",
        "signature_intersection",
    ]
    for stem in generated:
        for suffix in (".png", ".pdf", ".svg"):
            path = FIGURES / f"{stem}{suffix}"
            if not path.exists() or path.stat().st_size == 0:
                raise AssertionError(f"Missing publication figure: {path}")
    print(json.dumps({"figure_stems": generated, "renderings": len(generated) * 3}, indent=2))


if __name__ == "__main__":
    main()
