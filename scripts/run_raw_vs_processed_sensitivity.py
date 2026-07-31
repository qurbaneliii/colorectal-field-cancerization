from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.biology import classify_trajectories, field_candidates
from src.visualization.publication_figures import PALETTE, save_figure


COMPARISONS = {
    "adjacent_vs_healthy": "adjacent_normal_vs_healthy",
    "tumor_vs_healthy": "tumor_vs_healthy",
    "tumor_vs_adjacent_paired": "tumor_vs_adjacent_normal_paired",
}


def main() -> None:
    os.chdir(ROOT)
    raw_expression = pd.read_parquet(
        "data/processed/GSE44076_gene_expression_raw_cel_rma.parquet"
    ).set_index("gene_symbol")
    metadata = pd.read_csv("data/metadata/gse44076_samples.csv", dtype={"patient_id": str})
    raw_de = pd.read_csv("results/tables/supplementary_full_raw_cel_de_results.csv")
    processed_de = pd.read_csv(
        "results/differential_expression/processed_matrix_exploratory_all_comparisons.csv"
    )
    rows = []
    changed_frames = []
    plot_frames = []
    for raw_name, processed_name in COMPARISONS.items():
        raw = raw_de[raw_de["comparison"].eq(raw_name)][
            ["gene_symbol", "log2_fold_change", "adjusted_p_value"]
        ].rename(
            columns={
                "log2_fold_change": "raw_log2fc",
                "adjusted_p_value": "raw_fdr",
            }
        )
        processed = processed_de[processed_de["comparison"].eq(processed_name)][
            ["gene_symbol", "log2_fold_change", "adjusted_p_value"]
        ].rename(
            columns={
                "log2_fold_change": "deposited_log2fc",
                "adjusted_p_value": "deposited_fdr",
            }
        )
        merged = raw.merge(processed, on="gene_symbol", how="inner")
        merged["comparison"] = raw_name
        merged["raw_significant"] = merged["raw_fdr"].lt(0.05) & merged["raw_log2fc"].abs().ge(0.5)
        merged["deposited_significant"] = merged["deposited_fdr"].lt(0.05) & merged[
            "deposited_log2fc"
        ].abs().ge(0.5)
        merged["sign_agreement"] = np.sign(merged["raw_log2fc"]) == np.sign(
            merged["deposited_log2fc"]
        )
        merged["conclusion_changed_materially"] = (
            merged["raw_significant"] != merged["deposited_significant"]
        ) | ~merged["sign_agreement"]
        top_raw = set(merged.nsmallest(100, "raw_fdr")["gene_symbol"])
        top_processed = set(merged.nsmallest(100, "deposited_fdr")["gene_symbol"])
        rows.append(
            {
                "comparison": raw_name,
                "shared_tested_genes": len(merged),
                "pearson_log2fc_correlation": stats.pearsonr(
                    merged["raw_log2fc"], merged["deposited_log2fc"]
                ).statistic,
                "spearman_rank_correlation": stats.spearmanr(
                    merged["raw_log2fc"], merged["deposited_log2fc"]
                ).statistic,
                "sign_agreement_fraction": merged["sign_agreement"].mean(),
                "top_100_gene_overlap": len(top_raw & top_processed),
                "raw_significant_genes": int(merged["raw_significant"].sum()),
                "deposited_significant_genes": int(merged["deposited_significant"].sum()),
                "significant_gene_overlap": int(
                    (merged["raw_significant"] & merged["deposited_significant"]).sum()
                ),
                "materially_changed_genes": int(merged["conclusion_changed_materially"].sum()),
                "pathway_agreement": "reported after raw enrichment; no mechanism inferred",
            }
        )
        changed_frames.append(merged[merged["conclusion_changed_materially"]])
        plot_frames.append(merged)
    concordance = pd.DataFrame(rows)
    concordance.to_csv("results/tables/raw_vs_deposited_de_concordance.csv", index=False)
    pd.concat(changed_frames, ignore_index=True).to_csv(
        "results/tables/raw_vs_deposited_material_changes.csv", index=False
    )
    plot = pd.concat(plot_frames, ignore_index=True)
    figure, axes = plt.subplots(1, 3, figsize=(16, 5), sharex=False, sharey=False)
    for axis, (comparison, frame) in zip(axes, plot.groupby("comparison"), strict=True):
        axis.scatter(frame["deposited_log2fc"], frame["raw_log2fc"], s=6, alpha=0.25)
        limits = [
            min(frame["deposited_log2fc"].min(), frame["raw_log2fc"].min()),
            max(frame["deposited_log2fc"].max(), frame["raw_log2fc"].max()),
        ]
        axis.plot(limits, limits, linestyle="--", color="black", linewidth=0.8)
        axis.set(xlabel="Deposited-matrix log2 difference", ylabel="Raw-CEL RMA limma log2FC")
        axis.set_title(f"{comparison}\nshared genes n={len(frame):,}")
    figure.suptitle("Raw-CEL versus deposited-matrix differential-expression concordance")
    save_figure(figure, Path("results/figures/raw_vs_deposited_logfc_concordance"))

    trajectories = classify_trajectories(raw_expression, metadata, tolerance=0.15)
    trajectories["analysis_provenance"] = "raw_cel_rma"
    trajectories.to_csv("results/tables/final_expression_trajectory_genes.csv", index=False)
    adjacent_de = raw_de[raw_de["comparison"].eq("adjacent_vs_healthy")].copy()
    candidates = field_candidates(
        raw_expression,
        metadata,
        adjacent_de,
        trajectories,
        fdr=0.05,
        effect=0.5,
    )
    processed_adjacent = processed_de[
        processed_de["comparison"].eq("adjacent_normal_vs_healthy")
    ].set_index("gene_symbol")
    candidates["deposited_log2fc"] = candidates["gene_symbol"].map(
        processed_adjacent["log2_fold_change"]
    )
    candidates["deposited_adjusted_p_value"] = candidates["gene_symbol"].map(
        processed_adjacent["adjusted_p_value"]
    )
    candidates["deposited_sign_agreement"] = np.sign(candidates["log2_fold_change"]) == np.sign(
        candidates["deposited_log2fc"]
    )
    candidates["outlier_sensitivity_status"] = "no_qc_exclusions_applied"
    candidates["analysis_provenance"] = "raw_cel_rma_limma_primary"
    candidates["candidate_rule"] = (
        "raw FDR<0.05; |log2FC|>=0.5; direction consistency>=0.70; deposited sign agreement"
    )
    candidates = candidates[candidates["deposited_sign_agreement"]].copy()
    candidates.to_csv("results/tables/final_field_cancerization_candidates.csv", index=False)

    genes = candidates["gene_symbol"].head(30).tolist()
    if genes:
        ordered = metadata.sort_values(["tissue_class", "patient_id"])["geo_accession"].tolist()
        z = raw_expression.loc[genes, ordered]
        z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0)
        fig, ax = plt.subplots(figsize=(12, 8))
        sns.heatmap(z, cmap="vlag", center=0, xticklabels=False, yticklabels=True, ax=ax)
        ax.set_title(f"Raw-CEL field-cancerization candidates (n={len(genes)} genes)")
        ax.set_xlabel(f"GSE44076 samples ordered by tissue class (n={len(ordered)})")
        ax.set_ylabel("Gene symbol")
        save_figure(fig, Path("results/figures/final_field_signature_heatmap"))
        top = genes[:6]
        long = (
            raw_expression.loc[top]
            .T.reset_index(names="geo_accession")
            .merge(metadata[["geo_accession", "tissue_class"]], on="geo_accession")
            .melt(
                id_vars=["geo_accession", "tissue_class"],
                var_name="gene_symbol",
                value_name="expression",
            )
        )
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
        for axis, gene in zip(axes.ravel(), top, strict=False):
            subset = long[long["gene_symbol"].eq(gene)]
            sns.boxplot(
                data=subset,
                x="tissue_class",
                y="expression",
                order=["healthy", "adjacent_normal", "tumor"],
                palette=PALETTE,
                ax=axis,
            )
            axis.set_title(gene)
            axis.tick_params(axis="x", rotation=20)
            axis.set_ylabel("Raw-CEL RMA log2 expression")
        fig.suptitle("Healthy → adjacent-normal → tumor expression distributions")
        save_figure(fig, Path("results/figures/final_top_gene_trajectory_plots"))

    processed_de.to_csv(
        "results/tables/supplementary_full_processed_matrix_de_results.csv", index=False
    )
    report = f"""# Raw-CEL versus processed-matrix sensitivity report

Raw-CEL RMA limma is the primary biological analysis. GEO deposited-matrix
Welch and paired-t results are retained only as a sensitivity analysis.

{concordance.to_markdown(index=False)}

Field-cancerization candidates additionally required sample-direction
consistency and deposited-matrix sign agreement. Concordance does not establish
causality, and conclusions that materially change are listed in
`results/tables/raw_vs_deposited_material_changes.csv`.
"""
    Path("reports/raw_vs_processed_sensitivity_report.md").write_text(
        report, encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "comparisons": len(concordance),
                "field_candidates": len(candidates),
                "trajectory_genes": len(trajectories),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
