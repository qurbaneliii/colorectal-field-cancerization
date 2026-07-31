from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


PALETTE = {
    "healthy": "#2A9D8F",
    "adjacent_normal": "#E9C46A",
    "tumor": "#E76F51",
    "normal_colon": "#2A9D8F",
    "primary_tumor": "#E76F51",
}


def save_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def expression_boxplot(expression: pd.DataFrame, stem: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.boxplot(
        expression.to_numpy(),
        showfliers=False,
        patch_artist=True,
        boxprops={"facecolor": "#A8DADC", "linewidth": 0.3},
        medianprops={"color": "#1D3557", "linewidth": 0.5},
    )
    ax.set_title(title)
    ax.set_xlabel(f"Samples (n={expression.shape[1]})")
    ax.set_ylabel("Deposited normalized log2 expression")
    ax.set_xticks([])
    save_figure(fig, stem)


def pca_plot(
    expression: pd.DataFrame, metadata: pd.DataFrame, stem: Path, title: str
) -> pd.DataFrame:
    ordered = metadata.set_index("geo_accession").loc[expression.columns]
    variances = expression.var(axis=1).sort_values(ascending=False)
    use = variances.head(min(5000, len(variances))).index
    x = expression.loc[use].T
    model = PCA(n_components=2, random_state=44076).fit(x)
    xy = model.transform(x)
    pca = pd.DataFrame(
        {
            "geo_accession": expression.columns,
            "PC1": xy[:, 0],
            "PC2": xy[:, 1],
            "tissue_class": ordered["tissue_class"].to_numpy(),
            "patient_id": ordered["patient_id"].astype(str).to_numpy(),
        }
    )
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for label, group in pca.groupby("tissue_class"):
        ax.scatter(
            group["PC1"],
            group["PC2"],
            s=32,
            alpha=0.78,
            label=f"{label} (n={len(group)})",
            color=PALETTE.get(label, "#457B9D"),
            edgecolor="white",
            linewidth=0.3,
        )
    ax.set_xlabel(f"PC1 ({model.explained_variance_ratio_[0] * 100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({model.explained_variance_ratio_[1] * 100:.1f}% variance)")
    ax.set_title(title)
    ax.legend(frameon=False)
    save_figure(fig, stem)
    return pca


def correlation_heatmap(expression: pd.DataFrame, stem: Path, title: str) -> None:
    top_genes = expression.var(axis=1).nlargest(min(5000, len(expression))).index
    top = expression.loc[top_genes]
    corr = top.corr(method="pearson")
    fig, ax = plt.subplots(figsize=(9, 8))
    sns.heatmap(corr, cmap="vlag", center=0.8, xticklabels=False, yticklabels=False, ax=ax)
    ax.set_title(title)
    save_figure(fig, stem)


def top_variable_heatmap(
    expression: pd.DataFrame, metadata: pd.DataFrame, stem: Path, n_genes: int = 50
) -> None:
    top = expression.var(axis=1).nlargest(n_genes).index
    z = expression.loc[top]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, 1), axis=0)
    ordered_samples = (
        metadata.sort_values(["tissue_class", "patient_id"])["geo_accession"].tolist()
    )
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(
        z[ordered_samples],
        cmap="vlag",
        center=0,
        xticklabels=False,
        yticklabels=True,
        cbar_kws={"label": "Gene-wise z-score"},
        ax=ax,
    )
    ax.set_title(f"Top {n_genes} variable genes across GSE44076")
    ax.set_xlabel(f"Samples (n={len(ordered_samples)}; ordered by tissue class)")
    ax.set_ylabel("Gene symbol")
    save_figure(fig, stem)


def volcano_plot(de: pd.DataFrame, stem: Path, title: str, fdr: float, effect: float) -> None:
    plot = de.copy()
    plot["minus_log10_fdr"] = -np.log10(plot["adjusted_p_value"].clip(lower=1e-300))
    significant = plot["adjusted_p_value"].lt(fdr) & plot["log2_fold_change"].abs().ge(effect)
    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.scatter(
        plot.loc[~significant, "log2_fold_change"],
        plot.loc[~significant, "minus_log10_fdr"],
        s=8,
        alpha=0.35,
        color="#8D99AE",
        label="Other",
    )
    ax.scatter(
        plot.loc[significant, "log2_fold_change"],
        plot.loc[significant, "minus_log10_fdr"],
        s=11,
        alpha=0.7,
        color="#D62828",
        label=f"FDR < {fdr}, |log2FC| ≥ {effect}",
    )
    ax.axvline(effect, color="black", linewidth=0.7, linestyle="--")
    ax.axvline(-effect, color="black", linewidth=0.7, linestyle="--")
    ax.axhline(-np.log10(fdr), color="black", linewidth=0.7, linestyle="--")
    ax.set_xlabel("Log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(title)
    ax.legend(frameon=False)
    save_figure(fig, stem)
