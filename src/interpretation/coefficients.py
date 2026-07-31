from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from src.visualization.publication_figures import save_figure


def coefficient_plot(stability: pd.DataFrame, task: str, stem: Path, top_n: int = 30) -> None:
    subset = stability[stability["task"].eq(task)].copy()
    if subset.empty:
        return
    subset = subset.sort_values("median_absolute_coefficient").tail(top_n)
    fig, ax = plt.subplots(figsize=(8, max(5, 0.24 * len(subset))))
    colors = ["#D62828" if value > 0 else "#277DA1" for value in subset["median_coefficient"]]
    ax.barh(subset["gene_symbol"], subset["median_coefficient"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Median standardized Elastic Net coefficient")
    ax.set_ylabel("Gene symbol")
    ax.set_title(f"{task}: stable feature coefficients")
    save_figure(fig, stem)
