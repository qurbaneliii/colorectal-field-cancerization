from __future__ import annotations

import numpy as np
import pandas as pd


def feature_stability(
    coefficients: pd.DataFrame,
    total_outer_folds: int | None = None,
    strict_frequency: float = 0.65,
    strict_sign_consistency: float = 0.80,
    consensus_frequency: float = 0.50,
) -> pd.DataFrame:
    if coefficients.empty:
        return pd.DataFrame()
    selected = coefficients[coefficients["selected_nonzero"]].copy()
    selected["fold_key"] = (
        selected["repeat"].astype(str) + ":" + selected["outer_fold"].astype(str)
    )
    if total_outer_folds is None:
        total_outer_folds = coefficients[["repeat", "outer_fold"]].drop_duplicates().shape[0]
    total_repeats = coefficients["repeat"].nunique()
    summary = (
        selected.groupby(["task", "class", "gene_symbol"])
        .agg(
            selected_folds=("fold_key", "nunique"),
            positive_fraction=("coefficient", lambda x: float((x > 0).mean())),
            median_coefficient=("coefficient", "median"),
            median_absolute_coefficient=("coefficient", lambda x: float(np.median(np.abs(x)))),
            coefficient_standard_deviation=("coefficient", "std"),
            repeat_count=("repeat", "nunique"),
        )
        .reset_index()
    )
    summary["selection_frequency"] = summary["selected_folds"] / total_outer_folds
    summary["sign_consistency"] = np.maximum(
        summary["positive_fraction"], 1 - summary["positive_fraction"]
    )
    summary["total_unique_repeat_fold_count"] = int(total_outer_folds)
    summary["repeat_coverage"] = summary["repeat_count"] / int(total_repeats)
    summary["strictly_stable_gene"] = (
        summary["selection_frequency"].ge(strict_frequency)
        & summary["sign_consistency"].ge(strict_sign_consistency)
    )
    summary["compact_consensus_gene"] = (
        ~summary["strictly_stable_gene"]
        & summary["selection_frequency"].ge(consensus_frequency)
        & summary["sign_consistency"].ge(strict_sign_consistency)
    )
    summary["exploratory_selected_gene"] = ~(
        summary["strictly_stable_gene"] | summary["compact_consensus_gene"]
    )
    summary["stability_class"] = np.select(
        [summary["strictly_stable_gene"], summary["compact_consensus_gene"]],
        ["strictly_stable_gene", "compact_consensus_gene"],
        default="exploratory_selected_gene",
    )
    return summary.sort_values(
        ["task", "selection_frequency", "median_absolute_coefficient"],
        ascending=[True, False, False],
    )
