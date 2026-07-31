from __future__ import annotations

import numpy as np
import pandas as pd


def feature_stability(
    coefficients: pd.DataFrame, total_outer_folds: int | None = None
) -> pd.DataFrame:
    if coefficients.empty:
        return pd.DataFrame()
    selected = coefficients[coefficients["selected_nonzero"]].copy()
    selected["fold_key"] = (
        selected["repeat"].astype(str) + ":" + selected["outer_fold"].astype(str)
    )
    if total_outer_folds is None:
        total_outer_folds = coefficients[["repeat", "outer_fold"]].drop_duplicates().shape[0]
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
    return summary.sort_values(
        ["task", "selection_frequency", "median_absolute_coefficient"],
        ascending=[True, False, False],
    )
