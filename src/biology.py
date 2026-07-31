from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _welch_de(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    statistic, p_value = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="raise")
    se = np.sqrt(a.var(axis=1, ddof=1) / a.shape[1] + b.var(axis=1, ddof=1) / b.shape[1])
    return statistic, p_value, se


def _paired_de(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    statistic, p_value = stats.ttest_rel(a, b, axis=1, nan_policy="raise")
    delta = a - b
    se = delta.std(axis=1, ddof=1) / np.sqrt(delta.shape[1])
    return statistic, p_value, se


def exploratory_differential_expression(
    expression: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    """Processed-matrix sensitivity analysis; does not claim to replace limma."""
    meta = metadata.set_index("geo_accession").loc[expression.columns]
    comparisons = []
    for first, second, label in [
        ("adjacent_normal", "healthy", "adjacent_normal_vs_healthy"),
        ("tumor", "healthy", "tumor_vs_healthy"),
    ]:
        a = expression.loc[:, meta["tissue_class"].eq(first)].to_numpy()
        b = expression.loc[:, meta["tissue_class"].eq(second)].to_numpy()
        statistic, p_value, se = _welch_de(a, b)
        comparisons.append(
            _format_de(expression.index, a.mean(axis=1) - b.mean(axis=1), se, statistic, p_value, label)
        )

    pair_meta = meta[meta["tissue_class"].isin(["adjacent_normal", "tumor"])].copy()
    adjacent = pair_meta[pair_meta["tissue_class"].eq("adjacent_normal")]
    tumor = pair_meta[pair_meta["tissue_class"].eq("tumor")]
    common = sorted(set(adjacent["patient_id"]).intersection(tumor["patient_id"]))
    adjacent_samples = [
        adjacent.index[adjacent["patient_id"].eq(patient)][0] for patient in common
    ]
    tumor_samples = [tumor.index[tumor["patient_id"].eq(patient)][0] for patient in common]
    a = expression[tumor_samples].to_numpy()
    b = expression[adjacent_samples].to_numpy()
    statistic, p_value, se = _paired_de(a, b)
    comparisons.append(
        _format_de(
            expression.index,
            a.mean(axis=1) - b.mean(axis=1),
            se,
            statistic,
            p_value,
            "tumor_vs_adjacent_normal_paired",
        )
    )
    return pd.concat(comparisons, ignore_index=True)


def _format_de(
    genes: pd.Index,
    effect: np.ndarray,
    se: np.ndarray,
    statistic: np.ndarray,
    p_value: np.ndarray,
    label: str,
) -> pd.DataFrame:
    adjusted = multipletests(np.nan_to_num(p_value, nan=1.0), method="fdr_bh")[1]
    return pd.DataFrame(
        {
            "gene_symbol": genes,
            "entrez_id": "",
            "log2_fold_change": effect,
            "standard_error": se,
            "test_statistic": statistic,
            "p_value": p_value,
            "adjusted_p_value": adjusted,
            "comparison": label,
            "direction": np.where(effect >= 0, "up", "down"),
            "annotation_status": "unique_gene_symbol",
            "method": "processed_matrix_welch_t"
            if "paired" not in label
            else "processed_matrix_paired_t",
        }
    )


def classify_trajectories(
    expression: pd.DataFrame, metadata: pd.DataFrame, tolerance: float = 0.15
) -> pd.DataFrame:
    meta = metadata.set_index("geo_accession").loc[expression.columns]
    means = {
        label: expression.loc[:, meta["tissue_class"].eq(label)].mean(axis=1)
        for label in ("healthy", "adjacent_normal", "tumor")
    }
    result = pd.DataFrame(means)
    ha = result["adjacent_normal"] - result["healthy"]
    at = result["tumor"] - result["adjacent_normal"]
    strong = 0.5
    category = np.full(len(result), "non_monotonic", dtype=object)
    category[(ha >= tolerance) & (at >= tolerance)] = "monotonic_up"
    category[(ha <= -tolerance) & (at <= -tolerance)] = "monotonic_down"
    category[(np.abs(ha) >= strong) & (np.abs(at) < tolerance)] = "field_shift_then_plateau"
    category[(np.abs(ha) < tolerance) & (np.abs(at) >= strong)] = "tumor_specific"
    category[
        (np.abs(ha) >= strong)
        & (np.abs(result["adjacent_normal"] - result["tumor"]) >= strong)
        & ((result["adjacent_normal"] - result["healthy"])
           * (result["adjacent_normal"] - result["tumor"]) > 0)
    ] = "adjacent_specific"
    result["healthy_to_adjacent"] = ha
    result["adjacent_to_tumor"] = at
    result["trajectory_category"] = category
    result["rule_tolerance"] = tolerance
    result.index.name = "gene_symbol"
    return result.reset_index()


def field_candidates(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    de: pd.DataFrame,
    trajectories: pd.DataFrame,
    fdr: float,
    effect: float,
) -> pd.DataFrame:
    meta = metadata.set_index("geo_accession").loc[expression.columns]
    adjacent = expression.loc[:, meta["tissue_class"].eq("adjacent_normal")]
    healthy_median = expression.loc[:, meta["tissue_class"].eq("healthy")].median(axis=1)
    adjacent_de = de[de["comparison"].eq("adjacent_normal_vs_healthy")].copy()
    adjacent_de = adjacent_de[
        adjacent_de["adjusted_p_value"].lt(fdr)
        & adjacent_de["log2_fold_change"].abs().ge(effect)
    ]
    direction = adjacent_de.set_index("gene_symbol")["direction"]
    consistency = pd.Series(index=adjacent_de["gene_symbol"], dtype=float)
    for gene in consistency.index:
        if direction[gene] == "up":
            consistency[gene] = (adjacent.loc[gene] > healthy_median[gene]).mean()
        else:
            consistency[gene] = (adjacent.loc[gene] < healthy_median[gene]).mean()
    adjacent_de["adjacent_direction_consistency"] = adjacent_de["gene_symbol"].map(consistency)
    adjacent_de = adjacent_de[adjacent_de["adjacent_direction_consistency"].ge(0.70)]
    return adjacent_de.merge(trajectories, on="gene_symbol", how="left").sort_values(
        ["adjusted_p_value", "log2_fold_change"], ascending=[True, False]
    )
