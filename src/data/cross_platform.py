from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def common_gene_symbols(primary: pd.DataFrame, external: pd.DataFrame) -> list[str]:
    common = sorted(set(primary["gene_symbol"]).intersection(external["gene_symbol"]))
    if not common:
        raise ValueError("No common gene symbols between platforms")
    return common


def within_sample_percentile_rank(matrix: np.ndarray) -> np.ndarray:
    """Rank genes within each sample without using any labels or other samples."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("Expected samples x genes")
    n_features = matrix.shape[1]
    if n_features < 2:
        return np.ones_like(matrix, dtype=float)
    return np.vstack(
        [(rankdata(row, method="average") - 1.0) / (n_features - 1.0) for row in matrix]
    )
