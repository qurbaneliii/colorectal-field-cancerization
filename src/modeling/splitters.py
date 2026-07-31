from __future__ import annotations

import numpy as np
from sklearn.model_selection import StratifiedGroupKFold

from src.data.validation import assert_no_group_overlap


def stratified_group_splits(
    y: np.ndarray, groups: np.ndarray, n_splits: int, seed: int
) -> list[tuple[np.ndarray, np.ndarray]]:
    splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    dummy = np.zeros((len(y), 1))
    splits = list(splitter.split(dummy, y, groups))
    for train, test in splits:
        assert_no_group_overlap(train, test, groups)
    return splits
