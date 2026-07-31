from __future__ import annotations

import pandas as pd
import pytest

from src.data.platform import validate_chip_type
from src.modeling.stability import feature_stability


def test_platform_mismatch_raises():
    with pytest.raises(ValueError, match="chip mismatch"):
        validate_chip_type("HG-U133A", "HG-U219")


def test_repeat_aware_stability_denominator():
    coefficients = pd.DataFrame(
        {
            "task": ["task_c"] * 4,
            "repeat": [0, 0, 1, 1],
            "outer_fold": [0, 0, 0, 1],
            "class": ["tumor"] * 4,
            "gene_symbol": ["GENE1"] * 4,
            "coefficient": [1.0, 1.1, 0.9, 1.2],
            "selected_nonzero": [True] * 4,
        }
    )
    summary = feature_stability(coefficients, total_outer_folds=4)
    assert int(summary.loc[0, "selected_folds"]) == 3
    assert summary.loc[0, "selection_frequency"] == 0.75
