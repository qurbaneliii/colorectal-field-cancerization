from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.validation import validate_fold_assignments
from src.modeling.nested_cv import run_nested_cv


def main() -> None:
    rng = np.random.default_rng(44076)
    groups = np.repeat([f"patient_{index:02d}" for index in range(24)], 2)
    y = np.tile(np.array(["adjacent_normal", "tumor"]), 24)
    x = rng.normal(size=(48, 30))
    x[y == "tumor", :5] += 1.5
    modeling = {
        "variance_quantiles": [0.25],
        "max_features": [10],
        "elastic_net": {"C": [0.2], "l1_ratio": [0.5], "class_weight": ["balanced"]},
        "linear_svm": {"C": [0.2], "calibration_splits": 3},
        "random_forest": {
            "n_estimators": [20],
            "max_features": ["sqrt"],
            "min_samples_leaf": [1],
        },
    }
    predictions, metrics, assignments, coefficients = run_nested_cv(
        x,
        y,
        groups,
        np.array([f"sample_{index:02d}" for index in range(48)]),
        [f"GENE{index:03d}" for index in range(30)],
        "synthetic_task_c",
        ["elastic_net", "linear_svm"],
        [44076],
        4,
        3,
        modeling,
    )
    validate_fold_assignments(assignments)
    if not predictions.filter(like="probability_").notna().any().all():
        raise AssertionError("Synthetic smoke test did not produce calibrated probabilities")
    if metrics["f1_macro"].isna().any() or coefficients.empty:
        raise AssertionError("Synthetic smoke test outputs are incomplete")
    print(
        json.dumps(
            {
                "status": "PASS",
                "prediction_rows": len(predictions),
                "metric_rows": len(metrics),
                "assignment_rows": len(assignments),
                "coefficient_rows": len(coefficients),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
