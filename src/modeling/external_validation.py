from __future__ import annotations

import numpy as np
import pandas as pd

from src.modeling.evaluation import classification_metrics


def locked_binary_predictions(
    model,
    x: np.ndarray,
    sample_ids: np.ndarray,
    patient_ids: np.ndarray,
    y_true: np.ndarray,
    representation: str,
    threshold: float = 0.5,
) -> tuple[pd.DataFrame, dict[str, float]]:
    classes = np.asarray(model.named_steps["model"].classes_)
    probability = model.predict_proba(x)
    if "tumor" not in classes:
        raise ValueError(f"Locked model has unexpected classes: {classes}")
    tumor_index = int(np.flatnonzero(classes == "tumor")[0])
    tumor_probability = probability[:, tumor_index]
    predicted = np.where(tumor_probability >= threshold, "tumor", "adjacent_normal")
    frame = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "patient_id": patient_ids,
            "y_true": y_true,
            "y_pred": predicted,
            "representation": representation,
            "locked_threshold": threshold,
        }
    )
    for class_index, label in enumerate(classes):
        frame[f"probability_{label}"] = probability[:, class_index]
    metrics = classification_metrics(
        y_true, predicted, probability, classes, include_calibration=True
    )
    return frame, metrics
