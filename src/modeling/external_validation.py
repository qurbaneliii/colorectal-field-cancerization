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
    external_negative_label: str = "normal_colon",
    external_positive_label: str = "primary_tumor",
) -> tuple[pd.DataFrame, dict[str, float]]:
    classes = np.asarray(model.named_steps["model"].classes_)
    probability = model.predict_proba(x)
    if "tumor" not in classes:
        raise ValueError(f"Locked model has unexpected classes: {classes}")
    tumor_index = int(np.flatnonzero(classes == "tumor")[0])
    tumor_probability = probability[:, tumor_index]
    predicted = np.where(
        tumor_probability >= threshold, external_positive_label, external_negative_label
    )
    external_classes = np.asarray([external_negative_label, external_positive_label])
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
    frame[f"probability_{external_negative_label}"] = probability[:, 0]
    frame[f"probability_{external_positive_label}"] = probability[:, 1]
    metrics = classification_metrics(
        y_true, predicted, probability, external_classes, include_calibration=True
    )
    return frame, metrics
