from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import expit, softmax
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    log_loss,
    precision_recall_fscore_support,
    roc_auc_score,
)
from sklearn.preprocessing import label_binarize


def prediction_probabilities(model, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    estimator = model.named_steps["model"]
    classes = np.asarray(estimator.classes_)
    if hasattr(model, "predict_proba"):
        return classes, model.predict_proba(x)
    decision = model.decision_function(x)
    if decision.ndim == 1:
        positive = expit(decision)
        probability = np.column_stack([1 - positive, positive])
    else:
        probability = softmax(decision, axis=1)
    return classes, probability


def classification_metrics(
    y_true: np.ndarray, y_pred: np.ndarray, probability: np.ndarray, classes: np.ndarray
) -> dict[str, float]:
    probability = np.asarray(probability, dtype=float)
    probability = probability / probability.sum(axis=1, keepdims=True)
    log_probability = np.clip(probability, 1e-15, 1.0)
    log_probability = log_probability / log_probability.sum(axis=1, keepdims=True)
    result: dict[str, float] = {
        "f1_macro": f1_score(y_true, y_pred, average="macro"),
        "f1_weighted": f1_score(y_true, y_pred, average="weighted"),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "log_loss": log_loss(y_true, log_probability, labels=classes),
    }
    precision, recall, f1, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=classes, zero_division=0
    )
    for label, p, r, f in zip(classes, precision, recall, f1, strict=True):
        safe = str(label)
        result[f"precision_{safe}"] = float(p)
        result[f"recall_{safe}"] = float(r)
        result[f"f1_{safe}"] = float(f)
    if len(classes) == 2:
        positive = classes[1]
        binary = (y_true == positive).astype(int)
        result["roc_auc"] = (
            roc_auc_score(binary, probability[:, 1]) if np.unique(binary).size == 2 else np.nan
        )
        result["pr_auc"] = (
            average_precision_score(binary, probability[:, 1])
            if np.unique(binary).size == 2
            else np.nan
        )
        result["brier_score"] = brier_score_loss(binary, probability[:, 1])
        result["sensitivity"] = result[f"recall_{positive}"]
        negative = classes[0]
        result["specificity"] = result[f"recall_{negative}"]
    else:
        binary = label_binarize(y_true, classes=classes)
        result["roc_auc_ovr_macro"] = (
            roc_auc_score(binary, probability, average="macro", multi_class="ovr")
            if set(classes).issubset(set(y_true))
            else np.nan
        )
        result["brier_score_multiclass"] = float(
            np.mean(np.sum((probability - binary) ** 2, axis=1))
        )
    return result


def grouped_bootstrap_metrics(
    predictions: pd.DataFrame,
    metrics: list[str],
    iterations: int,
    seed: int,
) -> list[dict[str, float]]:
    rng = np.random.default_rng(seed)
    groups = predictions["patient_id"].astype(str).unique()
    classes = np.array(sorted(predictions["y_true"].unique()))
    probability_columns = [f"probability_{label}" for label in classes]
    group_values = predictions["patient_id"].astype(str).to_numpy()
    group_indices = {group: np.flatnonzero(group_values == group) for group in groups}
    y_true = predictions["y_true"].to_numpy()
    y_pred = predictions["y_pred"].to_numpy()
    probabilities = predictions[probability_columns].to_numpy()
    values: dict[str, list[float]] = {metric: [] for metric in metrics}
    for _ in range(iterations):
        sampled = rng.choice(groups, size=len(groups), replace=True)
        sampled_indices = np.concatenate([group_indices[group] for group in sampled])
        measured = classification_metrics(
            y_true[sampled_indices],
            y_pred[sampled_indices],
            probabilities[sampled_indices],
            classes,
        )
        for metric in metrics:
            value = measured.get(metric, np.nan)
            if np.isfinite(value):
                values[metric].append(float(value))
    rows = []
    for metric, observed in values.items():
        array = np.asarray(observed)
        rows.append(
            {
                "metric": metric,
                "mean": float(array.mean()),
                "median": float(np.median(array)),
                "standard_deviation": float(array.std(ddof=1)),
                "ci_lower": float(np.quantile(array, 0.025)),
                "ci_upper": float(np.quantile(array, 0.975)),
                "bootstrap_iterations": int(len(array)),
            }
        )
    return rows


def grouped_bootstrap_ci(
    predictions: pd.DataFrame,
    metric: str,
    iterations: int,
    seed: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    groups = predictions["patient_id"].astype(str).unique()
    classes = np.array(sorted(predictions["y_true"].unique()))
    probability_columns = [f"probability_{label}" for label in classes]
    group_values = predictions["patient_id"].astype(str).to_numpy()
    group_indices = {group: np.flatnonzero(group_values == group) for group in groups}
    y_true = predictions["y_true"].to_numpy()
    y_pred = predictions["y_pred"].to_numpy()
    probabilities = predictions[probability_columns].to_numpy()
    values = []
    for _ in range(iterations):
        sampled = rng.choice(groups, size=len(groups), replace=True)
        sampled_indices = np.concatenate([group_indices[group] for group in sampled])
        measured = classification_metrics(
            y_true[sampled_indices],
            y_pred[sampled_indices],
            probabilities[sampled_indices],
            classes,
        )
        if metric in measured and np.isfinite(measured[metric]):
            values.append(measured[metric])
    array = np.asarray(values)
    return {
        "metric": metric,
        "mean": float(array.mean()),
        "median": float(np.median(array)),
        "standard_deviation": float(array.std(ddof=1)),
        "ci_lower": float(np.quantile(array, 0.025)),
        "ci_upper": float(np.quantile(array, 0.975)),
        "bootstrap_iterations": int(len(array)),
    }
