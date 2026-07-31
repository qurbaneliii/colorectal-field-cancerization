from __future__ import annotations

import json
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV

from src.modeling.evaluation import classification_metrics, prediction_probabilities
from src.modeling.pipelines import build_pipeline, parameter_grid, selected_feature_names
from src.modeling.splitters import stratified_group_splits


def run_nested_cv(
    x: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    sample_ids: np.ndarray,
    feature_names: list[str],
    task: str,
    models: list[str],
    seeds: list[int],
    outer_splits: int,
    inner_splits: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prediction_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    assignment_rows: list[dict[str, object]] = []
    coefficient_rows: list[dict[str, object]] = []
    for repeat, seed in enumerate(seeds):
        outer = stratified_group_splits(y, groups, outer_splits, seed)
        for fold_number, (train_idx, test_idx) in enumerate(outer):
            for idx in train_idx:
                assignment_rows.append(
                    {
                        "task": task,
                        "repeat": repeat,
                        "outer_fold": fold_number,
                        "sample_id": sample_ids[idx],
                        "patient_id": groups[idx],
                        "split": "train",
                    }
                )
            for idx in test_idx:
                assignment_rows.append(
                    {
                        "task": task,
                        "repeat": repeat,
                        "outer_fold": fold_number,
                        "sample_id": sample_ids[idx],
                        "patient_id": groups[idx],
                        "split": "test",
                    }
                )
            inner = stratified_group_splits(
                y[train_idx], groups[train_idx], inner_splits, seed + 1000 + fold_number
            )
            for model_name in models:
                grid = GridSearchCV(
                    build_pipeline(model_name, seed),
                    parameter_grid(model_name),
                    scoring="f1_macro",
                    cv=inner,
                    n_jobs=-1,
                    refit=True,
                    error_score="raise",
                )
                grid.fit(x[train_idx], y[train_idx], groups=groups[train_idx])
                fitted = grid.best_estimator_
                predicted = fitted.predict(x[test_idx])
                classes, probability = prediction_probabilities(fitted, x[test_idx])
                metrics = classification_metrics(y[test_idx], predicted, probability, classes)
                metric_rows.append(
                    {
                        "task": task,
                        "model": model_name,
                        "repeat": repeat,
                        "outer_fold": fold_number,
                        "n_train": len(train_idx),
                        "n_test": len(test_idx),
                        "best_inner_f1_macro": grid.best_score_,
                        "best_parameters": json.dumps(grid.best_params_, sort_keys=True),
                        **metrics,
                    }
                )
                for row_index, sample_index in enumerate(test_idx):
                    row = {
                        "task": task,
                        "model": model_name,
                        "repeat": repeat,
                        "outer_fold": fold_number,
                        "sample_id": sample_ids[sample_index],
                        "patient_id": groups[sample_index],
                        "y_true": y[sample_index],
                        "y_pred": predicted[row_index],
                    }
                    row.update(
                        {
                            f"probability_{label}": float(probability[row_index, class_index])
                            for class_index, label in enumerate(classes)
                        }
                    )
                    prediction_rows.append(row)
                if model_name == "elastic_net":
                    genes = selected_feature_names(fitted, feature_names)
                    estimator = fitted.named_steps["model"]
                    coefficients = estimator.coef_
                    coefficient_classes = (
                        estimator.classes_
                        if coefficients.shape[0] > 1
                        else np.asarray([estimator.classes_[1]])
                    )
                    for class_index, label in enumerate(coefficient_classes):
                        for gene, coefficient in zip(genes, coefficients[class_index], strict=True):
                            coefficient_rows.append(
                                {
                                    "task": task,
                                    "repeat": repeat,
                                    "outer_fold": fold_number,
                                    "class": label,
                                    "gene_symbol": gene,
                                    "coefficient": float(coefficient),
                                    "selected_nonzero": bool(abs(coefficient) > 1e-12),
                                }
                            )
    return (
        pd.DataFrame(prediction_rows),
        pd.DataFrame(metric_rows),
        pd.DataFrame(assignment_rows),
        pd.DataFrame(coefficient_rows),
    )


def modal_best_parameters(metrics: pd.DataFrame, task: str, model: str) -> dict[str, object]:
    rows = metrics[(metrics["task"].eq(task)) & (metrics["model"].eq(model))]
    if rows.empty:
        raise ValueError(f"No metrics for {task}/{model}")
    value = Counter(rows["best_parameters"]).most_common(1)[0][0]
    return json.loads(value)
