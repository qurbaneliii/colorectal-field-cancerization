from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from joblib import Parallel, delayed
from sklearn.metrics import f1_score

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.modeling.nested_cv import modal_best_parameters
from src.modeling.pipelines import build_pipeline
from src.modeling.splitters import stratified_group_splits


TASKS = {
    "task_a_three_class": ["healthy", "adjacent_normal", "tumor"],
    "task_b_field_effect": ["healthy", "adjacent_normal"],
    "task_c_tumor_vs_adjacent": ["adjacent_normal", "tumor"],
}


def permute_labels(y: np.ndarray, groups: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    result = y.copy()
    unique_groups = np.unique(groups)
    group_indices = [np.flatnonzero(groups == group) for group in unique_groups]
    by_size: dict[int, list[np.ndarray]] = {}
    for indices in group_indices:
        by_size.setdefault(len(indices), []).append(indices)
    for size, blocks in by_size.items():
        if size == 1:
            labels = np.array([result[index[0]] for index in blocks], dtype=object)
            labels = rng.permutation(labels)
            for index, label in zip(blocks, labels, strict=True):
                result[index[0]] = label
        else:
            # Patient blocks remain intact; within paired blocks the tissue labels
            # are independently flipped. This preserves class counts and grouping.
            for index in blocks:
                result[index] = rng.permutation(result[index])
    return result


def evaluate_labels(
    x: np.ndarray,
    labels: np.ndarray,
    splits: list[tuple[np.ndarray, np.ndarray]],
    parameters: dict[str, object],
    seed: int,
) -> float:
    fold_scores = []
    for train, test in splits:
        model = build_pipeline("elastic_net", seed).set_params(**parameters, memory=None)
        model.fit(x[train], labels[train])
        fold_scores.append(f1_score(labels[test], model.predict(x[test]), average="macro"))
    return float(np.mean(fold_scores))


def one_permutation(
    permutation_seed: int,
    x: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    splits: list[tuple[np.ndarray, np.ndarray]],
    parameters: dict[str, object],
    model_seed: int,
) -> float:
    permuted = permute_labels(y, groups, np.random.default_rng(permutation_seed))
    return evaluate_labels(x, permuted, splits, parameters, model_seed)


def main() -> None:
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    metrics = pd.read_csv(ROOT / paths["metrics"] / "nested_cv_metrics.csv")
    expression = pd.read_parquet(
        ROOT / paths["processed"] / "GSE44076_gene_expression.parquet"
    ).set_index("gene_symbol")
    metadata = pd.read_csv(
        ROOT / paths["metadata"] / "gse44076_samples.csv", dtype={"patient_id": str}
    )
    iterations = int(config["modeling"]["permutation_iterations"])
    seed = int(config["project"]["random_seed"])
    rows = []
    for task, labels in TASKS.items():
        selected = metadata[
            metadata["inclusion_status"].eq("included") & metadata["tissue_class"].isin(labels)
        ]
        x = expression[selected["geo_accession"]].T.to_numpy(dtype=np.float32)
        y = selected["tissue_class"].to_numpy()
        groups = selected["donor_or_patient_group"].astype(str).to_numpy()
        splits = stratified_group_splits(y, groups, config["modeling"]["outer_splits"], seed)
        parameters = modal_best_parameters(metrics, task, "elastic_net")
        observed = evaluate_labels(x, y, splits, parameters, seed)
        null_scores = Parallel(n_jobs=-1, verbose=0)(
            delayed(one_permutation)(
                seed + 10000 * (len(rows) + 1) + iteration,
                x,
                y,
                groups,
                splits,
                parameters,
                seed,
            )
            for iteration in range(iterations)
        )
        p_value = (1 + int(np.sum(np.asarray(null_scores) >= observed))) / (iterations + 1)
        rows.append(
            {
                "task": task,
                "observed_mean_macro_f1": observed,
                "null_mean": float(np.mean(null_scores)),
                "null_standard_deviation": float(np.std(null_scores, ddof=1)),
                "null_95_percentile": float(np.quantile(null_scores, 0.95)),
                "permutation_p_value": p_value,
                "iterations": iterations,
                "permutation_unit": "patient/donor group preserving",
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(ROOT / paths["metrics"] / "group_preserving_permutation_tests.csv", index=False)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
