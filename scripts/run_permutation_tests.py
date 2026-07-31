from __future__ import annotations

import argparse
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

from src.data.provenance import expression_path, result_root
from src.modeling.nested_cv import modal_best_parameters
from src.modeling.pipelines import build_pipeline
from src.modeling.splitters import stratified_group_splits


TASKS = {
    "task_a_three_class": ["healthy", "adjacent_normal", "tumor"],
    "task_b_field_effect": ["healthy", "adjacent_normal"],
    "task_c_tumor_vs_adjacent": ["adjacent_normal", "tumor"],
}


def permute_labels(
    task: str, y: np.ndarray, groups: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    result = y.copy()
    unique_groups = np.unique(groups)
    group_indices = [np.flatnonzero(groups == group) for group in unique_groups]
    if task == "task_b_field_effect":
        if any(len(indices) != 1 for indices in group_indices):
            raise AssertionError("Task B permutation requires one eligible sample per donor group")
        return rng.permutation(result)
    if task == "task_c_tumor_vs_adjacent":
        if any(len(indices) != 2 for indices in group_indices):
            raise AssertionError("Task C permutation requires complete two-sample patient pairs")
        for indices in group_indices:
            if rng.integers(0, 2):
                result[indices] = result[indices[::-1]]
        return result
    if task == "task_a_three_class":
        # Conditional null: healthy labels remain fixed; tumor/adjacent labels are
        # independently flipped inside complete cancer-patient pairs.
        for indices in group_indices:
            if len(indices) == 2 and rng.integers(0, 2):
                result[indices] = result[indices[::-1]]
        return result
    raise ValueError(task)


def evaluate_labels(x, labels, splits, parameters, seed, modeling) -> float:
    fold_scores = []
    for train, test in splits:
        model = build_pipeline("elastic_net", seed, modeling).set_params(
            **parameters, memory=None
        )
        model.fit(x[train], labels[train])
        fold_scores.append(f1_score(labels[test], model.predict(x[test]), average="macro"))
    return float(np.mean(fold_scores))


def one_permutation(permutation_seed, task, x, y, groups, splits, parameters, model_seed, modeling):
    labels = permute_labels(task, y, groups, np.random.default_rng(permutation_seed))
    return evaluate_labels(x, labels, splits, parameters, model_seed, modeling)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--provenance", choices=["raw_cel_rma", "geo_deposited_series_matrix"])
    args = parser.parse_args()
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    modeling = config["modeling"]
    provenance = args.provenance or modeling["expression_provenance"]
    metrics_root = result_root(ROOT / paths["metrics"], provenance)
    metrics = pd.read_csv(metrics_root / "nested_cv_metrics.csv")
    expression = pd.read_parquet(
        expression_path(ROOT / paths["processed"], "GSE44076", "gene", provenance)
    ).set_index("gene_symbol")
    metadata = pd.read_csv(
        ROOT / paths["metadata"] / "gse44076_samples.csv", dtype={"patient_id": str}
    )
    iterations = int(modeling["permutation_iterations"])
    seed = int(config["project"]["random_seed"])
    null_hypotheses = {
        "task_a_three_class": "conditional null of no tumor-versus-adjacent information while healthy labels remain fixed",
        "task_b_field_effect": "global group-level null of no healthy-versus-adjacent association",
        "task_c_tumor_vs_adjacent": "paired null of exchangeable tumor/adjacent labels within patients",
    }
    rows = []
    for task_index, (task, labels) in enumerate(TASKS.items(), start=1):
        selected = metadata[
            metadata["inclusion_status"].eq("included")
            & metadata["tissue_class"].isin(labels)
        ]
        x = expression[selected["geo_accession"]].T.to_numpy(dtype=np.float32)
        y = selected["tissue_class"].to_numpy()
        groups = selected["donor_or_patient_group"].astype(str).to_numpy()
        splits = stratified_group_splits(y, groups, int(modeling["outer_splits"]), seed)
        parameters = modal_best_parameters(metrics, task, "elastic_net")
        observed = evaluate_labels(x, y, splits, parameters, seed, modeling)
        null_scores = Parallel(n_jobs=-1, verbose=0)(
            delayed(one_permutation)(
                seed + 10000 * task_index + iteration,
                task,
                x,
                y,
                groups,
                splits,
                parameters,
                seed,
                modeling,
            )
            for iteration in range(iterations)
        )
        p_value = (1 + int(np.sum(np.asarray(null_scores) >= observed))) / (iterations + 1)
        rows.append(
            {
                "task": task,
                "null_hypothesis": null_hypotheses[task],
                "observed_mean_macro_f1": observed,
                "null_mean": float(np.mean(null_scores)),
                "null_standard_deviation": float(np.std(null_scores, ddof=1)),
                "null_95_percentile": float(np.quantile(null_scores, 0.95)),
                "permutation_p_value": p_value,
                "iterations": iterations,
                "minimum_attainable_p_value": 1 / (iterations + 1),
                "analysis_provenance": provenance,
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(metrics_root / "group_preserving_permutation_tests.csv", index=False)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
