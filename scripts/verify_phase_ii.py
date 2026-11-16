"""Recompute held-out metrics from a hosted evidence bundle."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def binary_metrics(observed, predicted, scores):
    observed, predicted, scores = map(np.asarray, (observed, predicted, scores))
    n = len(scores)
    if n == 0 or not np.isfinite(scores).all() or ((scores < 0) | (scores > 1)).any():
        raise ValueError("Scores must be finite probabilities")
    tp = np.count_nonzero(observed & predicted)
    tn = np.count_nonzero(~observed & ~predicted)
    fp = np.count_nonzero(~observed & predicted)
    fn = np.count_nonzero(observed & ~predicted)
    positives, negatives = tp + fn, tn + fp
    if not positives or not negatives:
        raise ValueError("Both endpoint classes are required")
    ranks = rankdata(scores)
    order = np.argsort(-scores)
    group_end = np.r_[np.flatnonzero(np.diff(scores[order])), n - 1]
    cumulative = np.cumsum(observed[order])[group_end]
    probability = np.clip(np.where(observed, scores, 1 - scores), 1e-15, 1)
    return {
        "f1_macro": tp / (2 * tp + fp + fn) + tn / (2 * tn + fp + fn),
        "balanced_accuracy": (tp / positives + tn / negatives) / 2,
        "roc_auc": (ranks[observed].sum() - positives * (positives + 1) / 2)
        / (positives * negatives),
        "pr_auc": np.sum(np.diff(np.r_[0, cumulative]) * cumulative / (group_end + 1))
        / positives,
        "sensitivity": tp / positives,
        "specificity": tn / negatives,
        "brier_score": np.mean((observed - scores) ** 2),
        "log_loss": -np.log(probability).mean(),
        "tn": tn, "fp": fp, "fn": fn, "tp": tp,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.bundle / "results/metrics"
    rows = []

    def add(scope, key, values, reference, source):
        for metric, observed in values.items():
            expected = reference.get(metric)
            match = expected is not None
            if match:
                np.testing.assert_allclose(observed, expected, rtol=1e-7, atol=1e-10)
            rows.append({
                "scope": scope, "key": "|".join(map(str, key)), "metric": metric,
                "recalculated": observed, "reference": expected,
                "source": source.name, "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "comparison": "PASS" if match else "RECALCULATED",
            })

    keys = ["task", "repeat", "outer_fold"]
    assignments = pd.read_csv(root / "fold_assignments.csv", dtype={"patient_id": str})
    held_out = assignments.query("split == 'test'")
    assert len(assignments.groupby(keys)) == 45
    assert not held_out.duplicated(["task", "repeat", "sample_id"]).any()
    assert not assignments.duplicated(keys + ["sample_id"]).any()
    for _, group in assignments.groupby(keys):
        assert set(group.query("split == 'train'").patient_id).isdisjoint(
            group.query("split == 'test'").patient_id
        )
    source = root / "nested_cv_predictions.csv"
    predictions = pd.read_csv(source, dtype={"patient_id": str})
    assert len(predictions) == 5310
    assert not predictions.duplicated(keys + ["model", "sample_id"]).any()
    metrics = pd.read_csv(root / "nested_cv_metrics.csv").set_index(keys + ["model"])
    for key, fold in predictions.groupby(keys + ["model"]):
        expected = held_out.loc[(held_out.task == key[0]) & (held_out["repeat"] == key[1])
                               & (held_out.outer_fold == key[2])]
        assert set(zip(fold.sample_id, fold.patient_id)) == set(zip(expected.sample_id, expected.patient_id))
        classes = sorted(fold.y_true.unique())
        if len(classes) == 2:
            positive = classes[1]
            probabilities = fold[["probability_" + label for label in classes]].to_numpy(float)
            scores = probabilities[:, 1] / probabilities.sum(axis=1)
            add("nested_fold", key, binary_metrics(fold.y_true.eq(positive), fold.y_pred.eq(positive),
                                                  scores), metrics.loc[key], source)
        else:
            probabilities = fold[["probability_" + label for label in classes]].to_numpy(float)
            probabilities /= probabilities.sum(axis=1, keepdims=True)
            observed = fold.y_true.to_numpy()[:, None] == np.asarray(classes)
            predicted = fold.y_pred.to_numpy()[:, None] == np.asarray(classes)
            confusion = observed.astype(int).T @ predicted.astype(int)
            diagonal = confusion.diagonal()
            values = {
                "f1_macro": np.mean(2 * diagonal / (confusion.sum(axis=0) + confusion.sum(axis=1))),
                "balanced_accuracy": np.mean(diagonal / confusion.sum(axis=1)),
                "log_loss": -np.log(np.clip(probabilities[observed], 1e-15, 1)).mean(),
                "brier_score_multiclass": np.mean(((probabilities - observed) ** 2).sum(axis=1)),
                "roc_auc_ovr_macro": np.mean([
                    binary_metrics(observed[:, i], predicted[:, i], probabilities[:, i])["roc_auc"]
                    for i in range(len(classes))
                ]),
            }
            add("nested_fold", key, values, metrics.loc[key], source)
    source = root / "external_validation_predictions.csv"
    external = pd.read_csv(source, dtype={"patient_id": str})
    keys = ["representation", "evaluation_set", "threshold_policy"]
    metrics = pd.read_csv(root / "external_validation_metrics.csv").set_index(keys)
    for key, group in external.groupby(keys):
        scores = group.probability_primary_tumor
        assert group.locked_threshold.nunique() == 1
        decision = scores >= group.locked_threshold
        assert decision.eq(group.y_pred.eq("primary_tumor")).all()
        add("external", key, binary_metrics(group.y_true.eq("primary_tumor"), decision, scores),
            metrics.loc[key], source)
    primary = external.loc[
        external.representation.eq("within_sample_percentile_rank")
        & external.evaluation_set.eq("canonical_patient_tissue")
        & external.threshold_policy.eq("primary_gse44076_locked")
    ]
    assert len(primary) == 233 and primary.patient_id.nunique() == 190
    assert not primary.duplicated(["patient_id", "y_true"]).any()
    assert primary.locked_threshold.tolist() == [0.41] * 233
    groups = primary.patient_id.to_numpy()
    donors = pd.unique(groups)
    indices = {donor: np.flatnonzero(groups == donor) for donor in donors}
    observed = primary.y_true.eq("primary_tumor").to_numpy()
    probabilities = primary.probability_primary_tumor.to_numpy()
    selected = probabilities >= 0.41
    rng = np.random.default_rng(44076)
    draws = []
    for _ in range(1000):
        sample = np.concatenate([indices[donor] for donor in rng.choice(donors, len(donors), replace=True)])
        draws.append(binary_metrics(observed[sample], selected[sample], probabilities[sample]))
    reference = pd.read_csv(root / "external_validation_bootstrap_ci.csv").loc[lambda data:
        data.representation.eq("within_sample_percentile_rank")
        & data.evaluation_set.eq("canonical_patient_tissue")
        & data.threshold_policy.eq("primary_gse44076_locked")
    ].set_index("metric")
    for metric in ("roc_auc", "f1_macro", "specificity", "sensitivity", "balanced_accuracy",
                   "pr_auc", "brier_score", "log_loss"):
        lower, upper = np.quantile([row[metric] for row in draws], [0.025, 0.975])
        add("conditional_patient_bootstrap", [metric], {"ci_lower": lower, "ci_upper": upper},
            reference.loc[metric], source)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(f"Verified 45 fold assignments, {len(predictions)} held-out rows and {len(rows)} metric records")


if __name__ == "__main__":
    main()
