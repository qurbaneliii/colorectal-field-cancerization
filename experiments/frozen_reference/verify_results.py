"""Independent sklearn metric replay and reference/score provenance checks."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             brier_score_loss, confusion_matrix, f1_score, log_loss,
                             roc_auc_score)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.data.provenance import file_sha256
from src.preprocessing.frozen_quantile import FrozenQuantile


def metrics(y, scores, threshold):
    predicted = scores >= threshold
    tn, fp, fn, tp = confusion_matrix(y, predicted, labels=[False, True]).ravel()
    return {"f1_macro": f1_score(y, predicted, average="macro"),
            "balanced_accuracy": balanced_accuracy_score(y, predicted),
            "roc_auc": roc_auc_score(y, scores), "pr_auc": average_precision_score(y, scores),
            "brier_score": brier_score_loss(y, scores), "log_loss": log_loss(y, scores),
            "sensitivity": tp/(tp+fn), "specificity": tn/(tn+fp),
            "tn": tn, "fp": fp, "fn": fn, "tp": tp}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    execution = json.loads((args.output / "execution.json").read_text())
    for name, sha in execution["outputs"].items():
        assert file_sha256(args.output / name.replace("\\", "/")) == sha, name
    predictions = pd.read_csv(args.output / "heldout_predictions.csv")
    keys = ["task", "repeat", "fold", "method"]
    expected = pd.read_csv(args.output / "fold_metrics.csv").set_index(keys)
    rows = []
    for key, part in predictions.groupby(keys):
        for metric, value in metrics(part.observed, part.score, .5).items():
            np.testing.assert_allclose(value, expected.loc[key, metric], atol=1e-10, rtol=1e-8)
            rows.append({"scope": str(key), "metric": metric, "value": value, "status": "PASS"})
    for _, boundary in pd.read_csv(args.output / "reference_boundaries.csv").iterrows():
        train, test = boundary.training_ids.split("|"), boundary.test_ids.split("|")
        assert set(train).isdisjoint(test)
        path = args.output / "locks" / f"{boundary.task}_{boundary['repeat']}_{boundary.fold}_reference.json"
        reference = FrozenQuantile.load(path)
        assert reference.training_ids == tuple(train) and reference.digest() == boundary.reference_sha256
    external = pd.read_csv(args.output / "external_predictions.csv")
    expected = pd.read_csv(args.output / "external_metrics.csv").set_index("method")
    for method, part in external.groupby("method"):
        assert len(part) == 233 and part.patient_id.nunique() == 190
        for metric, value in metrics(part.observed, part.score, .41).items():
            np.testing.assert_allclose(value, expected.loc[method, metric], atol=1e-10, rtol=1e-8)
            rows.append({"scope": method, "metric": metric, "value": value, "status": "PASS"})
    pd.DataFrame(rows).to_csv(args.output / "independent_metric_replay.csv", index=False)
    print(f"PASS: {len(rows)} independently recalculated metrics and all execution hashes")


if __name__ == "__main__":
    main()
