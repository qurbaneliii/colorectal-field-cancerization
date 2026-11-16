"""Evaluate a locked linear predictor on an eligible single-specimen cohort."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import expit, logit

from scripts.verify_phase_ii import binary_metrics
from src.data.provenance import file_sha256

LABELS = {"healthy", "cancer_associated_normal", "tumor"}
ANNOTATION = {"sample_id", "donor_id", "tissue", "site", "acquisition", "source_tissue", "source_record"}


def evaluate(lock_path, model_path, mapping_path, matrix_path, metadata_path, *, synthetic=False):
    inputs = {"lock": lock_path, "model": model_path, "mapping": mapping_path,
              "matrix": matrix_path, "metadata": metadata_path}
    hashes = {name: file_sha256(path) for name, path in inputs.items()}
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if (lock["mode"] == "synthetic") != synthetic or lock["mode"] not in {"synthetic", "research"}:
        raise ValueError("Lock and explicit synthetic mode must agree")
    if lock["training_accession"] == lock["evaluation_accession"]:
        raise ValueError("Evaluation accession must differ from development")
    for name in ("model", "mapping"):
        if hashes[name] != lock[f"{name}_sha256"]:
            raise ValueError(f"{name} digest differs from the lock")
    if lock["normalization"] != "frozen_single_sample":
        raise ValueError("A frozen single-sample representation is required")
    model = json.loads(model_path.read_text(encoding="utf-8"))
    features = model["features"]
    coefficients = np.asarray(model["coefficients"], dtype=float)
    threshold = float(lock["threshold"])
    if (not features or len(features) != len(set(features)) or coefficients.shape != (len(features),)
            or not np.isfinite(coefficients).all() or not np.isfinite(model["intercept"])
            or not 0 < threshold < 1):
        raise ValueError("Invalid locked linear model or cutoff")
    mapping = pd.read_csv(mapping_path, dtype=str, keep_default_na=False)
    if ({"feature", "column", "annotation"}.difference(mapping)
            or mapping["feature"].duplicated().any() or mapping["column"].duplicated().any()
            or set(mapping.feature) != set(features) or mapping.eq("").any().any()):
        raise ValueError("Every locked feature requires a unique annotated input column")
    ordered = mapping.set_index("feature").loc[features, "column"]
    expression = pd.read_csv(matrix_path, dtype={"sample_id": str}).set_index("sample_id")
    if not expression.index.is_unique or expression.index.hasnans or "sample_id" in ordered.to_list():
        raise ValueError("Invalid expression sample or feature identifiers")
    features_missing = set(ordered).difference(expression.columns)
    if features_missing:
        raise ValueError(f"Locked input columns absent: {sorted(features_missing)}")
    values = expression[ordered.tolist()].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Non-finite locked feature values")
    probability = pd.Series(expit(values @ coefficients + model["intercept"]), index=expression.index)

    # Outcome metadata enters only after the complete fixed score vector exists.
    metadata = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
    if (ANNOTATION.difference(metadata) or metadata[list(ANNOTATION)].apply(
            lambda column: column.str.strip().eq("").any()).any()
            or metadata.sample_id.duplicated().any() or not set(metadata.tissue).issubset(LABELS)
            or set(metadata.sample_id) != set(probability.index)):
        raise ValueError("Metadata requires explicit sample, donor, tissue, site and source identity")
    negative, positive = (lock["negative_tissue"], lock["positive_tissue"])
    if not {negative, positive}.issubset(LABELS) or negative == positive:
        raise ValueError("Invalid locked endpoint")
    eligible = metadata.loc[metadata.tissue.isin([negative, positive])
                            & metadata.site.eq(lock["site"])
                            & metadata.acquisition.eq(lock["acquisition"])].copy()
    if eligible.donor_id.duplicated().any():
        raise ValueError("Primary protocol requires one prespecified specimen per donor")
    eligible = eligible.sort_values("donor_id")
    if eligible.groupby("tissue").size().reindex([negative, positive], fill_value=0).min() < 2:
        raise ValueError("Endpoint requires both classes with at least two independent donors each")
    scores = probability.loc[eligible.sample_id].to_numpy()
    observed = eligible.tissue.eq(positive).to_numpy()
    predicted = scores >= threshold
    metrics = binary_metrics(observed, predicted, scores)
    logits = logit(np.clip(scores, 1e-12, 1-1e-12))
    calibration_intercept = brentq(lambda b: expit(logits + b).mean() - observed.mean(), -100, 100)
    rng = np.random.default_rng(lock["bootstrap_seed"])
    strata = [np.flatnonzero(observed == outcome) for outcome in (False, True)]
    if type(lock["bootstrap_iterations"]) is not int or lock["bootstrap_iterations"] < 2:
        raise ValueError("At least two bootstrap iterations are required")
    intervals = []
    for _ in range(lock["bootstrap_iterations"]):
        sample = np.concatenate([rng.choice(group, len(group), replace=True) for group in strata])
        intervals.append(binary_metrics(observed[sample], predicted[sample], scores[sample]))
    if hashes != {name: file_sha256(path) for name, path in inputs.items()}:
        raise ValueError("An evaluation input changed during execution")
    return {
        "mode": lock["mode"], "evaluated_at_utc": datetime.now(timezone.utc).isoformat(),
        "training_accession": lock["training_accession"], "evaluation_accession": lock["evaluation_accession"],
        "normalization": lock["normalization"], "feature_columns": ordered.to_dict(),
        "threshold": threshold, "sha256": hashes, "eligible_donors": len(eligible),
        "excluded_samples": len(metadata) - len(eligible), "metrics": metrics,
        "calibration_in_large": calibration_intercept,
        "conditional_donor_bootstrap_intervals": {
            metric: np.quantile([row[metric] for row in intervals], [0.025, 0.975]).tolist()
            for metric in ("roc_auc", "pr_auc", "sensitivity", "specificity", "brier_score", "f1_macro")
        },
        "score_quantiles_by_tissue": {
            label: np.quantile(scores[observed == is_positive], [0.1, 0.5, 0.9]).tolist()
            for label, is_positive in ((negative, False), (positive, True))
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("lock", "model", "mapping", "matrix", "metadata", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--synthetic", action="store_true")
    args = parser.parse_args()
    report = evaluate(args.lock, args.model, args.mapping, args.matrix, args.metadata, synthetic=args.synthetic)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, default=lambda item: item.item(), allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
