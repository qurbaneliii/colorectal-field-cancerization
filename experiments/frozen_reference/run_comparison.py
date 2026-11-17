"""Run the prespecified fixed-panel real-CEL sensitivity; never tune on external labels."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.optimize import minimize
from scipy.special import expit, logit
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.verify_phase_ii import binary_metrics
from src.data.provenance import file_sha256
from src.preprocessing.frozen_quantile import FrozenQuantile


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False,
                              default=lambda x: x.item()) + "\n", encoding="utf-8")


def fit_predict(train, labels, test, model_path):
    model = make_pipeline(StandardScaler(), LogisticRegression(
        C=1, solver="lbfgs", max_iter=1000, class_weight="balanced", random_state=44076))
    model.fit(train, labels)
    assert list(model.classes_) == [False, True]
    scaler, logistic = model.steps[0][1], model.steps[1][1]
    lock = {"features": train.columns.tolist(), "mean": scaler.mean_.tolist(),
            "scale": scaler.scale_.tolist(), "coefficients": logistic.coef_[0].tolist(),
            "intercept": float(logistic.intercept_[0]), "positive_class": True}
    dump(model_path, lock)
    return model, model.predict_proba(test)[:, 1]


def calibration(observed, scores):
    x = logit(np.clip(scores, 1e-12, 1-1e-12))
    def loss(parameters):
        z = parameters[0] + parameters[1] * x
        return np.sum(np.logaddexp(0, z) - observed * z)
    result = minimize(loss, [0., 1.], method="BFGS")
    # Separation can make slope estimates unstable; retain the convergence flag.
    return {"calibration_intercept": float(result.x[0]), "calibration_slope": float(result.x[1]),
            "calibration_converged": bool(result.success),
            "mean_probability": float(expit(x).mean()), "observed_fraction": float(np.mean(observed))}


def paired_bootstrap(frame, threshold, iterations=1000):
    # Repeat-averaged predictions retain each patient as the sampling unit.
    frame = frame.groupby(["sample_id", "patient_id", "observed", "method"], as_index=False).score.mean()
    wide = frame.pivot(index=["sample_id", "patient_id", "observed"], columns="method", values="score")
    donors = wide.index.get_level_values("patient_id").to_numpy()
    observed = wide.index.get_level_values("observed").to_numpy(bool)
    unique = pd.unique(donors)
    groups = {donor: np.flatnonzero(donors == donor) for donor in unique}
    rng = np.random.default_rng(44076)
    draws = []
    for _ in range(iterations):
        idx = np.concatenate([groups[d] for d in rng.choice(unique, len(unique), replace=True)])
        if len(np.unique(observed[idx])) != 2:
            continue
        measured = {name: binary_metrics(observed[idx], wide[name].to_numpy()[idx] >= threshold,
                                         wide[name].to_numpy()[idx]) for name in wide.columns}
        draws.append({metric: measured["frozen"][metric] - measured["joint_rma"][metric]
                      for metric in measured["frozen"]})
    return [{"metric": metric, "delta_lower": np.quantile([d[metric] for d in draws], .025),
             "delta_upper": np.quantile([d[metric] for d in draws], .975),
             "valid_draws": len(draws), "estimand": "conditional_frozen_minus_joint_rma"}
            for metric in draws[0]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("archive_cache", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    locks = args.output / "locks"
    locks.mkdir()
    protocol_path = Path("experiments/frozen_reference/protocol.json")
    protocol = json.loads(protocol_path.read_text())
    pins = json.loads(Path("results/ai_audit/current_geo_archive_pins.json").read_text())
    for filename, expected in pins["archives"].items():
        assert file_sha256(args.archive_cache / filename) == expected, filename
    inputs, grids, ranks, maps, baseline, meta = [], {}, {}, {}, {}, {}
    for accession in ["GSE44076", "GSE41258"]:
        source = Path(f"data/metadata/{accession.lower()}_samples.csv")
        meta[accession] = pd.read_csv(source, dtype={"patient_id": str}).query(
            "inclusion_status == 'included'").set_index("geo_accession")
        mapping = Path(f"data/metadata/{accession}_probe_gene_mapping_raw_cel_rma.csv")
        maps[accession] = pd.read_csv(mapping, keep_default_na=False)
        members = pd.read_csv(f"results/ai_audit/{accession.lower()}_raw_integrity_members.csv")
        expected = members.set_index("geo_accession").gzip_sha256
        per_sample = []
        for gsm in meta[accession].index:
            prefix = args.cache / accession / gsm
            identity = pd.read_csv(str(prefix) + "_input.csv").iloc[0]
            assert identity.cel_sha256 == expected[gsm], gsm
            assert identity.code_sha256 == file_sha256(Path("experiments/frozen_reference/extract_single_array.R"))
            rank_path, grid_path = Path(str(prefix) + "_ranks.csv"), Path(str(prefix) + "_quantiles.bin")
            values = pd.read_csv(rank_path).set_index("probe_id").percentile.rename(gsm)
            assert values.index.is_unique and values.between(0, 1).all()
            per_sample.append(values)
            grids[gsm] = np.fromfile(grid_path, dtype="<f8")
            assert len(grids[gsm]) == protocol["quantile_grid_points"]
            inputs.append({"accession": accession, "sample_id": gsm,
                           "cel_sha256": identity.cel_sha256, "ranks_sha256": file_sha256(rank_path),
                           "quantiles_sha256": file_sha256(grid_path), "pm_count": int(identity.pm_count)})
        ranks[accession] = pd.concat(per_sample, axis=1)
        bp = args.baseline / f"{accession}_gene_expression_raw_cel_rma.parquet"
        baseline[accession] = pd.read_parquet(bp).set_index("gene_symbol").T
        inputs.append({"accession": accession, "sample_id": "baseline_and_mapping",
                       "baseline_sha256": file_sha256(bp), "mapping_sha256": file_sha256(mapping),
                       "metadata_sha256": file_sha256(source)})
    pd.DataFrame(inputs).to_csv(args.output / "inputs.csv", index=False)
    split_path = args.bundle / "results/metrics/fold_assignments.csv"
    assignments = pd.read_csv(split_path, dtype={"patient_id": str})
    tasks = {"task_b_field_effect": (protocol["task_b_panel"], "adjacent_normal"),
             "task_c_tumor_vs_adjacent": (protocol["task_c_panel"], "tumor")}
    # Explicitly fail if the archived task identifiers differ.
    assert set(tasks).issubset(assignments.task.unique()), assignments.task.unique()
    predictions, metrics, invariance, boundaries, stability = [], [], [], [], []
    for task, (panel, positive) in tasks.items():
        mapping = maps["GSE44076"].loc[lambda x: x.gene_symbol.isin(panel)]
        for (repeat, fold), split in assignments.loc[assignments.task.eq(task)].groupby(["repeat", "outer_fold"]):
            train = split.loc[split.split.eq("train"), "sample_id"].tolist()
            test = split.loc[split.split.eq("test"), "sample_id"].tolist()
            assert set(split.loc[split.split.eq("train"), "patient_id"]).isdisjoint(
                split.loc[split.split.eq("test"), "patient_id"])
            ref = FrozenQuantile.fit(grids, train, test)
            name = f"{task}_{repeat}_{fold}"
            ref.save(locks / f"{name}_reference.json")
            transformed = ref.transform(ranks["GSE44076"][train + test], mapping)[panel]
            y = meta["GSE44076"].loc[train, "tissue_class"].eq(positive)
            truth = meta["GSE44076"].loc[test, "tissue_class"].eq(positive)
            boundaries.append({"task": task, "repeat": repeat, "fold": fold,
                               "training_ids": "|".join(train), "test_ids": "|".join(test),
                               "reference_sha256": ref.digest()})
            for method, matrix in [("frozen", transformed), ("joint_rma", baseline["GSE44076"][panel])]:
                model_path = locks / f"{name}_{method}_model.json"
                model, scores = fit_predict(matrix.loc[train], y, matrix.loc[test], model_path)
                metrics.append({"task": task, "repeat": repeat, "fold": fold, "method": method,
                                **binary_metrics(truth.to_numpy(), scores >= .5, scores),
                                **calibration(truth.to_numpy(), scores)})
                for gsm, observed, score in zip(test, truth, scores):
                    predictions.append({"task": task, "repeat": repeat, "fold": fold,
                                        "sample_id": gsm, "patient_id": meta["GSE44076"].loc[gsm, "patient_id"],
                                        "observed": bool(observed), "method": method, "score": score,
                                        "model_sha256": file_sha256(model_path)})
                if method == "frozen":
                    # Every held-out specimen is independently transformed once again.
                    for i, gsm in enumerate(test):
                        alone = ref.transform(ranks["GSE44076"][[gsm]], mapping)[panel]
                        probability = model.predict_proba(alone)[0, 1]
                        difference = float(np.max(np.abs(alone.loc[gsm] - transformed.loc[gsm])))
                        np.testing.assert_allclose(probability, scores[i], atol=1e-12, rtol=1e-12)
                        assert difference <= 1e-12
                        invariance.append({"task": task, "repeat": repeat, "fold": fold,
                                           "sample_id": gsm, "gene_max_abs": difference,
                                           "probability_abs": abs(probability-scores[i]),
                                           "decision_equal": (probability >= .5) == (scores[i] >= .5)})
            for gene in panel:
                stability.append({"task": task, "repeat": repeat, "fold": fold, "gene": gene,
                                  "heldout_spearman": spearmanr(transformed.loc[test, gene],
                                                               baseline["GSE44076"].loc[test, gene]).statistic})
            print(task, repeat, fold, "PASS", flush=True)
    frame = pd.DataFrame(predictions)
    frame.to_csv(args.output / "heldout_predictions.csv", index=False)
    pd.DataFrame(metrics).to_csv(args.output / "fold_metrics.csv", index=False)
    pd.DataFrame(invariance).to_csv(args.output / "invariance.csv", index=False)
    pd.DataFrame(boundaries).to_csv(args.output / "reference_boundaries.csv", index=False)
    pd.DataFrame(stability).to_csv(args.output / "panel_stability.csv", index=False)
    ci = []
    for task, group in frame.groupby("task"):
        ci += [{"task": task, **row} for row in paired_bootstrap(group, .5)]
    pd.DataFrame(ci).to_csv(args.output / "paired_conditional_intervals.csv", index=False)

    # Lock full-development reference and numeric model before reading external outcome labels.
    dev = meta["GSE44076"].query("tissue_class in ['adjacent_normal', 'tumor']")
    ref = FrozenQuantile.fit(grids, dev.index.tolist(), meta["GSE41258"].index)
    ref.save(locks / "external_reference.json")
    common = sorted(set(maps["GSE44076"].query("mapping_status == 'mapped_unique'").gene_symbol)
                    & set(maps["GSE41258"].query("mapping_status == 'mapped_unique'").gene_symbol))
    assert set(protocol["task_c_transport_panel"]).issubset(common)
    assert "FOXQ1" not in common
    (locks / "rank_universe.txt").write_text("\n".join(common) + "\n", encoding="utf-8")
    panel = protocol["task_c_transport_panel"]
    external_predictions, external_metrics = [], []
    for method in ["frozen", "joint_rma"]:
        matrices = {}
        for accession in ["GSE44076", "GSE41258"]:
            expression = ref.transform(ranks[accession], maps[accession]) if method == "frozen" else baseline[accession]
            matrices[accession] = expression[common].rank(axis=1, pct=True)[panel]
        model, scores = fit_predict(matrices["GSE44076"].loc[dev.index], dev.tissue_class.eq("tumor"),
                                    matrices["GSE41258"], locks / f"external_{method}_model.json")
        score_table = pd.DataFrame({"sample_id": matrices["GSE41258"].index, "score": scores})
        score_table.to_csv(args.output / f"external_{method}_locked_scores.csv", index=False)
        if method == "frozen":
            digest_before = ref.digest()
            for i in np.linspace(0, len(scores)-1, 5, dtype=int):
                gsm = score_table.sample_id.iloc[i]
                alone = ref.transform(ranks["GSE41258"][[gsm]], maps["GSE41258"])
                gene_abs = float(np.max(np.abs(alone.loc[gsm] - expression.loc[gsm])))
                gene_relative = float(np.max(np.abs(alone.loc[gsm] - expression.loc[gsm])
                                             / np.maximum(np.abs(alone.loc[gsm]), 1e-12)))
                pair = ranks["GSE41258"].iloc[:, [i, (i+1) % len(scores)]]
                probes_alone = ref.probe_values(ranks["GSE41258"][gsm].to_numpy())
                probes_batch = ref.probe_values(pair.to_numpy())[:, 0]
                probe_abs = float(np.max(np.abs(probes_alone - probes_batch)))
                p = model.predict_proba(alone[common].rank(axis=1, pct=True)[panel])[0, 1]
                np.testing.assert_allclose(p, scores[i], atol=1e-12, rtol=1e-12)
                assert max(gene_abs, gene_relative, probe_abs) <= 1e-12
                invariance.append({"task": "external_task_c", "repeat": -1, "fold": -1,
                                   "sample_id": gsm, "gene_max_abs": gene_abs,
                                   "gene_max_relative": gene_relative, "probeset_max_abs": probe_abs,
                                   "probability_abs": abs(p-scores[i]),
                                   "decision_equal": (p >= .41) == (scores[i] >= .41)})
            assert ref.digest() == digest_before
        external = meta["GSE41258"].loc[score_table.sample_id]
        observed = external.tissue_class.eq("primary_tumor").to_numpy()
        external_metrics.append({"method": method, "threshold": .41,
                                 **binary_metrics(observed, scores >= .41, scores),
                                 **calibration(observed, scores)})
        for gsm, donor, label, score in zip(external.index, external.patient_id, observed, scores):
            external_predictions.append({"sample_id": gsm, "patient_id": donor, "observed": bool(label),
                                         "method": method, "score": score})
    external_frame = pd.DataFrame(external_predictions)
    external_frame.to_csv(args.output / "external_predictions.csv", index=False)
    pd.DataFrame(external_metrics).to_csv(args.output / "external_metrics.csv", index=False)
    pd.DataFrame(paired_bootstrap(external_frame, .41)).to_csv(args.output / "external_paired_intervals.csv", index=False)
    pd.DataFrame(invariance).to_csv(args.output / "invariance.csv", index=False)
    output_hashes = {str(p.relative_to(args.output)): file_sha256(p) for p in args.output.rglob("*") if p.is_file()}
    dump(args.output / "execution.json", {"analysis_commit": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True).strip(), "python": platform.python_version(),
        "protocol_sha256": file_sha256(protocol_path), "split_sha256": file_sha256(split_path),
        "source_archives": pins, "shared_genes": len(common), "outputs": output_hashes,
        "status": "PASS", "scope": "fixed-panel conditional sensitivity; no nested feature discovery"})


if __name__ == "__main__":
    main()
