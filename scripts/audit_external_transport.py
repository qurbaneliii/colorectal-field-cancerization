from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import brentq
from scipy.special import expit, logit
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    log_loss,
    roc_auc_score,
)

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.verify_external_results import ARTIFACTS, INTERVALS, verify
from src.data.cross_platform import within_sample_percentile_rank
from src.data.expression import aggregate_probe_expression, read_geo_series_matrix
from src.data.provenance import file_sha256
from src.modeling.evaluation import grouped_bootstrap_metrics
from src.modeling.external_result import EVALUATION_KEYS, PRIMARY_KEY


def replay(predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    canonical = predictions[predictions["evaluation_set"].eq("canonical_patient_tissue")]
    rows, strata = [], []
    for (representation, policy), frame in canonical.groupby(
        ["representation", "threshold_policy"]
    ):
        observed = frame["y_true"].eq("primary_tumor").to_numpy()
        positive = frame["y_pred"].eq("primary_tumor").to_numpy()
        p = frame["probability_primary_tumor"].to_numpy()
        logits = logit(np.clip(p, 1e-6, 1 - 1e-6))
        intercept = brentq(lambda offset: expit(logits + offset).mean() - observed.mean(), -50, 50)
        by_patient = frame.pivot(index="patient_id", columns="y_true", values="probability_primary_tumor")
        complete = by_patient.dropna()
        margin = complete["primary_tumor"] - complete["normal_colon"]
        for weighting in ("patient_tissue", "equal_patient"):
            weights = (
                np.ones(len(frame))
                if weighting == "patient_tissue"
                else 1 / frame.groupby("patient_id")["sample_id"].transform("size").to_numpy()
            )
            prevalence = np.average(observed, weights=weights)
            rows.append(
                dict(
                    representation=representation, policy=policy, weighting=weighting,
                    threshold=frame["locked_threshold"].iloc[0], samples=len(frame),
                    patients=frame["patient_id"].nunique(), prevalence=prevalence,
                    mean_probability=np.average(p, weights=weights),
                    true_negative=int((~observed & ~positive).sum()),
                    false_positive=int((~observed & positive).sum()),
                    false_negative=int((observed & ~positive).sum()),
                    true_positive=int((observed & positive).sum()),
                    specificity=np.average(~positive[~observed], weights=weights[~observed]),
                    sensitivity=np.average(positive[observed], weights=weights[observed]),
                    roc_auc=roc_auc_score(observed, p, sample_weight=weights),
                    average_precision=average_precision_score(observed, p, sample_weight=weights),
                    macro_f1=f1_score(observed, positive, average="macro", sample_weight=weights),
                    balanced_accuracy=balanced_accuracy_score(observed, positive, sample_weight=weights),
                    brier=brier_score_loss(observed, p, sample_weight=weights),
                    log_loss=log_loss(observed, p, sample_weight=weights),
                    prevalence_only_brier=prevalence * (1 - prevalence),
                    calibration_in_large=intercept if weighting == "patient_tissue" else None,
                    paired_patients=len(margin),
                    paired_concordance=((margin > 0).mean() + 0.5 * (margin == 0).mean()),
                )
            )
        if policy == "primary_gse44076_locked":
            for label, subset in frame.groupby("y_true"):
                scores = subset["probability_primary_tumor"]
                strata.append(
                    dict(
                        representation=representation, tissue=label, samples=len(subset),
                        minimum=scores.min(), q25=scores.quantile(0.25), median=scores.median(),
                        q75=scores.quantile(0.75), maximum=scores.max(),
                    )
                )
    return pd.DataFrame(rows), pd.DataFrame(strata)


def gene_shift(primary_path: Path, external_path: Path, genes: list[str]) -> pd.DataFrame:
    """Describe platform-level measurements without loading external outcomes."""
    gene_frames = {}
    probe_counts = {}
    for accession, path in (("GSE44076", primary_path), ("GSE41258", external_path)):
        mapping = pd.read_csv(
            ROOT / "data/metadata" / f"{accession}_probe_gene_mapping_raw_cel_rma.csv"
        ).query("mapping_status == 'mapped_unique'")
        probe = read_geo_series_matrix(path)
        gene_frames[accession] = aggregate_probe_expression(probe, mapping)
        probe_counts[accession] = mapping[mapping["probe_id"].isin(probe.index)]["gene_symbol"].value_counts()
    primary_meta = pd.read_csv(ROOT / "data/metadata/gse44076_samples.csv")
    development_samples = primary_meta.loc[
        primary_meta["tissue_class"].isin(["tumor", "adjacent_normal"]), "geo_accession"
    ]
    common = gene_frames["GSE44076"].index.intersection(gene_frames["GSE41258"].index)
    primary = gene_frames["GSE44076"].loc[common, development_samples]
    external = gene_frames["GSE41258"].loc[common]
    primary_ranks = within_sample_percentile_rank(primary.T.to_numpy())
    external_ranks = within_sample_percentile_rank(external.T.to_numpy())
    rows = []
    for gene in genes:
        column = common.get_loc(gene)
        for representation, before, after in (
            ("deposited_absolute", primary.loc[gene].to_numpy(), external.loc[gene].to_numpy()),
            ("within_sample_percentile_rank", primary_ranks[:, column], external_ranks[:, column]),
        ):
            rows.append(
                dict(
                    gene=gene, representation=representation,
                    primary_samples=len(before), external_samples=len(after),
                    primary_probes=probe_counts["GSE44076"][gene],
                    external_probes=probe_counts["GSE41258"][gene],
                    primary_median=np.median(before), external_median=np.median(after),
                    external_mean_shift_in_primary_sd=(after.mean() - before.mean()) / before.std(),
                    primary_sha256=file_sha256(primary_path),
                    external_sha256=file_sha256(external_path),
                )
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-matrix", type=Path)
    parser.add_argument("--external-matrix", type=Path)
    parser.add_argument("--bootstrap-replay", action="store_true")
    parser.add_argument("--output", type=Path, default=ROOT / "results/ai_audit")
    args = parser.parse_args()
    if bool(args.primary_matrix) != bool(args.external_matrix):
        parser.error("Supply both series matrices for the platform measurement diagnostic")
    summary = verify(ROOT)
    predictions = pd.read_csv(ROOT / ARTIFACTS["predictions"])
    results, strata = replay(predictions)
    args.output.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output / "external_prediction_replay.csv", index=False)
    strata.to_csv(args.output / "external_score_quantiles.csv", index=False)
    if args.bootstrap_replay:
        config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
        reference = pd.read_csv(ROOT / INTERVALS).set_index(EVALUATION_KEYS).loc[[PRIMARY_KEY]]
        canonical = predictions.groupby(EVALUATION_KEYS).get_group(PRIMARY_KEY)
        intervals = pd.DataFrame(grouped_bootstrap_metrics(
            canonical, reference["metric"].tolist(),
            config["external_validation"]["bootstrap_iterations"], config["project"]["random_seed"],
        ))
        joined = intervals.merge(reference.reset_index(drop=True), on="metric", suffixes=("", "_archived"))
        for metric in ("ci_lower", "ci_upper"):
            joined[f"{metric}_difference"] = joined[metric] - joined[f"{metric}_archived"]
        joined.to_csv(args.output / "external_rank_cluster_bootstrap_replay.csv", index=False)
    if args.primary_matrix:
        gene_shift(
            args.primary_matrix, args.external_matrix, summary["transport_signature_genes"]
        ).to_csv(args.output / "deposited_platform_gene_shift.csv", index=False)


if __name__ == "__main__":
    main()
