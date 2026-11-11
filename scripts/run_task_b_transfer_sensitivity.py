from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_geo_metadata import samples
from scripts.run_compact_panel_analysis import fit_panel_model
from scripts.run_development_audit import positive_scores
from src.data.cross_platform import within_sample_percentile_rank
from src.data.expression import aggregate_probe_expression, fetch_platform_mapping, read_geo_series_matrix
from src.data.provenance import artifact_sha256, file_sha256
from src.modeling.evaluation import classification_metrics, grouped_bootstrap_metrics


def early_onset_metadata(path: Path) -> pd.DataFrame:
    rows = []
    for sample_id, record in samples(path).items():
        title = record["title"].lower().replace("’", "'")
        source = record["source_name_ch1"].lower()
        donor = re.match(r"^(nm\d+|\d+m)\s", title)
        if donor is None or record["platform_id"] != "GPL570" or record["tissue"].lower() != "mucosa":
            raise ValueError(f"Unresolved GSE4107 source record {sample_id}")
        if "healthy control" in title and "colonoscopic" in source:
            label = "healthy"
        elif "patient's mucosa" in title and ("surgery" in source or "sugery" in source):
            label = "adjacent_normal"
        else:
            raise ValueError(f"Title and acquisition record disagree for {sample_id}")
        rows.append(dict(sample_id=sample_id, patient_id=donor[1], y_true=label))
    cohort = pd.DataFrame(rows)
    if cohort["patient_id"].duplicated().any() or cohort["y_true"].value_counts().to_dict() != {
        "adjacent_normal": 12, "healthy": 10,
    }:
        raise ValueError("GSE4107 cohort differs from the frozen eligibility definition")
    return cohort


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-matrix", required=True, type=Path)
    parser.add_argument("--external-matrix", required=True, type=Path)
    parser.add_argument("--external-annotation", required=True, type=Path)
    parser.add_argument("--external-geo-text", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "results/ai_audit")
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    config_path = ROOT / "config/analysis.yaml"
    card_path = ROOT / "models/task_b_final_model_card.json"
    primary_map = ROOT / "data/metadata/GSE44076_probe_gene_mapping_raw_cel_rma.csv"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    card = json.loads(card_path.read_text(encoding="utf-8"))
    with TemporaryDirectory() as temp:
        external_mapping = fetch_platform_mapping(
            args.external_annotation.resolve().as_uri(),
            "ID", "Gene symbol", "Gene ID", Path(temp) / "mapping.csv",
        )
    primary = aggregate_probe_expression(
        read_geo_series_matrix(args.primary_matrix),
        pd.read_csv(primary_map).query("mapping_status == 'mapped_unique'"),
    )
    primary_metadata_path = ROOT / "data/metadata/gse44076_samples.csv"
    metadata = pd.read_csv(primary_metadata_path)
    training = metadata[
        metadata["inclusion_status"].eq("included")
        & metadata["tissue_class"].isin(["healthy", "adjacent_normal"])
    ]
    common = primary.index.intersection(external_mapping["gene_symbol"]).tolist()
    genes = [gene for gene in card["feature_genes"] if gene in common]
    if len(genes) < 2:
        raise ValueError("GSE4107 annotation covers fewer than two prespecified Task B genes")
    columns = [common.index(gene) for gene in genes]
    values = primary.loc[common, training["geo_accession"]].T.to_numpy()
    seed = config["project"]["random_seed"]
    fitted = fit_panel_model(
        within_sample_percentile_rank(values)[:, columns],
        training["tissue_class"].to_numpy(), card["hyperparameters"], seed, config["modeling"],
    )
    threshold = 0.50
    lock = {
        "role": "external_cohort_association_sensitivity",
        "training_accession": "GSE44076",
        "evaluation_accession": "GSE4107",
        "expression_provenance": "geo_deposited_series_matrix",
        "representation": "within_sample_percentile_rank",
        "genes": genes,
        "missing_genes": [gene for gene in card["feature_genes"] if gene not in genes],
        "universe_size": len(common),
        "threshold": threshold,
        "positive_class": "adjacent_normal",
        "external_labels_used_for_fitting_or_threshold_selection": False,
        "seed": seed,
        "primary_matrix_sha256": file_sha256(args.primary_matrix),
        "external_annotation_sha256": file_sha256(args.external_annotation),
        "task_b_card_sha256": artifact_sha256(card_path),
        "development_metadata_sha256": artifact_sha256(primary_metadata_path),
        "primary_annotation_sha256": artifact_sha256(primary_map),
        "configuration_sha256": artifact_sha256(config_path),
    }
    (output / "gse4107_transfer_lock.json").write_text(
        json.dumps(lock, indent=2) + "\n", encoding="utf-8"
    )
    # Open external sample outcomes only after the model and threshold are fixed.
    cohort = early_onset_metadata(args.external_geo_text)
    external = aggregate_probe_expression(read_geo_series_matrix(args.external_matrix), external_mapping)
    values = external.loc[common, cohort["sample_id"]].T.to_numpy()
    scores = positive_scores(fitted, within_sample_percentile_rank(values)[:, columns], "adjacent_normal")
    cohort["probability_adjacent_normal"] = scores
    cohort["probability_healthy"] = 1 - scores
    cohort["y_pred"] = np.where(scores >= threshold, "adjacent_normal", "healthy")
    observed = cohort["y_true"].to_numpy()
    metrics = classification_metrics(
        observed, cohort["y_pred"].to_numpy(), np.column_stack([1 - scores, scores]),
        np.array(["healthy", "adjacent_normal"]),
    )
    intervals = grouped_bootstrap_metrics(
        cohort, ["roc_auc", "f1_macro", "balanced_accuracy", "brier_score", "log_loss"], 1000, seed
    )
    pd.DataFrame(intervals).to_csv(output / "gse4107_transfer_bootstrap.csv", index=False)
    pd.DataFrame([{
        "accession": "GSE4107", "patients": len(cohort), "threshold": threshold,
        "external_matrix_sha256": file_sha256(args.external_matrix),
        "external_metadata_sha256": file_sha256(args.external_geo_text), **metrics,
    }]).to_csv(output / "gse4107_transfer_metrics.csv", index=False)


if __name__ == "__main__":
    main()
