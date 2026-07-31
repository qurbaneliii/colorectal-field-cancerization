from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.biology import (
    classify_trajectories,
    exploratory_differential_expression,
    field_candidates,
)
from src.data.cross_platform import common_gene_symbols
from src.data.expression import (
    aggregate_probe_expression,
    fetch_platform_mapping,
    read_geo_series_matrix,
)
from src.data.metadata import discover_accession_files
from src.data.provenance import expression_path, sample_metadata_path
from src.data.validation import validate_expression
from src.visualization.publication_figures import (
    correlation_heatmap,
    expression_boxplot,
    pca_plot,
    top_variable_heatmap,
    volcano_plot,
)


def require_one(paths: list[Path], label: str) -> Path:
    if len(paths) != 1:
        raise FileNotFoundError(f"Expected one {label}, observed {len(paths)}: {paths}")
    return paths[0]


def save_expression(frame: pd.DataFrame, path: Path, id_name: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.reset_index(names=id_name).to_parquet(path, index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures-only", action="store_true")
    args = parser.parse_args()
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    raw_root = ROOT / paths["raw"]
    metadata_root = ROOT / paths["metadata"]
    processed_root = ROOT / paths["processed"]
    figures_root = ROOT / paths["figures"]
    tables_root = ROOT / paths["tables"]
    de_root = ROOT / paths["de"]
    for directory in (metadata_root, processed_root, figures_root, tables_root, de_root):
        directory.mkdir(parents=True, exist_ok=True)

    outputs: dict[str, pd.DataFrame] = {}
    for accession in ("GSE44076", "GSE41258"):
        metadata = pd.read_csv(metadata_root / f"{accession.lower()}_samples.csv", dtype={"patient_id": str})
        matrix_path = require_one(
            discover_accession_files(raw_root, accession)["series_matrix"],
            f"{accession} series matrix",
        )
        provenance = "geo_deposited_series_matrix"
        gene_path = expression_path(processed_root, accession, "gene", provenance)
        probe_path = expression_path(processed_root, accession, "probe", provenance)
        mapping_path = metadata_root / f"{accession}_probe_gene_mapping.csv"
        if args.figures_only and gene_path.exists():
            gene_table = pd.read_parquet(gene_path)
            gene = gene_table.set_index("gene_symbol")
        else:
            probe = read_geo_series_matrix(matrix_path)
            ann = config["annotation"][config["expected"][accession]["platform"]]
            mapping = (
                pd.read_csv(mapping_path, dtype=str)
                if mapping_path.exists()
                else fetch_platform_mapping(
                    ann["url"],
                    ann["id_column"],
                    ann["symbol_column"],
                    ann["entrez_column"],
                    mapping_path,
                )
            )
            gene = aggregate_probe_expression(probe, mapping)
            save_expression(probe, probe_path, "probe_id")
            save_expression(gene, gene_path, "gene_symbol")
        metadata = metadata.copy()
        metadata["expression_provenance"] = provenance
        metadata.to_csv(sample_metadata_path(processed_root, accession, provenance), index=False)
        validate_expression(gene.reset_index(), metadata, "gene_symbol")
        outputs[accession] = gene
        expression_boxplot(
            gene,
            figures_root / f"{accession}_deposited_normalized_expression_boxplot",
            f"{accession}: deposited normalized gene expression",
        )
        pca = pca_plot(
            gene,
            metadata,
            figures_root / f"{accession}_primary_pca",
            f"{accession}: PCA of deposited normalized expression",
        )
        pca.to_csv(tables_root / f"{accession}_pca_scores.csv", index=False)
        correlation_heatmap(
            gene,
            figures_root / f"{accession}_sample_correlation_heatmap",
            f"{accession}: sample correlation (top variable genes)",
        )
        if accession == "GSE44076":
            top_variable_heatmap(
                gene, metadata, figures_root / "GSE44076_top_variable_gene_heatmap"
            )
            med = np.median(gene.to_numpy(), axis=0)
            mad = np.median(np.abs(med - np.median(med)))
            outliers = pd.DataFrame(
                {
                    "geo_accession": gene.columns,
                    "sample_median": med,
                    "robust_z": 0.6745 * (med - np.median(med)) / (mad if mad else 1.0),
                }
            )
            outliers["outlier_candidate"] = outliers["robust_z"].abs().gt(5.0)
            outliers["automatic_exclusion"] = False
            outliers.to_csv(tables_root / "GSE44076_outlier_candidates.csv", index=False)
            pd.DataFrame(
                {
                    "missing_values": [int(gene.isna().sum().sum())],
                    "infinite_values": [int(np.isinf(gene.to_numpy()).sum())],
                    "genes": [len(gene)],
                    "samples": [gene.shape[1]],
                }
            ).to_csv(tables_root / "GSE44076_missingness_report.csv", index=False)

    common = common_gene_symbols(
        outputs["GSE44076"].reset_index(), outputs["GSE41258"].reset_index()
    )
    (processed_root / "common_genes_GSE44076_GSE41258_geo_deposited_series_matrix.txt").write_text(
        "\n".join(common) + "\n", encoding="utf-8"
    )

    primary_meta = pd.read_csv(metadata_root / "gse44076_samples.csv", dtype={"patient_id": str})
    primary_gene = outputs["GSE44076"]
    de = exploratory_differential_expression(primary_gene, primary_meta)
    de.to_csv(de_root / "processed_matrix_exploratory_all_comparisons.csv", index=False)
    for comparison, subset in de.groupby("comparison"):
        subset.to_csv(de_root / f"processed_matrix_{comparison}.csv", index=False)
    biology = config["biology"]
    trajectories = classify_trajectories(
        primary_gene, primary_meta, biology["trajectory_tolerance"]
    )
    trajectories.to_csv(tables_root / "expression_trajectory_genes.csv", index=False)
    candidates = field_candidates(
        primary_gene,
        primary_meta,
        de,
        trajectories,
        biology["fdr_threshold"],
        biology["log2fc_threshold"],
    )
    candidates.to_csv(tables_root / "field_cancerization_candidates.csv", index=False)
    adjacent_de = de[de["comparison"].eq("adjacent_normal_vs_healthy")]
    volcano_plot(
        adjacent_de,
        figures_root / "GSE44076_adjacent_vs_healthy_volcano_exploratory",
        "Adjacent-normal versus healthy (processed-matrix sensitivity analysis)",
        biology["fdr_threshold"],
        biology["log2fc_threshold"],
    )

    # Quantify where adjacent samples lie relative to healthy and tumor centroids.
    meta = primary_meta.set_index("geo_accession").loc[primary_gene.columns]
    variable = primary_gene.var(axis=1).nlargest(min(5000, len(primary_gene))).index
    x = primary_gene.loc[variable].T
    healthy_centroid = x.loc[meta["tissue_class"].eq("healthy")].mean(axis=0).to_numpy()
    tumor_centroid = x.loc[meta["tissue_class"].eq("tumor")].mean(axis=0).to_numpy()
    adjacent_ids = meta.index[meta["tissue_class"].eq("adjacent_normal")]
    adjacent_x = x.loc[adjacent_ids].to_numpy()
    distances = pd.DataFrame(
        {
            "geo_accession": adjacent_ids,
            "patient_id": meta.loc[adjacent_ids, "patient_id"].to_numpy(),
            "distance_to_healthy_centroid": cdist(
                adjacent_x, healthy_centroid.reshape(1, -1), metric="euclidean"
            ).ravel(),
            "distance_to_tumor_centroid": cdist(
                adjacent_x, tumor_centroid.reshape(1, -1), metric="euclidean"
            ).ravel(),
        }
    )
    distances["closer_centroid"] = np.where(
        distances["distance_to_healthy_centroid"] < distances["distance_to_tumor_centroid"],
        "healthy",
        "tumor",
    )
    distances.to_csv(tables_root / "adjacent_centroid_distances.csv", index=False)

    summary = {
        "GSE44076": {
            "genes": len(outputs["GSE44076"]),
            "samples": outputs["GSE44076"].shape[1],
        },
        "GSE41258": {
            "genes": len(outputs["GSE41258"]),
            "samples": outputs["GSE41258"].shape[1],
        },
        "common_genes": len(common),
        "field_candidates_processed_matrix_sensitivity": len(candidates),
        "provenance": "geo_deposited_series_matrix",
        "scale_note": {
            "GSE44076": "observed values are compatible with log2 expression",
            "GSE41258": "deposited values are not log2; no cross-sample scale claim is made",
        },
    }
    (ROOT / paths["reports"] / "processed_matrix_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
