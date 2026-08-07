from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results/tables"
METRICS = ROOT / "results/metrics"


def read(name: str, root: Path = TABLES) -> pd.DataFrame:
    return pd.read_csv(root / name)


def cohort_table() -> pd.DataFrame:
    frames = []
    for accession, filename in (
        ("GSE44076", "gse44076_samples.csv"),
        ("GSE41258", "gse41258_samples.csv"),
    ):
        metadata = pd.read_csv(ROOT / "data/metadata" / filename)
        included = metadata[metadata["inclusion_status"].eq("included")]
        for tissue, group in included.groupby("tissue_class"):
            frames.append(
                {
                    "accession": accession,
                    "tissue_class": tissue,
                    "arrays": len(group),
                    "unique_patients_or_donors": group["donor_or_patient_group"].nunique(),
                    "age_available": group["age"].notna().sum(),
                    "age_mean": group["age"].mean(),
                    "age_sd": group["age"].std(),
                    "female": int(group["sex"].eq("Female").sum()),
                    "male": int(group["sex"].eq("Male").sum()),
                    "left": int(group["location"].eq("Left").sum()),
                    "right": int(group["location"].eq("Right").sum()),
                }
            )
    return pd.DataFrame(frames)


def qc_table() -> pd.DataFrame:
    rows = []
    for accession in ("GSE44076", "GSE41258"):
        custom = read(f"{accession}_raw_cel_qc_diagnostics.csv")
        aqm = read(f"{accession}_array_quality_metrics_summary.csv")
        row: dict[str, object] = {
            "accession": accession,
            "arrays": len(custom),
            "non_finite_values": int(custom["non_finite_values"].sum()),
            "one_custom_failure_metric": int(
                custom["independent_failure_metrics"].eq(1).sum()
            ),
            "two_or_more_custom_failure_metrics": int(
                custom["independent_failure_metrics"].ge(2).sum()
            ),
            "aqm_flagged_any": int(aqm["aqm_flagged_any"].sum()),
            "aqm_two_or_more_flags": int(aqm["aqm_flag_count"].ge(2).sum()),
            "automatic_exclusions": int(custom["automatic_exclusion"].sum()),
            "aqm_report_files": len(
                list((ROOT / "results/qc" / f"{accession}_array_quality_metrics").glob("*"))
            ),
        }
        nuse_path = TABLES / f"{accession}_nuse_summary.csv"
        if nuse_path.exists():
            nuse = pd.read_csv(nuse_path)
            row.update(
                {
                    "nuse_method": nuse["method"].iloc[0],
                    "nuse_median_of_medians": nuse["median"].median(),
                    "nuse_minimum_median": nuse["median"].min(),
                    "nuse_maximum_median": nuse["median"].max(),
                }
            )
        else:
            row.update(
                {
                    "nuse_method": "not technically valid for this oligo platform",
                    "nuse_median_of_medians": np.nan,
                    "nuse_minimum_median": np.nan,
                    "nuse_maximum_median": np.nan,
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def de_table() -> pd.DataFrame:
    de = read("supplementary_full_raw_cel_de_results.csv")
    audit = read("covariate_design_audit.csv").set_index(["model", "contrast"])
    rows = []
    for (model, comparison), group in de.groupby(["model", "comparison"]):
        contrast = group["contrast"].iloc[0]
        design = audit.loc[(model, contrast)]
        significant = group["adjusted_p_value"].lt(0.05) & group[
            "log2_fold_change"
        ].abs().ge(0.5)
        rows.append(
            {
                "model": model,
                "comparison": comparison,
                "tested_genes": len(group),
                "fdr_lt_0_05_abs_log2fc_ge_0_5": int(significant.sum()),
                "up": int((significant & group["log2_fold_change"].gt(0)).sum()),
                "down": int((significant & group["log2_fold_change"].lt(0)).sum()),
                "design_rank": int(design["design_rank"]),
                "design_columns": int(design["design_columns"]),
                "full_rank": bool(design["full_rank"]),
                "analysis_provenance": group["analysis_provenance"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def sensitivity_table() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    task_b = read("task_b_confounding_sensitivity.csv")
    for _, row in task_b.iterrows():
        rows.append(
            {
                "domain": "Task B confounding",
                "analysis": row["scenario"],
                "key_metric": "macro_f1",
                "value": row["f1_macro"],
                "secondary_metric": f"ROC-AUC={row['roc_auc']:.4f}",
                "interpretation": "locked panel; grouped repeated OOF",
            }
        )
    high = read("final_high_confidence_field_signature.csv")
    rows.extend(
        [
            {
                "domain": "Tissue composition",
                "analysis": "high-confidence genes retaining U3 support",
                "key_metric": "gene_count",
                "value": int(high["composition_robust"].sum()),
                "secondary_metric": f"of {len(high)}",
                "interpretation": "MCP-counter PC-adjusted limma sensitivity",
            },
            {
                "domain": "Preanalytical stress",
                "analysis": "high-confidence genes retained after curated stress removal",
                "key_metric": "gene_count",
                "value": len(read("high_confidence_field_signature_without_stress_genes.csv")),
                "secondary_metric": f"of {len(high)}",
                "interpretation": "curated immediate-early/stress sensitivity",
            },
        ]
    )
    qc = read("qc_exclusion_sensitivity.csv")
    for _, row in qc.iterrows():
        rows.append(
            {
                "domain": "Borderline-sample QC",
                "analysis": row["analysis_component"],
                "key_metric": "sensitivity_value",
                "value": row["sensitivity_value"],
                "secondary_metric": f"primary={row['primary_value']}; {row['comparison']}",
                "interpretation": "primary included; sensitivity excluded",
            }
        )
    concordance = read("raw_vs_deposited_de_concordance.csv")
    for _, row in concordance.iterrows():
        rows.append(
            {
                "domain": "Raw versus deposited",
                "analysis": row.get("comparison", "concordance"),
                "key_metric": "spearman_rank_correlation",
                "value": row.get("spearman_rank_correlation", np.nan),
                "secondary_metric": "independently processed representations",
                "interpretation": "sensitivity; raw-CEL RMA is primary",
            }
        )
    return pd.DataFrame(rows)


def finalize_qc_exclusion_sensitivity() -> None:
    qc = read("qc_exclusion_sensitivity.csv")
    task_b = read("task_b_confounding_sensitivity.csv").set_index("scenario")
    primary = task_b.loc["primary_fixed_panel"]
    excluded = task_b.loc["borderline_sample_excluded_fixed_panel"]
    genes = read("qc_exclusion_model_gene_comparison.csv")
    shared = genes["primary_panel"] & genes["qc_excluded_reselected_panel"]
    enrichment = read("enrichment_top_terms_by_analysis_set.csv")
    primary_terms = set(
        enrichment.loc[
            enrichment["analysis_set"].eq("adjacent_vs_healthy_u1"), "Description"
        ].dropna().head(20)
    )
    sensitivity_terms = set(
        enrichment.loc[
            enrichment["analysis_set"].isin(
                ["qc_excluded_field_up", "qc_excluded_field_down"]
            ),
            "Description",
        ].dropna().head(20)
    )
    additions = pd.DataFrame(
        [
            {
                "analysis_component": "task_b_aggregated_oof_macro_f1",
                "primary_value": primary["f1_macro"],
                "sensitivity_value": excluded["f1_macro"],
                "comparison": "five-gene panel; repeated grouped OOF",
            },
            {
                "analysis_component": "task_b_aggregated_oof_balanced_accuracy",
                "primary_value": primary["balanced_accuracy"],
                "sensitivity_value": excluded["balanced_accuracy"],
                "comparison": "five-gene panel; repeated grouped OOF",
            },
            {
                "analysis_component": "task_b_reselected_panel_overlap",
                "primary_value": int(genes["primary_panel"].sum()),
                "sensitivity_value": int(shared.sum()),
                "comparison": (
                    "shared=" + ";".join(genes.loc[shared, "gene_symbol"].astype(str))
                    + "; sensitivity="
                    + ";".join(
                        genes.loc[genes["qc_excluded_reselected_panel"], "gene_symbol"].astype(str)
                    )
                ),
            },
            {
                "analysis_component": "top_ranked_pathway_term_overlap",
                "primary_value": len(primary_terms),
                "sensitivity_value": len(primary_terms & sensitivity_terms),
                "comparison": (
                    "primary ranked U1 versus QC-excluded over-representation; "
                    + ";".join(sorted(primary_terms & sensitivity_terms))
                ),
            },
        ]
    )
    qc = pd.concat(
        [
            qc[~qc["analysis_component"].isin(additions["analysis_component"])],
            additions,
        ],
        ignore_index=True,
    )
    qc.to_csv(TABLES / "qc_exclusion_sensitivity.csv", index=False)
    borderline = qc.loc[
        qc["analysis_component"].eq("borderline_samples"), "comparison"
    ].iloc[0]
    report = f"""# QC exclusion sensitivity

No GSE44076 array met the prespecified custom-QC exclusion rule of at least two
independent severe technical failures. `{borderline}` had one severe custom
IQR flag and was retained in the primary analysis; the sensitivity analysis
excluded it. arrayQualityMetrics flags remained review signals and did not
change the exclusion decision. PCA separation was never an exclusion rule.

{qc.to_markdown(index=False)}

The sensitivity reproduces the age/sex-adjusted differential-expression model,
the locked Task B classifier under repeated donor/patient-grouped out-of-fold
evaluation, and a separate full-development grouped-inner-CV panel selection.
Pathway overlap compares up to the top 20 exact descriptions available from
ranked U1 GO-BP analysis with the union of the QC-excluded up/down
over-representation results;
the different enrichment estimands are retained explicitly rather than treated
as interchangeable. Complete gene membership is in
`results/tables/qc_exclusion_model_gene_comparison.csv`.
"""
    (ROOT / "reports/qc_exclusion_sensitivity.md").write_text(
        report, encoding="utf-8"
    )


def main() -> None:
    os.chdir(ROOT)
    TABLES.mkdir(parents=True, exist_ok=True)
    cohort_table().to_csv(TABLES / "table_1_cohort_characteristics.csv", index=False)
    qc_table().to_csv(TABLES / "table_2_raw_cel_qc_summary.csv", index=False)
    de_table().to_csv(
        TABLES / "table_3_unadjusted_and_adjusted_de_summary.csv", index=False
    )
    read("final_high_confidence_field_signature.csv").to_csv(
        TABLES / "table_4_high_confidence_field_signature.csv", index=False
    )
    enrichment = read("enrichment_top_terms_by_analysis_set.csv")
    enrichment.to_csv(TABLES / "table_5_enrichment_summary.csv", index=False)
    high_enrichment = enrichment[
        enrichment["analysis_set"].str.startswith("high_confidence")
        | enrichment["analysis_set"].str.startswith("composition_")
    ]
    high_enrichment.to_csv(TABLES / "high_confidence_field_enrichment.csv", index=False)
    enrichment[
        enrichment["analysis_set"].isin(
            ["task_b_locked_panel", "task_c_locked_panel"]
        )
    ].to_csv(TABLES / "signature_enrichment_summary.csv", index=False)
    read("table_3_model_comparison.csv").to_csv(
        TABLES / "table_6_internal_model_comparison.csv", index=False
    )
    read("task_b_final_signature.csv").to_csv(
        TABLES / "table_7_task_b_compact_signature.csv", index=False
    )
    read("final_compact_signature.csv").to_csv(
        TABLES / "table_8_task_c_compact_signature.csv", index=False
    )
    pd.read_csv(METRICS / "external_patient_tissue_metrics.csv").to_csv(
        TABLES / "table_9_external_validation.csv", index=False
    )
    finalize_qc_exclusion_sensitivity()
    sensitivity_table().to_csv(TABLES / "table_10_sensitivity_analyses.csv", index=False)

    required_supplements = [
        "data/metadata/gse44076_samples.csv",
        "data/metadata/gse41258_samples.csv",
        "data/metadata/sample_exclusion_log.csv",
        "data/metadata/GSE44076_probe_gene_mapping_raw_cel_rma.csv",
        "data/metadata/GSE41258_probe_gene_mapping_raw_cel_rma.csv",
        "results/tables/supplementary_full_raw_cel_de_results.csv",
        "results/tables/field_gene_evidence_tiers.csv",
        "results/tables/supplementary_enrichment_results.csv",
        "results/metrics/fold_assignments.csv",
        "results/metrics/nested_cv_predictions.csv",
        "results/metrics/nested_cv_hyperparameters.csv",
        "results/tables/feature_stability.csv",
        "results/metrics/nested_panel_selection_decisions.csv",
        "results/metrics/group_preserving_permutation_tests.csv",
        "results/metrics/external_validation_predictions.csv",
        "results/metrics/external_patient_cluster_bootstrap_ci.csv",
    ]
    manifest_rows = []
    for relative in required_supplements:
        path = ROOT / relative
        digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""
        manifest_rows.append(
            {
                "artifact": relative.replace("\\", "/"),
                "exists": path.exists(),
                "bytes": path.stat().st_size if path.exists() else 0,
                "sha256": digest,
            }
        )
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(ROOT / "reports/result_artifact_manifest.csv", index=False)
    if not manifest["exists"].all() or not manifest["bytes"].gt(0).all():
        raise AssertionError("A required supplementary artifact is missing or empty")
    print(
        json.dumps(
            {
                "publication_tables": 10,
                "manifest_artifacts": len(manifest),
                "all_manifest_artifacts_nonempty": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
