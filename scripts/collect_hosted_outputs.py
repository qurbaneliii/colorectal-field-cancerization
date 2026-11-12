from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


MANIFEST = Path("reports/result_artifact_manifest.csv")
MAIN_TABLES = (
    "table_1_cohort_characteristics.csv",
    "table_2_raw_cel_qc_summary.csv",
    "table_3_unadjusted_and_adjusted_de_summary.csv",
    "table_4_high_confidence_field_signature.csv",
    "table_5_enrichment_summary.csv",
    "table_6_internal_model_comparison.csv",
    "table_7_task_b_compact_signature.csv",
    "table_8_task_c_compact_signature.csv",
    "table_9_external_validation.csv",
    "table_10_sensitivity_analyses.csv",
)
PUBLICATION_FIGURES = (
    "study_objectives_and_dataset_roles",
    "sample_count_flowchart",
    "raw_array_qc_summary",
    "unadjusted_adjacent_vs_healthy_volcano",
    "adjusted_adjacent_vs_healthy_volcano",
    "unadjusted_vs_adjusted_field_logfc",
    "paired_tumor_vs_adjacent_volcano",
    "final_field_signature_heatmap",
    "top_field_gene_group_distributions",
    "tissue_composition_by_group",
    "field_gene_evidence_tier_summary",
    "high_confidence_field_enrichment",
    "task_a_out_of_fold_confusion_matrix",
    "task_b_roc_pr_calibration",
    "task_b_panel_size_performance",
    "task_c_roc_pr_calibration",
    "task_c_panel_size_performance",
    "stable_feature_coefficient_plot_task_c",
    "task_c_oof_threshold_curve",
    "external_validation_roc_pr",
    "external_calibration",
    "external_confusion_matrix",
    "external_patient_structure_sensitivity",
    "signature_intersection",
    "biological_predictive_gene_intersection",
    "adjacent_vs_healthy_u1_ranked_go_bp",
    "tumor_vs_adjacent_paired_ranked_go_bp",
)
REQUIRED_OUTPUTS = {
    MANIFEST.as_posix(),
    "manuscript/manuscript.md",
    "reports/external_validation_report.json",
    "reports/publication_readiness_final.md",
    "reports/tier_specific_enrichment_report.md",
    "models/task_b_final_model_card.json",
    "models/task_c_primary_full_signature_model_card.json",
    "models/task_c_cross_platform_transport_model_card.json",
    "results/tables/supplementary_enrichment_results.csv",
    "results/tables/enrichment_analysis_set_audit.csv",
    "results/enrichment/adjacent_vs_healthy_u1_ranked_go_bp.csv",
    "results/enrichment/tumor_vs_adjacent_paired_ranked_go_bp.csv",
    "results/metrics/nested_cv_predictions.csv",
    "results/metrics/fold_assignments.csv",
    "results/metrics/external_validation_predictions.csv",
    "results/figures/internal_validation_figure_data.json",
    *(f"results/tables/{name}" for name in MAIN_TABLES),
    *(f"results/figures/{stem}.{ext}" for stem in PUBLICATION_FIGURES for ext in ("png", "pdf", "svg")),
}


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify and retain fresh hosted raw-CEL outputs.")
    parser.add_argument("--worktree", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--archive-identity", choices=["historical", "current-geo"], required=True)
    args = parser.parse_args()
    root = args.worktree
    started = (args.evidence / "pipeline-start").stat().st_mtime_ns
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()

    with (root / MANIFEST).open(encoding="utf-8", newline="") as handle:
        entries = {row["path"]: row for row in csv.DictReader(handle)}
    if {row["inventory_head_commit"] for row in entries.values()} != {revision}:
        raise ValueError("Publication inventory belongs to a different analysis commit")
    with (args.evidence / "stages.csv").open(encoding="utf-8", newline="") as handle:
        stages = list(csv.DictReader(handle))
    if len(stages) != 20 or any(row["status"] != "PASS" or row["exit_code"] != "0" for row in stages):
        raise ValueError("Publication evidence requires all 20 scientific stages to pass")

    workflow_revision = os.environ.get("GITHUB_SHA")
    stale = sorted(
        name for name in REQUIRED_OUTPUTS
        if not (root / name).is_file() or (root / name).stat().st_mtime_ns < started
    )
    if stale and workflow_revision in (None, revision):
        raise ValueError("Required publication outputs were not regenerated: " + ", ".join(stale))
    selected = REQUIRED_OUTPUTS.difference(stale)
    for pattern in (
        "manuscript/*.md",
        "models/*.json",
        "reports/*.md",
        "reports/*.json",
        "results/enrichment/*.csv",
        "results/metrics/*.csv",
        "results/figures/*.json",
    ):
        selected.update(
            path.relative_to(root).as_posix() for path in root.glob(pattern)
            if path.is_file() and path.stat().st_mtime_ns >= started
        )
    artifacts = []
    for name in sorted(selected):
        source = root / name
        if not source.is_file() or source.stat().st_mtime_ns < started or source.stat().st_size == 0:
            raise ValueError(f"Expected a fresh, nonempty publication output: {name}")
        digest = sha256(source)
        size = source.stat().st_size
        if name != MANIFEST.as_posix():
            entry = entries.get(name)
            if entry is None or entry["sha256"] != digest or int(entry["size_bytes"]) != size:
                raise ValueError(f"Publication inventory does not match generated bytes: {name}")
        target = args.output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        artifacts.append({"path": name, "bytes": size, "sha256": digest})

    for name in ("inputs.json", "stages.csv", "full-data-test.log"):
        shutil.copyfile(args.evidence / name, args.output / name)
    evidence = {
        "analysis_commit": revision,
        "workflow_commit": workflow_revision,
        "workflow_run_id": os.environ.get("GITHUB_RUN_ID"),
        "workflow_run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        "archive_identity": args.archive_identity,
        "inputs_sha256": sha256(args.output / "inputs.json"),
        "stages_sha256": sha256(args.output / "stages.csv"),
        "test_log_sha256": sha256(args.output / "full-data-test.log"),
        "missing_or_stale_required_outputs": stale,
        "artifacts": artifacts,
    }
    (args.output / "generated-output-manifest.json").write_text(
        json.dumps(evidence, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Verified and collected {len(artifacts)} fresh outputs from {revision}")


if __name__ == "__main__":
    main()
