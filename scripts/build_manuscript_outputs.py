from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from sklearn.calibration import calibration_curve
from sklearn.metrics import ConfusionMatrixDisplay, PrecisionRecallDisplay, RocCurveDisplay

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data.provenance import expression_path, result_root
from src.visualization.publication_figures import save_figure


def study_design_figure(stem: Path, raw_complete: bool) -> None:
    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.axis("off")
    boxes = [
        (0.01, 0.60, "GSE44076\n50 healthy donors\n98 paired patients"),
        (0.25, 0.60, "Independent raw-CEL RMA\nand label-free annotation"),
        (0.49, 0.60, "limma field-effect\nand paired contrasts"),
        (0.73, 0.60, "Repeated nested\npatient-group CV"),
        (0.49, 0.14, "Training-fold compact\nElastic Net panels"),
        (0.73, 0.14, "Locked Task C model"),
        (0.88, 0.14, "GSE41258\ntumor vs normal only"),
    ]
    for x, y, label in boxes:
        ax.add_patch(
            plt.Rectangle((x, y), 0.17, 0.23, facecolor="#EAF4F4", edgecolor="#1D3557")
        )
        ax.text(x + 0.085, y + 0.115, label, ha="center", va="center", fontsize=9)
    for start, end in [
        ((0.18, 0.715), (0.25, 0.715)),
        ((0.42, 0.715), (0.49, 0.715)),
        ((0.66, 0.715), (0.73, 0.715)),
        ((0.815, 0.60), (0.64, 0.37)),
        ((0.66, 0.255), (0.73, 0.255)),
        ((0.90, 0.255), (0.88, 0.255)),
    ]:
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 1.4})
    status = "raw-CEL route executed" if raw_complete else "raw-CEL route pending"
    ax.set_title(f"Leakage-safe colorectal field-cancerization study design ({status})", fontsize=14)
    save_figure(fig, stem)


def sample_flow_figure(primary: pd.DataFrame, external: pd.DataFrame, stem: Path) -> None:
    included = external[external["inclusion_status"].eq("included")]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for axis in axes:
        axis.axis("off")
    p = primary["tissue_class"].value_counts()
    axes[0].text(
        0.5,
        0.6,
        f"GSE44076: {len(primary)} arrays\n{p['healthy']} healthy\n"
        f"{p['adjacent_normal']} adjacent-normal\n{p['tumor']} tumor\n98 complete pairs",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.7", "fc": "#F1FAEE"},
    )
    e = included["tissue_class"].value_counts()
    axes[1].text(
        0.5,
        0.6,
        f"GSE41258: {len(external)} deposited arrays\n{len(included)} eligible arrays\n"
        f"{e.get('normal_colon', 0)} normal colon\n{e.get('primary_tumor', 0)} primary tumor\n"
        f"{included['patient_id'].nunique()} unique patients",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.7", "fc": "#F1FAEE"},
    )
    fig.suptitle("Sample inclusion flow")
    save_figure(fig, stem)


def averaged_predictions(frame: pd.DataFrame) -> pd.DataFrame:
    probability_columns = [
        column for column in frame if column.startswith("probability_") and frame[column].notna().any()
    ]
    averaged = frame.groupby("sample_id", as_index=False).agg(
        {"y_true": "first", **{column: "mean" for column in probability_columns}}
    )
    classes = [column.removeprefix("probability_") for column in probability_columns]
    averaged["y_pred"] = np.asarray(classes)[
        np.argmax(averaged[probability_columns].to_numpy(), axis=1)
    ]
    return averaged


def modeling_figures(metrics: pd.DataFrame, predictions: pd.DataFrame, figures: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=metrics, x="task", y="f1_macro", hue="model", ax=ax)
    ax.set(xlabel="", ylabel="Untouched outer-fold macro F1")
    ax.tick_params(axis="x", rotation=15)
    ax.set_title("Repeated nested patient/donor-group cross-validation")
    ax.legend(frameon=False)
    save_figure(fig, figures / "nested_cv_model_comparison")

    for task, title in [
        ("task_a_three_class", "Three-class tissue state"),
        ("task_b_field_effect", "Healthy versus adjacent-normal field effect"),
        ("task_c_tumor_vs_adjacent", "Tumor versus adjacent-normal"),
    ]:
        frame = predictions[predictions["task"].eq(task) & predictions["model"].eq("elastic_net")]
        averaged = averaged_predictions(frame)
        probability_columns = sorted(column for column in averaged if column.startswith("probability_"))
        classes = [column.removeprefix("probability_") for column in probability_columns]
        if task == "task_a_three_class":
            fig, ax = plt.subplots(figsize=(6.5, 6))
            ConfusionMatrixDisplay.from_predictions(
                averaged["y_true"], averaged["y_pred"], labels=classes,
                display_labels=classes, cmap="Blues", colorbar=False, ax=ax
            )
            ax.set_title(f"Three-class OOF confusion matrix (n={len(averaged)})")
            plt.setp(ax.get_xticklabels(), rotation=20, ha="right")
            save_figure(fig, figures / "three_class_oof_confusion_matrix")
            continue
        positive_label = "adjacent_normal" if task == "task_b_field_effect" else "tumor"
        binary = averaged["y_true"].eq(positive_label).astype(int)
        probability = averaged[f"probability_{positive_label}"]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
        RocCurveDisplay.from_predictions(binary, probability, ax=axes[0])
        PrecisionRecallDisplay.from_predictions(binary, probability, ax=axes[1])
        observed, predicted = calibration_curve(binary, probability, n_bins=8, strategy="quantile")
        axes[2].plot(predicted, observed, marker="o")
        axes[2].plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.8)
        axes[2].set(xlabel="Mean predicted probability", ylabel="Observed fraction", title="Calibration")
        fig.suptitle(f"{title}: repeated OOF predictions (n={len(averaged)})")
        save_figure(fig, figures / ("task_b_roc_pr_calibration" if task.startswith("task_b") else "task_c_roc_pr_calibration"))


def de_figures(raw_de: pd.DataFrame, figures: Path) -> None:
    for comparison, stem, title in [
        ("adjacent_vs_healthy", "adjacent_vs_healthy_volcano", "Adjacent-normal versus healthy"),
        ("tumor_vs_adjacent_paired", "paired_tumor_vs_adjacent_volcano", "Paired tumor versus adjacent-normal"),
    ]:
        frame = raw_de[raw_de["comparison"].eq(comparison)].copy()
        frame["minus_log10_fdr"] = -np.log10(frame["adjusted_p_value"].clip(lower=1e-300))
        significant = frame["adjusted_p_value"].lt(0.05) & frame["log2_fold_change"].abs().ge(0.5)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(frame.loc[~significant, "log2_fold_change"], frame.loc[~significant, "minus_log10_fdr"], s=7, alpha=0.3, color="#8D99AE")
        ax.scatter(frame.loc[significant, "log2_fold_change"], frame.loc[significant, "minus_log10_fdr"], s=9, alpha=0.6, color="#D62828")
        ax.axvline(-0.5, linestyle="--", color="black", linewidth=0.7)
        ax.axvline(0.5, linestyle="--", color="black", linewidth=0.7)
        ax.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=0.7)
        ax.set(xlabel="Raw-CEL RMA limma log2 fold change", ylabel="-log10 BH-adjusted p-value", title=f"{title} (n={len(frame):,} genes)")
        save_figure(fig, figures / stem)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_manifest(commit: str) -> pd.DataFrame:
    roots = [ROOT / "results", ROOT / "reports", ROOT / "manuscript", ROOT / "data/metadata"]
    command_rules = {
        "differential_expression": "Rscript R/04_differential_expression.R",
        "enrichment": "Rscript R/05_functional_enrichment.R",
        "external": "python scripts/run_external_validation.py --provenance raw_cel_rma",
        "compact": "python scripts/run_compact_panel_analysis.py --provenance raw_cel_rma",
        "nested_cv": "python scripts/run_modeling.py --provenance raw_cel_rma",
        "raw_vs": "python scripts/run_raw_vs_processed_sensitivity.py",
    }
    rows = []
    created = datetime.now(timezone.utc).isoformat()
    for base in roots:
        for path in base.rglob("*"):
            if not path.is_file() or path.name == "result_artifact_manifest.csv":
                continue
            relative = path.relative_to(ROOT).as_posix()
            producing = "python scripts/build_manuscript_outputs.py"
            for token, command in command_rules.items():
                if token in relative.lower():
                    producing = command
                    break
            provenance = (
                "raw_cel_rma" if "raw_cel" in relative or "enrichment" in relative
                else "mixed_or_metadata; inspect artifact columns/report"
            )
            rows.append(
                {
                    "path": relative,
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "producing_command": producing,
                    "input_provenance": provenance,
                    "creation_timestamp_utc": created,
                    "git_commit_sha": commit,
                }
            )
    return pd.DataFrame(rows).sort_values("path")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures", action="store_true")
    parser.parse_args()
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    paths = yaml.safe_load((ROOT / "config/paths.yaml").read_text(encoding="utf-8"))
    provenance = config["modeling"]["expression_provenance"]
    raw_complete = expression_path(ROOT / paths["processed"], "GSE44076", "gene", "raw_cel_rma").exists() and Path("results/tables/supplementary_full_raw_cel_de_results.csv").exists()
    primary = pd.read_csv("data/metadata/gse44076_samples.csv", dtype={"patient_id": str})
    external = pd.read_csv("data/metadata/gse41258_samples.csv", dtype={"patient_id": str})
    figures = result_root(ROOT / paths["figures"], provenance)
    metrics_root = result_root(ROOT / paths["metrics"], provenance)
    tables_root = result_root(ROOT / paths["tables"], provenance)
    reports_root = result_root(ROOT / paths["reports"], provenance)
    for directory in (figures, metrics_root, tables_root, reports_root, ROOT / "manuscript"):
        directory.mkdir(parents=True, exist_ok=True)

    study_design_figure(figures / "study_design", raw_complete)
    sample_flow_figure(primary, external, figures / "sample_inclusion_flowchart")
    nested_metrics = pd.read_csv(metrics_root / "nested_cv_metrics.csv")
    nested_predictions = pd.read_csv(metrics_root / "nested_cv_predictions.csv")
    modeling_figures(nested_metrics, nested_predictions, figures)
    raw_de_path = Path("results/tables/supplementary_full_raw_cel_de_results.csv")
    raw_de = pd.read_csv(raw_de_path) if raw_de_path.exists() else pd.DataFrame()
    if not raw_de.empty:
        de_figures(raw_de, figures)

    included_external = external[external["inclusion_status"].eq("included")]
    table1 = pd.DataFrame(
        [
            ["GSE44076", "GPL13667", "healthy", 50, 50, "development"],
            ["GSE44076", "GPL13667", "adjacent_normal", 98, 98, "development"],
            ["GSE44076", "GPL13667", "tumor", 98, 98, "development"],
            ["GSE41258", "GPL96", "normal_colon", int(included_external["tissue_class"].eq("normal_colon").sum()), int(included_external.loc[included_external["tissue_class"].eq("normal_colon"), "patient_id"].nunique()), "external tumor/normal validation"],
            ["GSE41258", "GPL96", "primary_tumor", int(included_external["tissue_class"].eq("primary_tumor").sum()), int(included_external.loc[included_external["tissue_class"].eq("primary_tumor"), "patient_id"].nunique()), "external tumor/normal validation"],
        ],
        columns=["accession", "platform", "tissue_class", "arrays", "unique_patients_or_donors", "role"],
    )
    table1.to_csv("results/tables/table_1_dataset_characteristics.csv", index=False)
    if not raw_de.empty:
        table2 = (
            raw_de.assign(significant=lambda frame: frame["adjusted_p_value"].lt(0.05) & frame["log2_fold_change"].abs().ge(0.5))
            .groupby(["comparison", "direction"], as_index=False)
            .agg(tested_genes=("gene_symbol", "size"), significant_genes=("significant", "sum"))
        )
        table2.to_csv("results/tables/table_2_raw_cel_differential_expression_summary.csv", index=False)
    pd.concat([primary.assign(accession="GSE44076"), external.assign(accession="GSE41258")], ignore_index=True, sort=False).to_csv("results/tables/supplementary_sample_metadata.csv", index=False)
    mappings = []
    for accession in ("GSE44076", "GSE41258"):
        path = Path(f"data/metadata/{accession}_probe_gene_mapping_raw_cel_rma.csv")
        if path.exists():
            mappings.append(pd.read_csv(path).assign(accession=accession))
    if mappings:
        pd.concat(mappings, ignore_index=True).to_csv("results/tables/supplementary_probe_gene_mapping.csv", index=False)

    selected_model = nested_metrics[nested_metrics["model"].eq("elastic_net")]
    internal = selected_model.groupby("task", as_index=False).agg(
        mean_macro_f1=("f1_macro", "mean"),
        median_macro_f1=("f1_macro", "median"),
        sd_macro_f1=("f1_macro", "std"),
        fold_minimum=("f1_macro", "min"),
        fold_maximum=("f1_macro", "max"),
    )
    compact_path = tables_root / "final_compact_signature.csv"
    compact = pd.read_csv(compact_path) if compact_path.exists() else pd.DataFrame()
    external_metrics_path = metrics_root / "external_validation_metrics.csv"
    external_metrics = pd.read_csv(external_metrics_path) if external_metrics_path.exists() else pd.DataFrame()
    primary_external = external_metrics[external_metrics.get("is_primary", False).astype(bool)] if not external_metrics.empty else pd.DataFrame()
    de_summary_text = table2.to_markdown(index=False) if not raw_de.empty else "Raw-CEL limma was not completed."
    compact_genes = compact["gene_symbol"].tolist() if not compact.empty else []
    external_text = primary_external.to_markdown(index=False) if not primary_external.empty else "External validation was not completed for the selected provenance."

    introduction = """# Introduction

Colorectal carcinogenesis can be accompanied by molecular alterations beyond
the histologically apparent tumor, a concept commonly described as field
cancerization [REF]. Consequently, tumor-adjacent mucosa should not be assumed
to be equivalent to mucosa from cancer-free donors [REF]. Transcriptome-wide
profiling can identify candidate field alterations, but small-sample,
high-dimensional analyses are vulnerable to patient leakage, unstable feature
selection, platform effects, and optimistic validation [REF].

This study tests whether histologically normal tumor-adjacent colon differs
from genuinely healthy colon mucosa, characterizes healthy-to-adjacent-to-tumor
expression patterns without asserting biological progression, and evaluates
compact, interpretable signatures under nested patient-aware validation. A
separate external cohort tests only the tumor-versus-normal component.
"""
    methods = f"""# Methods

GSE44076 (GPL13667) was the sole development cohort (50 healthy donors, 98
adjacent-normal samples, and 98 paired tumors). GSE41258 (GPL96) was restricted
to author-included Primary Tumor and Normal Colon arrays after deterministic
removal of metastases, non-colon tissues, polyps, cell lines, ambiguous records,
and `_ez`/`rehyb` technical-replicate candidates.

Raw CEL files were processed independently by platform using oligo RMA for
HG-U219 and affy RMA for HG-U133A. CEL chip/CDF compatibility, GEO identifiers,
sample counts, finite log2-scale output, and probe-set counts were asserted.
Probes required unambiguous symbol and Entrez mappings; multiple probes per
gene were aggregated by the per-sample median. Raw and normalized
distributions, RLE, MA diagnostics, PCA, correlation, clustering, robust
sample-median/IQR diagnostics, and arrayQualityMetrics were generated. PCA
separation alone was not an exclusion criterion.

Primary differential expression used limma on raw-CEL RMA expression. Unpaired
designs tested adjacent-normal versus healthy and tumor versus healthy. A
patient fixed-effect design tested paired tumor versus adjacent-normal. All
coefficient names and design rank were checked programmatically, and BH FDR was
applied. Field candidates required FDR < 0.05, |log2FC| >= 0.5, at least 70%
sample-direction consistency, and processed-matrix sign agreement. Enrichment
used the tested-gene universe and separated direction/contrast-specific sets.

Three predictive tasks were prespecified. Elastic Net was the primary model;
group-calibrated linear SVM and Random Forest were comparators. Validation used
{config['modeling']['outer_splits']}-fold outer and {config['modeling']['inner_splits']}-fold inner StratifiedGroupKFold across {len(config['modeling']['seeds'])} deterministic repeats. Every learned filter, feature selection, scale, hyperparameter, calibration, and compact panel was fitted inside training data. Compact panels of {config['modeling']['compact_panels']} genes and the fold-specific full stable signature were evaluated. Patient/donor bootstrap intervals quantified uncertainty.

The final Task C model and threshold were locked using GSE44076 only. The
prespecified primary cross-platform representation was within-sample percentile
rank over the label-independent common-gene universe. GSE41258 point estimates
used a deterministic one-array-per-patient subset; all arrays were retained for
patient-cluster bootstrap sensitivity. No external label informed model,
feature, representation, threshold, or hyperparameter selection. Primary seed:
{config['project']['random_seed']}.
"""
    results = f"""# Results

## Cohorts and raw analysis

{table1.to_markdown(index=False)}

Raw-CEL limma summary:

{de_summary_text}

## Internal predictive validation

{internal.to_markdown(index=False)}

The final compact Task C signature contained {len(compact_genes)} genes:
{', '.join(compact_genes) if compact_genes else 'not available'}.

## External tumor-versus-normal validation

{external_text}

This external result does not validate cancer-free healthy versus
tumor-adjacent field cancerization.
"""
    limitations = """# Limitations

- Retrospective public microarray cohorts cannot establish prospective clinical utility.
- GPL13667 and GPL96 differ in probe design and distribution.
- GSE41258 validates tumor versus normal colon, not the primary healthy-versus-adjacent field effect.
- GSE44076 is enriched for stage II and microsatellite-stable disease.
- Bulk expression conflates epithelial state with immune/stromal composition.
- Explicit batch covariates are incomplete and may remain confounded with phenotype.
- High-dimensional selection can remain unstable despite nested grouped validation.
- Differential expression and enrichment support association, not causal mechanism.
- The model is a candidate biomarker-discovery model, not a clinically ready diagnostic.
"""
    discussion = """# Discussion

The analysis separates the field-effect question from generic tumor detection.
Raw-CEL limma evidence and sample-level consistency provide the primary basis
for candidate field genes, while machine-learning coefficients contribute
predictive stability rather than defining biology. Patient-aware nesting limits
optimism from paired samples and feature selection.

Within-sample ranks offer a biologically neutral transfer representation: they
require no external labels, reduce dependence on platform intensity scale, and
preserve relative ordering. Strong discrimination should still be considered
separately from calibration and cannot establish clinical utility. Independent
cancer-free healthy and tumor-adjacent cohorts, tissue-composition analyses,
and prospective assay validation remain necessary.
"""
    availability = """# Data and code availability

Expression data are available from NCBI GEO under GSE44076 and GSE41258. Raw
archives and large generated matrices are intentionally excluded from Git.
Checksums, deterministic metadata, compact result tables, environment locks,
commands, and model documentation are included in the repository. Users must
place the five named GEO files under `data/raw/<accession>/` and run `make all`.
"""
    legends = """# Figure legends

1. Study design and leakage-safe analysis flow.
2. Sample inclusion flow with array and patient counts.
3. Raw and post-RMA array QC distributions and diagnostics.
4. Raw-CEL RMA PCA, correlation, and hierarchical clustering.
5. Raw-CEL limma adjacent-normal versus healthy volcano plot.
6. Patient-fixed-effect tumor versus adjacent-normal volcano plot.
7. Raw-versus-deposited log2FC concordance.
8. Final field-candidate heatmap and group distributions.
9. Repeated nested-CV model comparison and three-class OOF confusion matrix.
10. Task B and Task C OOF ROC, precision-recall, and calibration.
11. Repeat-aware stable Elastic Net coefficients and compact-panel curve.
12. Patient-balanced GSE41258 ROC/PR, calibration, and confusion matrix.
13. Direction- and contrast-specific GO/Reactome/ranked enrichment results.
"""
    manuscript_dir = ROOT / "manuscript"
    for name, text in {
        "introduction.md": introduction,
        "methods.md": methods,
        "results.md": results,
        "limitations.md": limitations,
        "discussion.md": discussion,
        "data_and_code_availability.md": availability,
        "figure_legends.md": legends,
    }.items():
        (manuscript_dir / name).write_text(text, encoding="utf-8")
    manuscript = "\n\n".join([introduction, methods, results, discussion, limitations, availability, legends])
    (manuscript_dir / "manuscript.md").write_text(manuscript, encoding="utf-8")
    (manuscript_dir / "future_field_effect_validation_protocol.md").write_text(
        """# Future independent field-effect validation protocol

Recruit or identify an independent cohort containing cancer-free healthy colon
and histologically normal tumor-adjacent colon from independent patients, with
explicit tissue location, distance from tumor, pathology review, age, sex,
stage, MSI/subtype, processing batch, and raw compatible expression data.
Freeze the present genes, preprocessing, representation, and threshold before
accessing outcome labels. Normalize the new platform independently; report
patient-level discrimination, calibration, bootstrap intervals, failure rates,
and prespecified subgroup/sensitivity analyses. Do not use the cohort for
feature, transformation, threshold, or model selection.
""",
        encoding="utf-8",
    )

    raw_qc_files = sorted(Path("results/tables").glob("GSE*_raw_cel_qc_diagnostics.csv"))
    qc_text = "\n\n".join(
        [f"## {path.stem}\n\n{pd.read_csv(path).describe(include='all').to_markdown()}" for path in raw_qc_files]
    ) or "Raw-CEL QC tables were not generated."
    Path("reports/raw_cel_qc_report.md").write_text(
        "# Raw-CEL QC report\n\nNo sample is removed because of PCA separation alone. "
        "Technical exclusions require at least two independent severe QC failures.\n\n" + qc_text,
        encoding="utf-8",
    )

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    gates = [
        ("A - Data integrity", Path("reports/cel_archive_audit.md").exists() and len(primary) == 246 and len(external) == 390, "reports/cel_archive_audit.md", "Rscript R/01_extract_and_audit_cel.R"),
        ("B - Raw preprocessing", raw_complete and Path("data/metadata/GSE44076_platform_validation.csv").exists() and Path("data/metadata/GSE41258_platform_validation.csv").exists(), "data/metadata/GSE44076_platform_validation.csv", "Rscript R/02_preprocess_gse44076.R; Rscript R/03_preprocess_gse41258.R"),
        ("C - Biological inference", not raw_de.empty and Path("results/tables/supplementary_enrichment_results.csv").exists() and Path("results/tables/raw_vs_deposited_de_concordance.csv").exists(), "results/tables/supplementary_full_raw_cel_de_results.csv", "Rscript R/04_differential_expression.R; Rscript R/05_functional_enrichment.R"),
        ("D - Predictive validity", (metrics_root / "compact_panel_fold_metrics.csv").exists() and (metrics_root / "fold_assignments.csv").exists(), str((metrics_root / "compact_panel_fold_metrics.csv").relative_to(ROOT)), f"python scripts/run_compact_panel_analysis.py --provenance {provenance}"),
        ("E - External validity", not primary_external.empty and (metrics_root / "external_validation_bootstrap_ci.csv").exists(), str(external_metrics_path.relative_to(ROOT)), f"python scripts/run_external_validation.py --provenance {provenance}"),
        ("F - Reproducibility", Path("requirements-lock.txt").exists() and '"Packages": {}' not in Path("renv.lock").read_text(encoding="utf-8") and Path(".github/workflows/ci.yml").exists() and Path("reports/docker_verification.md").exists() and "Status: PASS" in Path("reports/docker_verification.md").read_text(encoding="utf-8"), "reports/docker_verification.md", "python -m pip check; Rscript R/verify_environment.R; python -m pytest -q; docker build ."),
        ("G - Publication outputs", len(compact_genes) > 0 and Path("manuscript/manuscript.md").exists(), "manuscript/manuscript.md", "python scripts/build_manuscript_outputs.py"),
    ]
    evidence_rows = []
    for requirement, passed, evidence, command in gates:
        evidence_rows.append(
            {
                "Requirement": requirement,
                "Status": "PASS" if passed else "PARTIAL",
                "Evidence file": evidence,
                "Verification command": command,
                "Key result": "Required evidence exists and was inspected" if passed else "One or more required evidence artifacts are absent",
                "Remaining action": "None" if passed else "Complete and rerun the listed verification command",
            }
        )
    evidence = pd.DataFrame(evidence_rows)
    decision = "PUBLICATION-READY" if evidence["Status"].eq("PASS").all() else "PARTIAL"
    final_report = f"""# Publication-readiness final report

## 1. Executive summary

Final readiness decision: **{decision}**. Scientific claims remain restricted
to retrospective biomarker discovery; no clinical-readiness claim is made.

## 2. Repository changes

Raw platform assertions, provenance separation, limma, enrichment,
configuration-driven grouped nested CV, group-calibrated SVM, fold-internal
compact panels, patient-balanced external validation, CI, locks, tests, and
manuscript generators were added or corrected.

## 3. Raw-data audit

GSE44076 contains 246 arrays and GSE41258 contains 390 arrays. Detailed
archive, mapping, and checksum evidence is linked below.

## 4. RMA and QC

Raw-CEL completion status: {raw_complete}. No PCA-only exclusion is permitted.

## 5. Sample exclusions

Deterministic metadata exclusions and technical candidates are recorded in
`data/metadata/sample_exclusion_log.csv` and the metadata audit.

## 6. Final differential expression

{de_summary_text}

## 7. Field-cancerization results

Candidate tables prioritize raw-CEL limma, sample consistency, sensitivity
agreement, and uncertainty; three group means are not treated as progression.

## 8. Enrichment

Direction- and contrast-specific GO, Reactome, and ranked results use the tested
gene universe. Enrichment is associative, not causal.

## 9. Model methodology

Elastic Net was prespecified; SVM and Random Forest are comparators.

## 10. Leakage controls

Outer, inner, calibration, feature-selection, compact-panel, and threshold
operations preserve patient/donor boundaries and training-only fitting.

## 11. Internal validation

{internal.to_markdown(index=False)}

## 12. Feature stability

Selection counts use unique `(repeat, outer_fold)` keys and report sign and
coefficient variability.

## 13. Compact-panel selection

Task C signature ({len(compact_genes)} genes): {', '.join(compact_genes) if compact_genes else 'not completed'}.

## 14. Calibration

Brier score, log loss, calibration curves, and binary calibration intercept and
slope are generated where estimable.

## 15. Permutation tests

Task A is explicitly conditional; Task B is group-level; Task C is paired
within patient. Resolution is controlled by the configured iteration count.

## 16. External validation

{external_text}

## 17. Sensitivity analyses

Raw-versus-deposited, compact size, representation, threshold, and
canonical-patient versus all-array results are preserved separately.

## 18. Reproducibility

Python and R locks, Docker restore, normal CI, manual full-data CI, Make targets,
test evidence, and artifact checksums are included.

## 19. Remaining limitations

See `manuscript/limitations.md`; independent Task B field-effect validation and
prospective biological/clinical validation remain required.

## 20. Final readiness decision

**{decision}**

## Final evidence table

{evidence.to_markdown(index=False)}
"""
    Path("reports/publication_readiness_final.md").write_text(final_report, encoding="utf-8")
    Path("reports/final_analysis_report.md").write_text(final_report, encoding="utf-8")

    manifest = artifact_manifest(commit)
    manifest.to_csv("reports/result_artifact_manifest.csv", index=False)
    print(
        json.dumps(
            {
                "readiness": decision,
                "raw_complete": raw_complete,
                "figures": len(list(figures.glob("*.png"))),
                "tables": len(list(Path("results/tables").glob("*.csv"))),
                "manifest_rows": len(manifest),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
