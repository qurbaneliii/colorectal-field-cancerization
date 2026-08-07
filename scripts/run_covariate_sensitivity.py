from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.visualization.publication_figures import save_figure


def read_de(name: str) -> pd.DataFrame:
    path = ROOT / "results/differential_expression" / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path).set_index("gene_symbol", drop=False)


def significant(frame: pd.DataFrame, fdr: float, effect: float) -> pd.Series:
    return frame["adjusted_p_value"].lt(fdr) & frame["log2_fold_change"].abs().ge(effect)


def concordance_row(
    reference: pd.DataFrame,
    adjusted: pd.DataFrame,
    reference_model: str,
    adjusted_model: str,
    fdr: float,
    effect: float,
) -> tuple[dict[str, object], pd.DataFrame]:
    merged = reference[["gene_symbol", "log2_fold_change", "adjusted_p_value"]].reset_index(
        drop=True
    ).merge(
        adjusted[["gene_symbol", "log2_fold_change", "adjusted_p_value"]].reset_index(
            drop=True
        ),
        on="gene_symbol",
        suffixes=("_reference", "_adjusted"),
        validate="one_to_one",
    )
    merged["reference_significant"] = (
        merged["adjusted_p_value_reference"].lt(fdr)
        & merged["log2_fold_change_reference"].abs().ge(effect)
    )
    merged["adjusted_significant"] = (
        merged["adjusted_p_value_adjusted"].lt(fdr)
        & merged["log2_fold_change_adjusted"].abs().ge(effect)
    )
    merged["sign_agreement"] = np.sign(merged["log2_fold_change_reference"]) == np.sign(
        merged["log2_fold_change_adjusted"]
    )
    merged["lost_significance"] = merged["reference_significant"] & ~merged[
        "adjusted_significant"
    ]
    merged["gained_significance"] = ~merged["reference_significant"] & merged[
        "adjusted_significant"
    ]
    merged["changed_direction"] = ~merged["sign_agreement"]
    reference_hits = merged["reference_significant"].sum()
    retained = (merged["reference_significant"] & merged["adjusted_significant"]).sum()
    reference_effect = merged.loc[merged["reference_significant"], "log2_fold_change_reference"].abs()
    adjusted_effect = merged.loc[merged["reference_significant"], "log2_fold_change_adjusted"].abs()
    return (
        {
            "reference_model": reference_model,
            "adjusted_model": adjusted_model,
            "shared_tested_genes": len(merged),
            "pearson_log2fc_correlation": stats.pearsonr(
                merged["log2_fold_change_reference"], merged["log2_fold_change_adjusted"]
            ).statistic,
            "spearman_log2fc_correlation": stats.spearmanr(
                merged["log2_fold_change_reference"], merged["log2_fold_change_adjusted"]
            ).statistic,
            "sign_agreement_fraction": merged["sign_agreement"].mean(),
            "reference_significant_genes": int(reference_hits),
            "adjusted_significant_genes": int(merged["adjusted_significant"].sum()),
            "significant_gene_retention_count": int(retained),
            "significant_gene_retention_fraction": float(retained / reference_hits)
            if reference_hits
            else np.nan,
            "median_effect_size_retention": float(np.median(adjusted_effect / reference_effect))
            if len(reference_effect)
            else np.nan,
            "genes_losing_significance": int(merged["lost_significance"].sum()),
            "genes_gaining_significance": int(merged["gained_significance"].sum()),
            "genes_changing_direction": int(merged["changed_direction"].sum()),
            "pathway_concordance": "evaluated after tier-specific enrichment",
        },
        merged,
    )


def subgroup_support(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: pd.Index,
    expected_direction: pd.Series,
    minimum_n: int,
) -> pd.DataFrame:
    selected = metadata[metadata["tissue_class"].isin(["healthy", "adjacent_normal"])].copy()
    selected["age"] = pd.to_numeric(selected["age"], errors="coerce")
    selected["age_stratum"] = np.where(
        selected["age"] <= selected["age"].median(), "age_lower", "age_upper"
    )
    strata = {
        "sex": selected["sex"].astype(str),
        "location": selected["location"].astype(str),
        "age_stratum": selected["age_stratum"].astype(str),
    }
    effects: dict[str, pd.Series] = {}
    for variable, values in strata.items():
        for level in sorted(values.dropna().unique()):
            index = selected.index[values.eq(level)]
            group = selected.loc[index]
            counts = group["tissue_class"].value_counts()
            if (
                counts.get("healthy", 0) < minimum_n
                or counts.get("adjacent_normal", 0) < minimum_n
            ):
                continue
            healthy_ids = group.loc[
                group["tissue_class"].eq("healthy"), "geo_accession"
            ]
            adjacent_ids = group.loc[
                group["tissue_class"].eq("adjacent_normal"), "geo_accession"
            ]
            effects[f"{variable}={level}"] = (
                expression.loc[genes, adjacent_ids].mean(axis=1)
                - expression.loc[genes, healthy_ids].mean(axis=1)
            )
    matrix = pd.DataFrame(effects, index=genes)
    agreement = np.sign(matrix).eq(expected_direction.loc[genes], axis=0)
    return pd.DataFrame(
        {
            "subgroups_tested": matrix.notna().sum(axis=1),
            "subgroup_labels": ";".join(matrix.columns),
            "subgroup_sign_agreement_fraction": agreement.mean(axis=1),
            "subgroup_effect_min": matrix.min(axis=1),
            "subgroup_effect_max": matrix.max(axis=1),
            "major_demographic_subgroup_dependence": ~agreement.all(axis=1),
        },
        index=genes,
    )


def main() -> None:
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    biology = config["biology"]
    fdr = float(biology["fdr_threshold"])
    effect = float(biology["log2fc_threshold"])
    high_fdr = float(biology["high_confidence_fdr_threshold"])
    high_effect = float(biology["high_confidence_log2fc_threshold"])
    consistency_threshold = float(biology["direction_consistency_threshold"])
    high_consistency = float(biology["high_confidence_direction_consistency"])
    minimum_n = int(biology["subgroup_minimum_n"])

    u0 = read_de("raw_cel_adjacent_vs_healthy_unadjusted.csv")
    u1 = read_de("raw_cel_adjacent_vs_healthy_age_sex_adjusted.csv")
    u2 = read_de("raw_cel_adjacent_vs_healthy_age_sex_location_adjusted.csv")
    qc = read_de("raw_cel_adjacent_vs_healthy_age_sex_adjusted_qc_excluded.csv")
    deposited = pd.read_csv(
        ROOT / "results/differential_expression/processed_matrix_adjacent_normal_vs_healthy.csv"
    ).set_index("gene_symbol", drop=False)
    expression = pd.read_parquet(
        ROOT / "data/processed/GSE44076_gene_expression_raw_cel_rma.parquet"
    ).set_index("gene_symbol")
    metadata = pd.read_csv(
        ROOT / "data/metadata/gse44076_samples.csv", dtype={"patient_id": str}
    )
    trajectories = pd.read_csv(ROOT / "results/tables/final_expression_trajectory_genes.csv").set_index(
        "gene_symbol"
    )

    concordance_rows = []
    sensitive_frames = []
    for adjusted, adjusted_name in [(u1, "U1_age_sex_adjusted"), (u2, "U2_age_sex_location_adjusted")]:
        row, merged = concordance_row(u0, adjusted, "U0_unadjusted", adjusted_name, fdr, effect)
        concordance_rows.append(row)
        changed = merged[merged["lost_significance"] | merged["changed_direction"]].copy()
        changed["adjusted_model"] = adjusted_name
        sensitive_frames.append(changed)
    concordance = pd.DataFrame(concordance_rows)
    concordance.to_csv(ROOT / "results/tables/covariate_adjusted_de_concordance.csv", index=False)
    pd.concat(sensitive_frames, ignore_index=True).to_csv(
        ROOT / "results/tables/covariate_sensitive_field_genes.csv", index=False
    )

    common = u0.index.intersection(u1.index).intersection(u2.index).intersection(qc.index)
    common = common.intersection(deposited.index).intersection(expression.index)
    evidence = pd.DataFrame(index=common)
    evidence.index.name = "gene_symbol"
    for prefix, frame in [("u0", u0), ("u1", u1), ("u2", u2), ("qc", qc)]:
        evidence[f"{prefix}_log2fc"] = frame.loc[common, "log2_fold_change"]
        evidence[f"{prefix}_fdr"] = frame.loc[common, "adjusted_p_value"]
    evidence["deposited_log2fc"] = deposited.loc[common, "log2_fold_change"]
    evidence["deposited_fdr"] = deposited.loc[common, "adjusted_p_value"]
    expected_sign = np.sign(evidence["u1_log2fc"])

    field_meta = metadata[metadata["tissue_class"].isin(["healthy", "adjacent_normal"])]
    healthy = field_meta.loc[field_meta["tissue_class"].eq("healthy"), "geo_accession"]
    adjacent = field_meta.loc[field_meta["tissue_class"].eq("adjacent_normal"), "geo_accession"]
    healthy_median = expression.loc[common, healthy].median(axis=1)
    adjacent_values = expression.loc[common, adjacent]
    up_consistency = adjacent_values.gt(healthy_median, axis=0).mean(axis=1)
    down_consistency = adjacent_values.lt(healthy_median, axis=0).mean(axis=1)
    evidence["sample_direction_consistency"] = np.where(expected_sign >= 0, up_consistency, down_consistency)

    subgroup = subgroup_support(expression, metadata, common, expected_sign, minimum_n)
    evidence = evidence.join(subgroup)
    evidence = evidence.join(
        trajectories[[
            "healthy", "adjacent_normal", "tumor", "healthy_to_adjacent",
            "adjacent_to_tumor", "trajectory_category"
        ]],
        how="left",
    )

    evidence["u0_screen"] = significant(u0.loc[common], fdr, effect)
    evidence["u1_adjusted_support"] = significant(u1.loc[common], high_fdr, high_effect)
    evidence["u2_location_support"] = significant(u2.loc[common], fdr, effect)
    evidence["qc_robust"] = significant(qc.loc[common], fdr, effect)
    evidence["unadjusted_adjusted_sign_agreement"] = np.sign(evidence["u0_log2fc"]) == expected_sign
    evidence["location_adjusted_sign_agreement"] = np.sign(evidence["u2_log2fc"]) == expected_sign
    evidence["qc_sign_agreement"] = np.sign(evidence["qc_log2fc"]) == expected_sign
    evidence["deposited_support"] = (
        (np.sign(evidence["deposited_log2fc"]) == expected_sign)
        & evidence["deposited_fdr"].lt(fdr)
        & evidence["deposited_log2fc"].abs().ge(effect)
    )
    evidence["direction_consistent"] = evidence["sample_direction_consistency"].ge(
        consistency_threshold
    )
    evidence["high_direction_consistency"] = evidence["sample_direction_consistency"].ge(
        high_consistency
    )
    evidence["subgroup_robust"] = (
        evidence["subgroups_tested"].ge(4)
        & evidence["subgroup_sign_agreement_fraction"].eq(1.0)
    )

    components = pd.DataFrame(
        {
            "adjusted_statistical_evidence": np.minimum(
                -np.log10(evidence["u1_fdr"].clip(lower=1e-300)) / 20.0, 1.0
            ),
            "adjusted_effect_size": np.minimum(evidence["u1_log2fc"].abs() / 2.0, 1.0),
            "sample_direction_consistency": evidence["sample_direction_consistency"],
            "preprocessing_concordance": evidence["deposited_support"].astype(float),
            "qc_robustness": evidence["qc_robust"].astype(float),
            "demographic_robustness": evidence["subgroup_robust"].astype(float),
            "trajectory_field_support": evidence["trajectory_category"].isin(
                ["field_shift_then_plateau", "monotonic_up", "monotonic_down"]
            ).astype(float),
        }
    )
    evidence["evidence_score"] = components.mean(axis=1)
    tier1 = (
        evidence["u1_adjusted_support"]
        & evidence["unadjusted_adjusted_sign_agreement"]
        & evidence["u2_location_support"]
        & evidence["location_adjusted_sign_agreement"]
        & evidence["high_direction_consistency"]
        & evidence["deposited_support"]
        & evidence["qc_robust"]
        & evidence["qc_sign_agreement"]
        & evidence["subgroup_robust"]
        & evidence["trajectory_category"].isin(
            ["field_shift_then_plateau", "monotonic_up", "monotonic_down"]
        )
    )
    tier2 = evidence["u0_screen"] & ~tier1 & (
        (
            evidence["u1_fdr"].lt(fdr)
            & evidence["unadjusted_adjusted_sign_agreement"]
        )
        | (evidence["direction_consistent"] & evidence["deposited_support"])
    )
    tier3 = evidence["u0_screen"] & ~tier1 & ~tier2
    evidence["evidence_tier"] = np.select(
        [tier1, tier2, tier3],
        ["Tier 1 - high-confidence", "Tier 2 - provisional", "Tier 3 - exploratory"],
        default="not_in_initial_screen",
    )
    evidence = evidence.reset_index().sort_values(
        ["evidence_tier", "evidence_score", "u1_fdr"], ascending=[True, False, True]
    )
    evidence.to_csv(ROOT / "results/tables/field_gene_evidence_tiers.csv", index=False)
    high = evidence[evidence["evidence_tier"].eq("Tier 1 - high-confidence")].copy()
    provisional = evidence[evidence["evidence_tier"].eq("Tier 2 - provisional")].copy()
    exploratory = evidence[evidence["evidence_tier"].eq("Tier 3 - exploratory")].copy()
    high.to_csv(ROOT / "results/tables/final_high_confidence_field_signature.csv", index=False)
    provisional.to_csv(ROOT / "results/tables/provisional_field_associated_genes.csv", index=False)
    exploratory.to_csv(ROOT / "results/tables/exploratory_field_gene_universe.csv", index=False)

    plot = u0[["log2_fold_change"]].rename(columns={"log2_fold_change": "unadjusted"}).join(
        u1[["log2_fold_change"]].rename(columns={"log2_fold_change": "age_sex_adjusted"}),
        how="inner",
    )
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.scatter(plot["unadjusted"], plot["age_sex_adjusted"], s=7, alpha=0.28, color="#457B9D")
    limits = [plot.to_numpy().min(), plot.to_numpy().max()]
    ax.plot(limits, limits, linestyle="--", color="black", linewidth=0.8)
    ax.set(
        xlabel="U0 unadjusted log2 fold change",
        ylabel="U1 age/sex-adjusted log2 fold change",
        title=(
            "Adjacent-normal versus healthy covariate sensitivity\n"
            f"Pearson r={concordance.iloc[0]['pearson_log2fc_correlation']:.3f}"
        ),
    )
    save_figure(fig, ROOT / "results/figures/unadjusted_vs_adjusted_field_logfc")

    qc_row, qc_merged = concordance_row(
        u1, qc, "U1_age_sex_adjusted", "U1_age_sex_adjusted_qc_sensitivity", fdr, effect
    )
    qc_diagnostics = pd.read_csv(ROOT / "results/tables/GSE44076_raw_cel_qc_diagnostics.csv")
    borderline = qc_diagnostics.loc[
        qc_diagnostics["independent_failure_metrics"].eq(1), "geo_accession"
    ].tolist()
    qc_table = pd.DataFrame(
        [
            {"analysis_component": "borderline_samples", "primary_value": 0, "sensitivity_value": len(borderline), "comparison": ";".join(borderline)},
            {"analysis_component": "significant_DE_genes", "primary_value": int(significant(u1, fdr, effect).sum()), "sensitivity_value": int(significant(qc, fdr, effect).sum()), "comparison": "FDR<0.05 and |log2FC|>=0.5"},
            {"analysis_component": "top_100_gene_overlap", "primary_value": 100, "sensitivity_value": len(set(u1.nsmallest(100, "adjusted_p_value").index) & set(qc.nsmallest(100, "adjusted_p_value").index)), "comparison": "ranked by adjusted p-value"},
            {"analysis_component": "log2FC_pearson", "primary_value": 1.0, "sensitivity_value": qc_row["pearson_log2fc_correlation"], "comparison": "all shared tested genes"},
        ]
    )
    qc_table.to_csv(ROOT / "results/tables/qc_exclusion_sensitivity.csv", index=False)

    report = f"""# Covariate-adjustment sensitivity

The primary biological comparison is adjacent-normal versus genuinely healthy
colon mucosa. U0 is unadjusted, U1 adjusts for continuous age and sex, and U2
additionally adjusts for left/right tumor location. All three designs were
full-rank; stage was not used because it is structurally undefined for healthy
donors.

{concordance.to_markdown(index=False)}

At FDR < {fdr} and |log2FC| >= {effect}, U0 detected
{int(significant(u0, fdr, effect).sum()):,} genes, U1 detected
{int(significant(u1, fdr, effect).sum()):,}, and U2 detected
{int(significant(u2, fdr, effect).sum()):,}. Covariate sensitivity is reported
as concordance and retention, not as proof that either model is causally
correct.

The evidence-tier rules were declared in `config/analysis.yaml`. Tier 1 requires
strong U1 evidence (FDR < {high_fdr}, |adjusted log2FC| >= {high_effect}), U0/U1
direction agreement, U2 support, at least {high_consistency:.0%} sample-direction
consistency, deposited-matrix support, one-sample QC sensitivity robustness,
consistent direction in all adequately represented age/sex/location strata,
and a field-compatible trajectory. Machine-learning coefficients do not define
this biological tier.

- Tier 1 high-confidence genes: {len(high)}
- Tier 2 provisional genes: {len(provisional)}
- Tier 3 exploratory genes: {len(exploratory)}
"""
    (ROOT / "reports/covariate_adjustment_sensitivity.md").write_text(report, encoding="utf-8")
    qc_report = f"""# QC exclusion sensitivity

No GSE44076 sample met the prespecified exclusion rule of at least two
independent severe technical failures. {', '.join(borderline) if borderline else 'No sample'}
had exactly one severe custom-QC flag and was retained in the primary analysis.
A secondary U1 age/sex-adjusted limma model omitted that sample.

{qc_table.to_markdown(index=False)}

Pathway and Task B model comparisons are populated after the enrichment and
model workflows run. PCA separation was not used to exclude any sample.
"""
    (ROOT / "reports/qc_exclusion_sensitivity.md").write_text(qc_report, encoding="utf-8")
    print(
        json.dumps(
            {
                "u0_significant": int(significant(u0, fdr, effect).sum()),
                "u1_significant": int(significant(u1, fdr, effect).sum()),
                "u2_significant": int(significant(u2, fdr, effect).sum()),
                "tier_1": len(high),
                "tier_2": len(provisional),
                "tier_3": len(exploratory),
                "borderline_qc_samples": borderline,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
