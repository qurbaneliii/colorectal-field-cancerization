from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    os.chdir(ROOT)
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text(encoding="utf-8"))
    curated_path = ROOT / config["biology"]["stress_gene_set"]
    curated = pd.read_csv(curated_path, sep="\t")
    tiers = pd.read_csv(ROOT / "results/tables/field_gene_evidence_tiers.csv")
    high = tiers[tiers["evidence_tier"].eq("Tier 1 - high-confidence")].copy()
    top = high.sort_values(["evidence_score", "u1_fdr"], ascending=[False, True]).head(100)
    audit = curated.merge(
        tiers[
            [
                "gene_symbol",
                "evidence_tier",
                "evidence_score",
                "u1_log2fc",
                "u1_fdr",
                "composition_adjusted_log2fc",
                "composition_adjusted_fdr",
                "composition_robust",
                "trajectory_category",
            ]
        ],
        on="gene_symbol",
        how="left",
        validate="one_to_one",
    )
    audit["present_in_tested_matrix"] = audit["u1_fdr"].notna()
    audit["in_high_confidence_field_signature"] = audit["evidence_tier"].eq(
        "Tier 1 - high-confidence"
    )
    audit["in_top_100_high_confidence"] = audit["gene_symbol"].isin(top["gene_symbol"])
    audit["removal_policy"] = "flag_only_not_automatically_removed"
    audit.to_csv(ROOT / "results/tables/immediate_early_stress_gene_audit.csv", index=False)

    stress = set(curated["gene_symbol"])
    without_stress = high[~high["gene_symbol"].isin(stress)].copy()
    without_stress.to_csv(
        ROOT / "results/tables/high_confidence_field_signature_without_stress_genes.csv",
        index=False,
    )
    high_stress = high[high["gene_symbol"].isin(stress)]
    report = f"""# Preanalytical stress and immediate-early sensitivity

The prespecified flag set contains {len(curated)} immediate-early, AP-1,
feedback, and acute-response genes documented from the curated immediate-early
response analysis of Tullai et al. (Journal of Biological Chemistry 2007;
282:23981-23995; DOI: 10.1074/jbc.M703225200). The set is versioned at
`{curated_path.relative_to(ROOT).as_posix()}`. These genes were flagged rather
than automatically removed because acute signaling can reflect tissue handling,
true mucosal biology, or both.

- Curated genes present in the tested matrix: {int(audit['present_in_tested_matrix'].sum())}
- Flagged genes in the {len(high)}-gene high-confidence field set: {len(high_stress)}
- Prevalence in the high-confidence field set: {len(high_stress) / len(high):.1%}
- Flagged genes among the top 100 high-confidence candidates: {int(audit['in_top_100_high_confidence'].sum())}
- High-confidence set after flag-only removal sensitivity: {len(without_stress)} genes

Flagged high-confidence genes: {', '.join(high_stress['gene_symbol']) if len(high_stress) else 'none'}.

Tier-specific enrichment is repeated both with and without these flagged genes.
Conclusions are treated as stress-dependent only when the nonredundant pathway
themes materially disappear from the exclusion sensitivity; individual genes
remain reported in the full results.
"""
    (ROOT / "reports/preanalytical_stress_sensitivity.md").write_text(
        report, encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "curated_genes": len(curated),
                "tested_genes": int(audit["present_in_tested_matrix"].sum()),
                "high_confidence_flagged": len(high_stress),
                "high_confidence_without_flags": len(without_stress),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
