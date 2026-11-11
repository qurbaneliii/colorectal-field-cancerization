# Statistical genomics audit

## F11 | Differential expression estimands

`R/04_differential_expression.R:67–184` implements limma models with BH
adjustment within each contrast and model. Unpaired comparisons use separate
groups of 50 and 98, with age and sex in U1 and location added in U2. The paired
tumor-versus-distant-mucosa design includes patient fixed effects and a tissue
term. The tracked paired design is rank 99 of 99 on 196 arrays. Static review and
the tracked design audit support the paired estimand. The case-control estimand
is a contrast between sampling protocols as well as between populations.

`results/tables/covariate_adjusted_de_concordance.csv` reports 18,490 tested
genes in the U0/U1/U2 comparisons. At FDR < 0.05 and absolute log2 fold change at
least 0.5, U1 retained 1,549 of 1,662 U0 hits and U2 retained 1,488. Effect
correlations with U0 were 0.996 and 0.995, respectively. Adjustment for measured
demographics therefore preserves much of the original association. U1/U2 design
condition numbers are approximately 401 and 371 on the uncentered age scale.
Centering age is a numerical sensitivity. It leaves the fitted model space
unchanged. The documented full-rank matrices do not establish overlap in
unobserved handling variables. These are archived counts from the primary route.

## Deposited-expression limma replication

`scripts/audit_deposited_limma.R` reran U0, U1, U2 and the patient-fixed-effect
contrast on the verified GSE44076 deposited matrix. It used the same 18,490
gene symbols, pinned probe annotation and limma 3.64.3, and centered age
in adjusted models. Pearson effect correlations with the archived raw-CEL
results exceeded 0.99998 in all four contrasts. For U1, 1,578 of the 1,580
archived hits recurred with the same direction. For the paired contrast,
4,496 of 4,502 recurred. Both routes used the same arrays. This supports
agreement across preprocessing sources within GSE44076.

When limma formally tested an absolute effect larger than 0.5 with
`treat`, U1 retained 537 genes at BH FDR below 0.05. The corresponding
paired contrast retained 3,184. These are distinct from testing against a
zero effect followed by an estimated effect-size filter. The data and
sample acquisition still limit interpretation. Concordance and panel
effects are in `results/ai_audit/deposited_limma_concordance.csv` and
`results/ai_audit/deposited_limma_panel_effects.csv`.

## F12 | Interpretation of stability and enrichment

Separate the FDR and effect-size thresholds from the biological evidence tiers.
Each limma BH result refers to its own tested contrast. Tier membership combines
selection, significance, magnitude and direction, so it has no standalone
experiment-wide FDR. Likewise, GO and Reactome enrichments are adjusted within
calls, not across the collection of alternative tier and sensitivity analyses.
`R/05_functional_enrichment.R:159–167` can explicitly combine old
over-representation artifacts with a new ranked run. This is a reproducibility
risk unless both inputs have compatible provenance. `R/05_functional_enrichment.R:232–245`
collapses duplicate Entrez IDs by the largest absolute moderated statistic,
which can favor multiply mapped genes. Quantify the difference against a
predeclared representative-mapping rule before treating marginal pathways as stable.

## F13 | Composition and repeated measures

`R/06_tissue_composition_sensitivity.R:45–64,126–151` uses marker-expression
scores and two derived PCs. The study lacks measured cell fractions. Estimated
scores can track correlated cell abundance, epithelial state, RNA handling or
downstream response to disease. A composition-adjusted coefficient is a
conditional association. Loss of statistical support does not apportion these
mechanisms. Matched tumor/distant comparisons absorb patient-level covariates,
but cannot establish a trajectory or eliminate tissue-specific handling effects.

The next analysis should prespecify contrasts and tested universes, retain
effect estimates with uncertainty, test probe aggregation and verify stability
with participant-level resampling and independent site-matched mucosa. Pathway
claims require independent direction and gene-set support.
See [Ritchie et al., limma](https://doi.org/10.1093/nar/gkv007) and
[Brazma et al., MIAME](https://doi.org/10.1038/ng1201-365).
