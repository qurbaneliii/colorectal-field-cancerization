# Covariate-adjustment sensitivity

The primary biological comparison is adjacent-normal versus genuinely healthy
colon mucosa. U0 is unadjusted, U1 adjusts for continuous age and sex, and U2
additionally adjusts for left/right tumor location. All three designs were
full-rank; stage was not used because it is structurally undefined for healthy
donors.

| reference_model   | adjusted_model               |   shared_tested_genes |   pearson_log2fc_correlation |   spearman_log2fc_correlation |   sign_agreement_fraction |   reference_significant_genes |   adjusted_significant_genes |   significant_gene_retention_count |   significant_gene_retention_fraction |   median_effect_size_retention |   genes_losing_significance |   genes_gaining_significance |   genes_changing_direction | pathway_concordance                      |
|:------------------|:-----------------------------|----------------------:|-----------------------------:|------------------------------:|--------------------------:|------------------------------:|-----------------------------:|-----------------------------------:|--------------------------------------:|-------------------------------:|----------------------------:|-----------------------------:|---------------------------:|:-----------------------------------------|
| U0_unadjusted     | U1_age_sex_adjusted          |                 18490 |                     0.996376 |                      0.99207  |                  0.970092 |                          1662 |                         1580 |                               1549 |                              0.93201  |                       0.965816 |                         113 |                           31 |                        553 | evaluated after tier-specific enrichment |
| U0_unadjusted     | U2_age_sex_location_adjusted |                 18490 |                     0.995487 |                      0.989865 |                  0.964413 |                          1662 |                         1521 |                               1488 |                              0.895307 |                       0.947896 |                         174 |                           33 |                        658 | evaluated after tier-specific enrichment |

At FDR < 0.05 and |log2FC| >= 0.5, U0 detected
1,662 genes, U1 detected
1,580, and U2 detected
1,521. Covariate sensitivity is reported
as concordance and retention, not as proof that either model is causally
correct.

The evidence-tier rules were declared in `config/analysis.yaml`. Tier 1 requires
strong U1 evidence (FDR < 0.01, |adjusted log2FC| >= 1.0), U0/U1
direction agreement, U2 support, at least 80% sample-direction
consistency, deposited-matrix support, one-sample QC sensitivity robustness,
consistent direction in all adequately represented age/sex/location strata,
and a field-compatible trajectory. Machine-learning coefficients do not define
this biological tier.

- Tier 1 high-confidence genes: 101
- Tier 2 provisional genes: 1370
- Tier 3 exploratory genes: 7
