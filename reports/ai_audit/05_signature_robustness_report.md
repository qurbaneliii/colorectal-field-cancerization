# Signature robustness

`results/tables/task_b_final_signature.csv` distinguishes nonzero selection in
the wider Elastic Net ranking model from inclusion in the final compact panel
across outer folds. CLC, DYNC1H1, FOS and VIP were present in every wider-model
selection. The corresponding compact counts were 15, 7, 12 and 14 out of 15.
SNORA12 was included in 4 of 15 compact panels and 8 of 15 wider models.
The five-gene full-development panel is a fitted model specification.

The GEO-deposited GSE44076 sensitivity reselected up to five genes separately
in each outer training fold with the same prespecified Elastic Net parameters.
CLC, DYNC1H1, FOS and VIP each occurred in 15 of 15 unadjusted folds. After
removing the stress and MCP-counter marker lists and residualizing on
demographics and two training-fitted marker-score PCs, CLC and LY96 occurred in
15 folds, VIP in 14 and DYNC1H1 in 11. Frequencies reflect overlapping training
partitions in the same donors. See `results/ai_audit/task_b_development_stability.csv`.

The historical biological tier contains 101 high-confidence field-associated
genes. Seventy-four retained composition-adjusted support and three were flagged
by the curated stress list. These counts are properties of the historical tier
definition and marker panel. FOS has a credible ischemia response. Even after
excluding curated genes, correlated transcripts can still carry the same
acquisition-related information.

The Task C compact model contains FOXQ1, CEMIP and ETV4. Archived compact fold
counts are 12, 11 and 10 of 15. The GPL96 mapping excludes FOXQ1. Only the
CEMIP/ETV4 refit was tested in the reported cross-platform result. Similar
within-study selection frequencies do not guarantee concordance after platform
change.
