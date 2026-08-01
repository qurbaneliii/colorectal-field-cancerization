# QC exclusion sensitivity

No GSE44076 array met the prespecified custom-QC exclusion rule of at least two
independent severe technical failures. `GSM1077736` had one severe custom
IQR flag and was retained in the primary analysis; the sensitivity analysis
excluded it. arrayQualityMetrics flags remained review signals and did not
change the exclusion decision. PCA separation was never an exclusion rule.

| analysis_component                      |   primary_value |   sensitivity_value | comparison                                                                                  |
|:----------------------------------------|----------------:|--------------------:|:--------------------------------------------------------------------------------------------|
| borderline_samples                      |               0 |            1        | GSM1077736                                                                                  |
| significant_DE_genes                    |            1580 |         1584        | FDR<0.05 and |log2FC|>=0.5                                                                  |
| top_100_gene_overlap                    |             100 |           97        | ranked by adjusted p-value                                                                  |
| log2FC_pearson                          |               1 |            0.999881 | all shared tested genes                                                                     |
| top_20_pathway_term_overlap             |               5 |            1        | primary ranked U1 versus QC-excluded over-representation; extracellular matrix organization |
| task_b_aggregated_oof_macro_f1          |               1 |            1        | five-gene panel; repeated grouped OOF                                                       |
| task_b_aggregated_oof_balanced_accuracy |               1 |            1        | five-gene panel; repeated grouped OOF                                                       |
| task_b_reselected_panel_overlap         |               5 |            4        | shared=CLC;DYNC1H1;FOS;VIP; sensitivity=APOLD1;CLC;DYNC1H1;FOS;VIP                          |
| top_ranked_pathway_term_overlap         |               5 |            1        | primary ranked U1 versus QC-excluded over-representation; extracellular matrix organization |

The sensitivity reproduces the age/sex-adjusted differential-expression model,
the locked Task B classifier under repeated donor/patient-grouped out-of-fold
evaluation, and a separate full-development grouped-inner-CV panel selection.
Pathway overlap compares up to the top 20 exact descriptions available from
ranked U1 GO-BP analysis with the union of the QC-excluded up/down
over-representation results;
the different enrichment estimands are retained explicitly rather than treated
as interchangeable. Complete gene membership is in
`results/tables/qc_exclusion_model_gene_comparison.csv`.
