# Historical external experiment

`external_validation_19_gene_experiment.json` and
`external_validation_19_gene_experiment_metrics.csv` preserve the JSON and
Table 5 from the initial Git revision `62cce06`. At the audit baseline
`c2a36a1d2f80fd25e1565b543c909eebef64ca09`, these files were still present as
`reports/external_validation_report.json` and
`results/tables/table_5_external_validation.csv`.

The figure trio `GSE41258_external_validation_performance.{pdf,png,svg}`
also depicts the 14-common-gene analysis, as shown by the SVG title and the
51/1/1/180 confusion counts. Git traces the figure to `a5a53d7`. The active
figure is `results/figures/external_validation_primary.svg`.

The original experiment used 19 development genes, 14 genes in common with
GPL96 and a threshold of 0.50. Its metrics use the earlier label mapping and
evaluation definition. The active report describes the later three-gene
development model and two-gene transport refit.
