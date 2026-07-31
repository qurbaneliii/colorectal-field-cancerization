# Raw-CEL versus processed-matrix sensitivity report

Raw-CEL RMA limma is the primary biological analysis. GEO deposited-matrix
Welch and paired-t results are retained only as a sensitivity analysis.

| comparison               |   shared_tested_genes |   pearson_log2fc_correlation |   spearman_rank_correlation |   sign_agreement_fraction |   top_100_gene_overlap |   raw_significant_genes |   deposited_significant_genes |   significant_gene_overlap |   materially_changed_genes | pathway_agreement                                    |
|:-------------------------|----------------------:|-----------------------------:|----------------------------:|--------------------------:|-----------------------:|------------------------:|------------------------------:|---------------------------:|---------------------------:|:-----------------------------------------------------|
| adjacent_vs_healthy      |                 15879 |                     0.99791  |                    0.995523 |                  0.99408  |                     67 |                    1478 |                          1483 |                       1466 |                        123 | reported after raw enrichment; no mechanism inferred |
| tumor_vs_healthy         |                 15879 |                     0.997807 |                    0.997295 |                  0.995592 |                     59 |                    4403 |                          4422 |                       4384 |                        124 | reported after raw enrichment; no mechanism inferred |
| tumor_vs_adjacent_paired |                 15879 |                     0.99771  |                    0.997048 |                  0.994521 |                     99 |                    3958 |                          3980 |                       3943 |                        135 | reported after raw enrichment; no mechanism inferred |

Field-cancerization candidates additionally required sample-direction
consistency and deposited-matrix sign agreement. Concordance does not establish
causality, and conclusions that materially change are listed in
`results/tables/raw_vs_deposited_material_changes.csv`.
