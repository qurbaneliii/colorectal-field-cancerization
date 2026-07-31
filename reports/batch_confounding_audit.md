# Batch and confounding audit

Only deposited covariates were used; absent fields were not invented. Screening
associations do not trigger automatic ComBat or sample exclusion.

## Missingness and usability

| accession   | covariate         |   included_samples |   available_values |   missing_values |   unique_nonmissing_values | usable_for_association   | source                  |
|:------------|:------------------|-------------------:|-------------------:|-----------------:|---------------------------:|:-------------------------|:------------------------|
| GSE44076    | batch             |                246 |                  0 |              246 |                          0 | False                    | deposited metadata only |
| GSE44076    | processing_date   |                246 |                  0 |              246 |                          0 | False                    | deposited metadata only |
| GSE44076    | scan_date         |                246 |                  0 |              246 |                          0 | False                    | deposited metadata only |
| GSE44076    | sex               |                246 |                246 |                0 |                          2 | True                     | deposited metadata only |
| GSE44076    | age               |                246 |                246 |                0 |                         45 | True                     | deposited metadata only |
| GSE44076    | location          |                246 |                246 |                0 |                          2 | True                     | deposited metadata only |
| GSE44076    | stage             |                246 |                246 |                0 |                          3 | True                     | deposited metadata only |
| GSE44076    | msi_status        |                246 |                  0 |              246 |                          0 | False                    | deposited metadata only |
| GSE44076    | molecular_subtype |                246 |                  0 |              246 |                          0 | False                    | deposited metadata only |
| GSE44076    | center            |                246 |                  0 |              246 |                          0 | False                    | deposited metadata only |
| GSE41258    | batch             |                233 |                  0 |              233 |                          0 | False                    | deposited metadata only |
| GSE41258    | processing_date   |                233 |                  0 |              233 |                          0 | False                    | deposited metadata only |
| GSE41258    | scan_date         |                233 |                  0 |              233 |                          0 | False                    | deposited metadata only |
| GSE41258    | sex               |                233 |                224 |                9 |                          2 | True                     | deposited metadata only |
| GSE41258    | age               |                233 |                224 |                9 |                         53 | True                     | deposited metadata only |
| GSE41258    | location          |                233 |                224 |                9 |                          9 | True                     | deposited metadata only |
| GSE41258    | stage             |                233 |                224 |                9 |                          4 | True                     | deposited metadata only |
| GSE41258    | msi_status        |                233 |                208 |               25 |                          3 | True                     | deposited metadata only |
| GSE41258    | molecular_subtype |                233 |                  0 |              233 |                          0 | False                    | deposited metadata only |
| GSE41258    | center            |                233 |                  0 |              233 |                          0 | False                    | deposited metadata only |

## Tissue-label association screens

| accession   | covariate   | method         |   statistic |    p_value |   samples | interpretation                                            |
|:------------|:------------|:---------------|------------:|-----------:|----------:|:----------------------------------------------------------|
| GSE44076    | sex         | chi-square     |   6.30563   | 0.0427317  |       246 | screening association only; no automatic batch correction |
| GSE44076    | age         | Kruskal-Wallis |  13.7302    | 0.00104358 |       246 | screening association only; no automatic batch correction |
| GSE44076    | location    | chi-square     |   3.79379   | 0.150034   |       246 | screening association only; no automatic batch correction |
| GSE44076    | stage       | chi-square     | 246         | 4.7337e-52 |       246 | screening association only; no automatic batch correction |
| GSE41258    | sex         | chi-square     |   0         | 1          |       224 | screening association only; no automatic batch correction |
| GSE41258    | age         | Kruskal-Wallis |   0.0226799 | 0.880292   |       224 | screening association only; no automatic batch correction |
| GSE41258    | location    | chi-square     |   8.30964   | 0.403823   |       224 | screening association only; no automatic batch correction |
| GSE41258    | stage       | chi-square     |   1.4908    | 0.684395   |       224 | screening association only; no automatic batch correction |
| GSE41258    | msi_status  | chi-square     |   1.19373   | 0.550534   |       208 | screening association only; no automatic batch correction |

Any estimable covariate must be incorporated through a prespecified inferential
design or learned inside training folds. GSE41258 labels must never be used to
harmonize the external cohort with GSE44076.
