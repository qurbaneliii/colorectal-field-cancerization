# Abstract

## Background

Histologically non-neoplastic mucosa adjacent to colorectal cancer may differ
from colon mucosa of cancer-free donors, but molecular field-effect studies are
vulnerable to demographic imbalance, tissue-composition confounding, patient
leakage, and overstatement of tumor-normal transfer as field validation.

## Methods

We reprocessed raw CEL files from GSE44076 (246
arrays: 50 healthy, 98 adjacent-normal, and 98 paired tumor) by robust multi-array
average and analyzed adjacent-normal versus healthy expression with unadjusted,
age/sex-adjusted, and age/sex/location-adjusted limma models. Candidate genes were
tiered by adjusted evidence, direction consistency, preprocessing concordance,
quality-control sensitivity, demographic strata, and expression trajectory.
MCP-counter and a curated immediate-early/stress set assessed microenvironmental
and preanalytical sensitivity. Elastic Net panels were evaluated with three
repeats of five-fold donor/patient-grouped nested cross-validation. A separate
tumor-versus-adjacent model was locked in GSE44076 and evaluated on GSE41258 with
patient-cluster bootstrap uncertainty.

## Results

Age/sex adjustment identified 1,580
genes at Benjamini-Hochberg FDR <0.05 and absolute log2 fold change >=0.5;
101 formed the high-confidence field signature and 1,377
were provisional. Of the high-confidence genes, 74/101
retained direction and significance after composition-PC adjustment and
3 overlapped the curated stress set. The primary Task B panel was
CLC, DYNC1H1, FOS, VIP, SNORA12. Its aggregated repeated out-of-fold ROC-AUC was
1.000 (1.000-1.000) and Brier score was
0.027 (0.023-0.031); there is no independent healthy-versus-adjacent
cohort. The secondary Task C panel was FOXQ1, CEMIP, ETV4. Only
CEMIP, ETV4 were represented on GPL96, so the external analysis
used a separately labeled common-gene rank-transport refit, not direct transport
of the exact serialized model. In 233 canonical patient-tissue
arrays from 190 patients, ROC-AUC was
0.985 (0.957-1.000), macro-F1 was
0.831 (0.755-0.889), and specificity was
0.561 (0.417-0.692) at the GSE44076-locked threshold.

## Conclusions

GSE44076 supports a reproducible, internally validated field-effect classifier
and a composition-sensitive candidate field signature. The GSE41258 analysis is
external evaluation of a separate tumor-normal transfer model, not validation of
Task B. These retrospective findings are biomarker-discovery evidence and do not
establish clinical utility or causal mechanisms.
