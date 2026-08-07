# Results

## Cohort and raw-CEL audit

All 246 expected GSE44076 CEL files were reconciled to 50 healthy donors and 98
complete tumor-adjacent patient pairs. All 390 expected GSE41258 CEL files were
reconciled to metadata and GPL96. The external primary analysis included
233 canonical patient-tissue arrays from 190
patients; the patient identifier, tissue role, and exclusions remain available
for every row.

## Raw-array quality control

All-array pre/post-RMA densities replaced the previous single-array diagnostic.
arrayQualityMetrics completed for both cohorts and produced
21 and
21 retained report files,
respectively. It flagged 8/246 and
48/390 arrays on at least one criterion.
No array met two severe custom failure criteria, so no primary exclusion changed.
GSE41258 NUSE medians ranged from 0.979
to 1.113. The single custom-QC
borderline GSE44076 array was retained; excluding it changed U1 significant genes
from 1,580 to
1,584, with
the full comparison in `results/tables/qc_exclusion_sensitivity.csv`.

## Differential expression and covariate sensitivity

At FDR <0.05 and absolute log2 fold change >=0.5, U0 identified
1,662 adjacent-versus-healthy
genes, U1 identified 1,580,
and U2 identified 1,521.
All designs were full rank. The patient-fixed-effect tumor-versus-adjacent
contrast identified 4,502
genes. These contrasts describe cross-sectional group differences and paired
tumor differences; they do not establish temporal progression.

Evidence integration yielded 101 Tier 1 high-confidence and
1,370 Tier 2 provisional field-associated genes. The complete
gene-level table preserves adjusted estimates, subgroup direction, preprocessing
agreement, QC sensitivity, trajectory, composition status, and evidence score.

## Tissue composition and stress sensitivity

74/101 (73.3%) high-confidence
genes remained composition robust and 27 were composition
sensitive. MCP-counter scores differed most strongly for fibroblasts
(adjusted adjacent-minus-healthy 1.376)
and endothelial cells (0.360),
with BH FDR values 4.93e-22
and 2.83e-25.
This attenuation indicates substantial microenvironmental contribution but does
not distinguish altered cell abundance from altered expression within cells.
The curated set contained 24 genes, 22
were measured, and 3 overlapped Tier 1; 98 Tier 1
genes remained after removal.

## Biological enrichment

High-confidence upregulated genes were enriched for extracellular matrix
organization, ECM degradation, proteoglycans, collagen degradation, and
regulation of insulin-like growth-factor transport. Composition-robust and
composition-sensitive subsets both retained ECM themes, supporting a mixed
epithelial/microenvironmental interpretation. Provisional downregulated genes
were enriched for fatty-acid, organic-acid, and small-molecule catabolism,
whereas provisional upregulated genes included angiogenesis and cell-substrate
adhesion. The three-gene Task C panel produced exploratory hyaluronan-related
terms; these small-set results are descriptive, not mechanistic evidence.

## Internal predictive performance

Task A was supporting. Elastic Net patient/donor-bootstrap aggregated OOF
macro-F1 was 0.959 (0.929-0.982) and multiclass ROC-AUC was
0.995 (0.988-0.999). Its conditional permutation null was
correctly restricted to tumor-versus-adjacent information and gave p=
0.0010.

Task B was primary. The nested panel policy selected a median of
5 genes and produced mean
outer-fold macro-F1 1.000 (fold SD
0.000). The final panel was
CLC, DYNC1H1, FOS, VIP, SNORA12. Aggregated repeated OOF ROC-AUC was
1.000 (1.000-1.000), macro-F1 1.000 (1.000-1.000),
Brier score 0.027 (0.023-0.031), and log loss
0.164 (0.153-0.178). The group-level permutation p value was
0.0010.
Four panel genes met strict stability; SNORA12
did not and is explicitly labeled exploratory.

Task B remained strong after age/sex/location residualization (macro-F1
0.977) and
in the age/sex-matched subset (1.000).
The demographics-only model reached macro-F1
0.626. By contrast,
composition-plus-demographic residualization reduced macro-F1 to
0.850
(ROC-AUC 0.939),
which is a central limitation and is consistent with a substantial stromal or
immune contribution to the classifier signal.

Task C was secondary. The nested policy selected a median of
3 genes and mean outer-fold
macro-F1 was 0.990 (fold SD
0.016). The final panel was
FOXQ1, CEMIP, ETV4; FOXQ1: strictly_stable_gene, selection frequency 1.000, sign consistency 1.000; CEMIP: strictly_stable_gene, selection frequency 1.000, sign consistency 1.000; ETV4: strictly_stable_gene, selection frequency 1.000, sign consistency 1.000. Aggregated repeated OOF ROC-AUC was
0.991 (0.963-1.000), Brier score 0.021 (0.014-0.033),
and paired-permutation p=
0.0010.

## Locked threshold and external evaluation

The GSE44076-only grouped OOF rank-transport analysis selected threshold
0.41 with internal balanced accuracy
0.995; GSE41258 labels were not accessed.
The common platform universe contained 12,039 genes, but only
CEMIP, ETV4 of the three panel genes were present on GPL96.
Accordingly, external results refer to the separately serialized two-gene
transport refit.

In the primary external 233-array, 190-patient
estimand, ROC-AUC was 0.983, PR-AUC
0.994, macro-F1 0.828,
balanced accuracy 0.779, sensitivity
1.000, specificity
0.558, Brier score
0.053, and log loss
0.195. Patient-cluster 95% bootstrap intervals
were 0.985 (0.957-1.000) for ROC-AUC,
0.831 (0.755-0.889) for macro-F1, and
0.561 (0.417-0.692) for specificity. Calibration
intercept was -2.066 and slope was
2.620. Thus ranking remained strong,
but the locked threshold produced incomplete specificity and non-ideal
calibration transport. Threshold 0.5 and one-array-per-patient results remain
labeled sensitivities rather than replacements for the prespecified primary
estimand.
