# Final scientific analysis report

## Executive summary

GSE44076 supports reproducible transcriptomic differences between cancer-free
healthy colon and tumor-adjacent histologically normal mucosa. Adjustment for
age, sex, location, preprocessing, one-sample QC sensitivity, and demographic
subgroups preserves most of the field signal. Composition adjustment retains
74/101 high-confidence genes, showing that the result is not
wholly explained by estimated cell mixture while also preventing a purely
epithelial interpretation. Task B is internally strong but has no independent
healthy-versus-adjacent validation. Task C is a separate tumor-transition model:
its cross-platform evaluation shows strong discrimination with incomplete
threshold and probability transport. No clinical-readiness claim is made.

## Scientific objectives

The primary objective is healthy versus adjacent-normal field cancerization.
Task A supports three-state tissue discrimination. Task C is secondary and
supports tumor-versus-normal-colon cross-platform evaluation. GSE41258 cannot
validate the primary field-effect endpoint because it lacks cancer-free healthy
versus tumor-adjacent sampling.

## Cohorts

GSE44076 contains 246 HG-U219 arrays: 50 healthy, 98 adjacent-normal, and 98
tumor arrays, including 98 complete tumor-adjacent pairs. GSE41258 contains 390
audited HG-U133A arrays. Its locked external cohort contains 233
canonical patient-tissue arrays from 190 patients, including
43 patients with both normal-colon and primary-tumor tissue.

## Data integrity

GSM identifiers are unique, every GSE44076 CEL maps to metadata, healthy donors
have distinct donor IDs, and all 98 tumor-adjacent pairs are complete. GSE41258
eligibility excludes non-colon tissues, metastases, polyps, unresolved tissues,
author-excluded arrays, and named technical-replicate variants. The canonical
external rule retains at most one array per patient and tissue; in this dataset
all 233 eligible arrays are already canonical.

## Raw preprocessing

GSE44076 was independently RMA-normalized with `oligo` and the HG-U219 design;
GSE41258 was independently RMA-normalized with `affy` after automatic HG-U133A
CDF confirmation. The resulting probe-set counts were 49,386 and 22,283,
respectively. Probe-to-gene mappings were filtered to unambiguous symbols and
multiple probe sets were aggregated by a deterministic median. Platforms were
never jointly normalized.

## Quality control

All-array pre/post-RMA densities, RLE, MA diagnostics, PCA, correlations,
hierarchical clustering, sample medians/IQRs, and non-finite checks were
executed. `arrayQualityMetrics` completed for both cohorts and affyPLM NUSE was
used only for GPL96. No GSE44076 array met the prespecified requirement of at
least two independent severe technical failures; biological PCA separation was
not an exclusion criterion.

## Covariate structure

Age and sex were associated with GSE44076 tissue class in the screening audit;
location was available and estimable. Stage was structurally inappropriate for
healthy-versus-adjacent adjustment and was excluded. U0 (`~ tissue`), U1
(`~ tissue + age + sex`), and U2 (`~ tissue + age + sex + location`) were all
full rank.

## Differential expression

At BH FDR <0.05 and absolute log2 fold change >=0.5, U0 identified
1,662 adjacent-versus-healthy
genes. The paired patient-fixed-effect tumor-versus-adjacent analysis identified
4,502 genes. These are
association estimates from bulk tissue, not causal or longitudinal effects.

## Adjusted field-effect analysis

U1 identified 1,580 genes
and U2 identified 1,521.
U0 versus U1 effect correlation was Pearson
0.996
and Spearman
0.992;
0.932
of U0 significant genes were retained. Per-gene effect differences, FDRs,
direction agreement, and significance retention are exported separately.

## Field-gene evidence tiers

Prespecified evidence integration produced 101 high-confidence and
1,370 provisional field-associated genes. The broad exploratory
screen is explicitly separate. Tiering uses adjusted statistical evidence,
effect size, sample-direction consistency, subgroup agreement, raw-versus-
deposited concordance, QC sensitivity, and trajectory; ML importance alone does
not define biological significance.

## Tissue-composition sensitivity

MCP-counter 1.2.0 was run on the independently normalized gene matrix using a
pinned human marker definition. Endothelial and fibroblast scores were higher
in adjacent-normal tissue, while several immune scores were lower. U3 adjusts
for age, sex, and two composition PCs; 74/101 Tier 1 genes
retained FDR, effect-size, and direction support. Task B composition sensitivity
also derives two PCs strictly within each outer-training fold. Because these
scores come from the same bulk expression data, this is an overadjustment-prone
sensitivity rather than causal cell-mixture decomposition.

## Stress-response sensitivity

3 high-confidence genes overlap the curated immediate-early/
stress list, and 98 high-confidence genes remain after flagging
them out. Genes were not automatically removed. Pathway summaries without the
flagged set preserve the dominant field themes, although preanalytical stress
cannot be excluded from this retrospective dataset.

## Functional enrichment

Adjusted and tier-specific GO Biological Process and Reactome analyses support
extracellular-matrix organization, collagen/proteoglycan turnover, stromal and
vascular remodeling, and more provisional immune/metabolic programs. Ranked
results include leading-edge membership. Redundant terms are pruned for summary
figures; full BH-adjusted results and tested-set provenance remain available.

## Task A

Task A is supporting three-state discrimination. Elastic Net aggregated OOF
macro-F1 was 0.959 (0.929-0.982); the fully nested compact-panel
policy produced mean outer-fold macro-F1
0.948.
It is not the study's primary endpoint.

## Task B

Task B is the primary predictive endpoint. Grouped repeated nested CV keeps all
feature filtering, selection, scaling, hyperparameter search, calibration, and
panel-size choice within the relevant training boundary. The nested-policy mean
outer-fold macro-F1 was 1.000;
aggregated OOF ROC-AUC was 1.000 (1.000-1.000). Near-perfect
performance is treated as a reason for additional confounding scrutiny, not as
proof that confounding is absent.

## Task B confounding analysis

The demographic-only baseline achieved macro-F1
0.626. The fixed
five-gene sensitivity was 0.977
after fold-local demographic residualization and
0.850
after demographics plus two fold-local composition PCs. The age/sex-matched
subset achieved 1.000.
These fixed-panel sensitivities characterize robustness; the unbiased primary
performance estimate remains the nested-policy outer result.

## Task B final signature

The final full-development Elastic Net panel is CLC, DYNC1H1, FOS, VIP, SNORA12.
Per-gene selection frequency, sign consistency, coefficient direction, median
coefficient, repeat coverage, and biological tier are reported in the signature
table. It is internally validated only, retrospective, and not clinically ready.

## Task C

Task C models tumor versus adjacent-normal in GSE44076. The fully nested policy
achieved mean outer-fold macro-F1
0.990. Task C remains secondary to
the field-effect question.

## Task C final signature

The locked three-gene panel is FOXQ1, CEMIP, ETV4. Each gene meets the
configured strict stability frequency and sign-consistency thresholds. Only
CEMIP, ETV4 occur on GPL96; FOXQ1 is absent and is not silently
imputed.

## Threshold locking

The transport threshold 0.41 maximized the
prespecified balanced-accuracy criterion on patient-aggregated repeated grouped
GSE44076 OOF rank probabilities. Ties were resolved toward 0.5 and then the
lower threshold. GSE41258 labels were never accessed for this choice; 0.5 is
reported only as a sensitivity.

## External cohort structure

The external primary set contains 233 arrays from
190 patients: 52
normal-colon and 181 primary-tumor
arrays. 43 patients contribute
both tissues. Patient-cluster bootstrap resamples patients and retains all
canonical tissue observations for each sampled patient.

## External validation

The exact three-gene primary model was not applied unchanged. The signature
specification, Elastic Net family, rank representation, hyperparameter strategy,
and threshold were locked in GSE44076; a distinct transport model was refitted
there using CEMIP and ETV4 and then evaluated without tuning on GSE41258. Point
estimates were ROC-AUC 0.983, PR-AUC
0.994, macro-F1
0.828, balanced accuracy
0.779, sensitivity
1.000, specificity
0.558, Brier score
0.053, and log loss
0.195.

## Calibration

External calibration intercept was
-2.066 and slope was
2.620. Perfect sensitivity coexists
with specificity 0.558; therefore the result
is strong discrimination with incomplete threshold/probability transport, not
clinical excellence.

## Permutation tests

Each task used 1,000 group-preserving permutations. Task A used its conditional
null, Task B the global donor-level label null, and Task C paired within-patient
exchangeability. All three attained p=
0.000999, the minimum
attainable value 0.000999;
this is reported as finite Monte Carlo resolution, not p=0.

## Uncertainty

Mean outer-fold scores, repeat variability, aggregated repeated OOF scores,
group-bootstrap intervals, and external patient-cluster intervals are distinct
estimands. Task B aggregated OOF macro-F1 was
1.000 (1.000-1.000). External ROC-AUC was
0.985 (0.957-1.000) and macro-F1 was
0.831 (0.755-0.889).

## Raw-vs-deposited sensitivity

For adjacent versus healthy, raw-CEL and deposited-matrix log2 fold changes had
Pearson correlation
0.998,
sign agreement
0.994,
and 1466
overlapping significant genes. Raw-CEL analysis remains primary; the deposited
matrix is a sensitivity, not a substitute provenance route.

## Biological interpretation

The most defensible interpretation is a bulk-tissue field-associated program
with prominent ECM, stromal, and vascular components plus provisional immune,
metabolic, epithelial, and stress-response contributions. The evidence cannot
separate cell abundance from within-cell regulation or establish a longitudinal
healthy-to-adjacent-to-tumor mechanism. Biological and predictive gene sets are
compared explicitly rather than conflated.

## Limitations

The study is retrospective and single-cohort for Task B; healthy and adjacent
samples differ demographically; bulk expression is composition-sensitive;
processing metadata are incomplete; no independent field-effect cohort exists;
one Task C gene is absent on GPL96; external threshold/calibration transport is
imperfect; and no prospective clinical-utility analysis is available. Docker
execution is blocked by the unavailable host Linux-container engine.

## Final scientific conclusions

The primary field-effect evidence is strong within GSE44076 and survives several
technical and demographic checks, but its marked composition sensitivity and
lack of independent Task B validation require a cohort-specific, noncausal
interpretation. Task B is ready for internal research use, not clinical use.
Task C supplies strong cross-platform ranking evidence but only moderate
threshold/probability transport. Publication readiness remains partial because
current-revision remote CI and container execution are blocked in this host
session and administrative manuscript metadata remain to be supplied; clinical
readiness is not claimed.
