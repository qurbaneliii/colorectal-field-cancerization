# Identification of transcriptomic field-cancerization signatures in colorectal cancer using explainable and leakage-safe machine learning

[AUTHOR NAME]

[AFFILIATION]

Corresponding author: [AUTHOR NAME]

Article type: Original research

Running title: Colorectal transcriptomic field cancerization

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
101 formed the high-confidence field signature and 1,370
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

# Keywords

colorectal cancer; field cancerization; adjacent-normal mucosa; microarray;
nested cross-validation; Elastic Net; tissue composition; external validation

# Introduction

Field cancerization describes molecular alteration extending beyond a visible
neoplasm into apparently non-neoplastic tissue. Genomic and transcriptomic
differences have been demonstrated in colorectal mucosa surrounding tumors
[1], and GSE44076 was specifically generated to compare healthy-donor mucosa,
tumor-adjacent mucosa, and matched colorectal tumors [2]. Tumor-adjacent normal
tissue is therefore a biologically informative compartment rather than a
universally exchangeable healthy control [3]. Pan-cancer analyses likewise show
that adjacent tissue can contain inflammation, stromal response, and other
signals distinct from both healthy tissue and tumor [4].

Three questions must remain separate. The primary biological and predictive
question is whether healthy-donor mucosa differs from adjacent-normal mucosa
(Task B). A supporting three-class analysis asks whether healthy, adjacent, and
tumor tissue states can be discriminated (Task A), without interpreting group
means as a longitudinal sequence. A secondary question asks whether a
GSE44076-derived tumor-versus-adjacent classifier can transfer across Affymetrix
platforms to GSE41258 tumor and normal-colon specimens (Task C). GSE41258 cannot
independently validate healthy-versus-adjacent field cancerization because it
does not contain the required healthy and adjacent-normal endpoint definition.

Microarray biomarker discovery creates several opportunities for optimism:
normalization may be dataset-specific, multiple samples from one patient may
cross validation folds, feature and panel selection may use test outcomes, and
thresholds may be tuned on external labels. Nested cross-validation is required
when tuning is part of the model-building procedure [10], and grouped resampling
is necessary when observations share a patient or donor [11]. This study was
therefore designed around raw-data provenance, prespecified covariate and
composition sensitivity, training-only transformations, grouped nested
validation, task-specific permutation tests, and explicit separation of the
primary model from the cross-platform transport refit.

We aimed to identify candidate field-cancerization genes with graded evidence;
quantify demographic, tissue-composition, stress, and QC sensitivity; derive an
interpretable internally validated Task B classifier; and evaluate a distinct
Task C tumor-normal signal across platforms while preserving honest boundaries
on calibration, transport, and clinical readiness.

# Methods

## Study design and datasets

This retrospective secondary analysis used public Affymetrix data. GSE44076 on
GPL13667 contained 246 audited arrays: 50 healthy-donor mucosa, 98
adjacent-normal mucosa, and 98 tumors from the same 98 cancer patients [2].
Healthy donors received unique donor IDs; paired tumor and adjacent samples
shared patient IDs. GSE41258 on GPL96 contained 390 deposited arrays [3]. The
external endpoint retained one canonical sample per patient and tissue after
excluding non-primary lesions and designated technical replicates, yielding
233 arrays from 190 patients (52
normal-colon and 181 primary-tumor
arrays); 43 patients contributed both tissues. No cross-platform
normalization or pooled fitting was performed.

## Raw-data audit, preprocessing, and quality control

CEL membership was reconciled to GEO accessions and platform CDFs before
analysis. GSE44076 was processed with the `oligo` route appropriate to HG-U219;
GSE41258 used the `affy` route appropriate to HG-U133A. Each cohort was
background corrected, quantile normalized, and summarized independently using
robust multi-array average [5]. Probe sets were mapped to gene symbols using
platform-specific annotation; ambiguous mappings and collapse rules are retained
in supplementary mapping tables.

Quality assessment comprised pre/post-RMA all-array density distributions,
boxplots, relative log expression, MA diagnostics, PCA, sample correlations,
hierarchical clustering, robust sample median/IQR diagnostics, and executed
arrayQualityMetrics reports [14]. For GPL96, probe-level models supplied NUSE;
no surrogate NUSE was generated for the oligo platform. A sample could be
excluded only after at least two independent severe custom technical failures;
PCA separation alone was never sufficient. Borderline arrays remained in the
primary analysis and were removed only in a labeled sensitivity analysis.

## Differential expression and field-gene evidence

Gene-wise linear models were fit with limma empirical Bayes moderation [6]. For
adjacent-normal versus healthy tissue, U0 included tissue only, U1 included
tissue, continuous age, and sex, and U2 additionally included left/right
location. Stage was not added because it is structurally undefined for healthy
donors. The paired tumor-versus-adjacent model used patient fixed effects and a
tissue coefficient. All design matrices and coefficient names were audited for
full rank. P values were adjusted by the Benjamini-Hochberg procedure [7]; the
main reporting threshold was FDR <0.05 and absolute log2 fold change >=0.5.

Tier 1 high-confidence field genes required strong U1 evidence (FDR <0.01 and
absolute adjusted log2 fold change >=1), consistent U0/U1 direction, U2 support,
at least 80% sample-direction consistency, deposited-matrix support, QC-exclusion
robustness, consistent direction in adequately represented age/sex/location
strata, and a field-compatible three-group trajectory. Tier 2 genes were labeled
provisional. Model coefficients did not define biological evidence tiers.

## Composition, stress, and enrichment sensitivity

MCP-counter 1.2.0 estimated eight immune and two stromal population scores from
human gene symbols [15]. The U3 sensitivity model added two composition principal
components learned from these scores. Genes retaining significance and direction
were called composition robust; attenuation was reported rather than interpreted
as proof of cell-intrinsic or cell-extrinsic causality. A 24-gene literature-based
immediate-early/stress set [18] was audited against measured symbols and removed
in a sensitivity enrichment analysis.

GO Biological Process and Reactome analyses used `clusterProfiler` [16] and the
Reactome knowledgebase [17]. Over-representation analyses used the tested mapped
gene universe, direction-specific gene sets, and BH correction. Ranked analyses
used signed limma statistics. Redundant GO terms were reduced by semantic
similarity; every output records the set, direction, database, universe, and
analysis type.

## Predictive modeling and internal validation

Task B (healthy versus adjacent-normal) was the primary predictive endpoint;
Task C (tumor versus adjacent-normal) was secondary; Task A was supporting.
Elastic Net logistic regression [8] was prespecified as the interpretable model,
with calibrated linear SVM and random forest comparators. All scaling, variance
filtering, univariate selection, residualization, hyperparameter tuning, feature
ranking, and compact-panel selection occurred within training partitions.

Internal performance used three deterministic repetitions of five-fold outer
cross-validation with donor/patient groups kept intact. Four-fold grouped inner
cross-validation selected hyperparameters and panel size from candidates 3, 5,
10, 15, and 20 plus the fold-specific nonzero model. The smallest panel within
the configured macro-F1 tolerance, balanced-accuracy floor, and Brier-loss limit
was selected. Outer-fold means and variability describe model-building-policy
performance; patient/donor bootstrap intervals of sample-aggregated repeated OOF
predictions describe a separate estimand. Calibration was assessed with curves,
Brier score, log loss, intercept, and slope where estimable [12].

Stability used all repeat-fold opportunities in its denominator. A gene was
strictly stable only with selection frequency >=0.65 and sign consistency >=0.80;
compact-consensus and exploratory labels were reported separately. Null tests
were task-specific: Task A held healthy labels fixed while permuting the
tumor/adjacent component, Task B permuted group-level healthy/adjacent labels,
and Task C exchanged tumor/adjacent labels within patients. Each used 1,000
permutations.

## Locked Task C threshold and external evaluation

The Task C panel and probability threshold were derived without GSE41258 labels.
Within-sample percentile ranks were computed over the label-independent common
gene universe. Repeated grouped GSE44076 OOF probabilities were averaged by
patient and the threshold maximizing balanced accuracy was selected, with ties
resolved toward 0.5 and then the lower threshold. Because FOXQ1 was
absent from GPL96, the exact three-gene serialized primary model was not directly
transportable. A distinct two-gene refit using CEMIP, ETV4 and
GSE44076 outcomes was serialized and labeled as a transport refit. External
labels were used only for final evaluation. The primary external estimand used
all canonical patient-tissue rows and patient-cluster bootstrap intervals (1,000
replicates); one-array-per-patient and threshold 0.5 were sensitivities.

## Reproducibility and reporting

Python and R dependencies were locked in `requirements-lock.txt`/`pyproject.toml`
and `renv.lock`, respectively. Seeds, grids, thresholds, panel sizes, bootstrap
counts, and permutation counts are centralized in `config/analysis.yaml`.
Reporting follows TRIPOD principles where applicable [13], while recognizing
that this is retrospective molecular biomarker discovery rather than a clinical
prediction-model validation study.

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

# Discussion

This analysis supports reproducible transcriptomic differences between
adjacent-normal colorectal mucosa and mucosa from cancer-free donors. The signal
was not an artifact of the unadjusted contrast: 1,580
genes remained after age/sex adjustment, effect estimates were highly concordant
after adding tumor location, and 101 genes met a deliberately stricter,
multicomponent high-confidence definition. These are candidate field-associated
genes, not proof that a tumor caused each alteration or that each alteration
pre-dated cancer.

The dominant enrichment pattern involved extracellular-matrix organization,
collagen and proteoglycan turnover, stromal activation, vascular development,
and angiogenic signaling. MCP-counter analysis materially qualifies that result:
fibroblast and endothelial scores were elevated in adjacent tissue, and only
74/101 Tier 1 genes retained composition-adjusted support.
Both robust and sensitive sets still showed ECM themes, suggesting that the
bulk-tissue signal plausibly combines epithelial transcription with altered
microenvironmental abundance or state. Bulk arrays cannot separate those
components. The reduced Task B performance after composition-plus-demographic
residualization reinforces this concern and makes the near-perfect unadjusted
classification unsuitable as evidence of a clinically transportable epithelial
biomarker.

Metabolic attenuation was concentrated in the broader provisional tier,
including fatty-acid and organic-acid catabolism, while strong Tier 1 evidence
was more ECM-centered. This difference argues against treating every significant
gene as equivalently robust. Immediate-early and stress-response genes were a
small but visible part of the Tier 1 list. Removing the curated set left
98 genes and preserved the main ECM interpretation, but sample
handling remains incompletely recorded and cannot be ruled out as a contributor.

The Task B panel (CLC, DYNC1H1, FOS, VIP, SNORA12) directly addresses the paper's
healthy-versus-adjacent question. Grouped nested validation, training-only panel
selection, calibration assessment, bootstrap uncertainty, and a group-level
permutation null reduce common sources of leakage. Nonetheless, very high
within-cohort performance can reflect cohort recruitment, demographics, tissue
handling, anatomic sampling, and composition in addition to biological field
effects. Residualization and matching support—but cannot prove—robustness.
Independent healthy and adjacent-normal cohorts processed under harmonized
protocols are required before generalization can be claimed.

Task C answers a different question. Its tumor-versus-adjacent ranking
transferred strongly to GSE41258 tumor versus normal-colon tissue, but FOXQ1 was
absent on GPL96. The reported external evaluation therefore belongs to a
two-gene rank-transport refit, not the exact three-gene primary artifact. The
contrast between ROC-AUC 0.983 and balanced
accuracy 0.779 illustrates why
threshold-free ranking and threshold-dependent classification must be reported
together. Perfect sensitivity with limited specificity at the locked threshold,
plus calibration intercept -2.066
and slope 2.620, indicates incomplete
threshold and calibration transport despite strong discrimination.

The study's strongest contribution is methodological alignment: biological
field genes, the primary Task B classifier, the supporting Task A model, and the
secondary Task C transfer model are no longer conflated. Prospective translation
would require independent field-effect cohorts, standardized sampling distance
and processing, cell-resolved validation, assay conversion, prespecified locked
models and thresholds, external calibration, clinical utility analysis, and
prospective evaluation in the intended population.

# Limitations

This study is retrospective and uses two public microarray cohorts. Healthy and
cancer-bearing participants were not randomized or prospectively matched, so
unmeasured demographic, clinical, anatomic, medication, and processing factors
may remain. Recorded covariates were incomplete, particularly in GSE41258.

Bulk tissue conflates epithelial, stromal, endothelial, and immune signals.
MCP-counter is a marker-based sensitivity method, not a cell-count gold standard,
and composition-PC residualization may remove biologically meaningful field
signals as well as confounding. The immediate-early/stress audit is limited by a
curated gene list and incomplete preanalytical metadata.

Task B has internal grouped nested validation only. No independent cohort with
both genuinely healthy donor mucosa and tumor-adjacent mucosa was available, so
there is no external field-effect validation. Very high internal discrimination
may still reflect cohort-specific structure.

Task C uses a different normal comparator in GSE41258 and one of three primary
genes is absent on GPL96. Its two-gene transport refit is not the exact serialized
primary model. Strong ROC-AUC coexists with imperfect locked-threshold specificity
and calibration. Patient-cluster bootstrap accounts for repeated tissues but
does not correct dataset shift.

The models are retrospective biomarker-discovery tools. They have not undergone
assay validation, decision-curve analysis, prospective calibration, clinical
impact testing, regulatory review, or evaluation for individual patient care.
No causal biological mechanism or longitudinal progression is inferred from
the three cross-sectional tissue groups.

# Conclusion

Raw-CEL reanalysis of GSE44076 identified a graded colorectal field-associated
transcriptomic signal and a five-gene, internally validated healthy-versus-adjacent
classifier. Composition sensitivity shows that stromal and vascular context is
part of that signal. A separate three-gene tumor-versus-adjacent model supported
cross-platform tumor-normal ranking in GSE41258 through a clearly labeled
two-gene transport refit, with incomplete threshold and calibration transport.
The repository supports retrospective biomarker discovery and transparent
reproduction; it does not support an externally validated clinical
field-cancerization diagnostic.

# Data and code availability

GSE44076 and GSE41258 expression data and metadata are publicly available from
NCBI Gene Expression Omnibus at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44076
and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41258. Analysis code,
configuration, machine-readable provenance, result tables, model cards, and
manuscript sources are available at
https://github.com/qurbaneliii/colorectal-field-cancerization. Multi-gigabyte raw
archives and generated intermediate matrices are intentionally excluded from
Git; their expected checksums and reconstruction paths are documented in the
metadata manifests. No restricted patient-level clinical data were used.

# Author contributions

[AUTHOR NAME]: conceptualization, methodology, software, formal analysis,
validation, visualization, data curation, writing—original draft, and
writing—review and editing.

# Funding

[FUNDING INFORMATION]

# Conflicts of interest

The author declares no conflicts of interest.

# Acknowledgments

The author thanks the investigators who generated and publicly deposited
GSE44076 and GSE41258 and the maintainers of GEO, Bioconductor, and the open-source
scientific Python ecosystem.

# References

1. Hawthorn L, Lan L, Mojica W. Evidence for field effect cancerization in colorectal cancer. Genomics. 2014;103(2-3):211-221. doi:10.1016/j.ygeno.2013.11.003. PMID:24316131.
2. Sanz-Pamplona R, Berenguer A, Cordero D, et al. Aberrant gene expression in mucosa adjacent to tumor reveals a molecular crosstalk in colon cancer. Molecular Cancer. 2014;13:46. doi:10.1186/1476-4598-13-46. PMID:24597571.
3. Sheffer M, Bacolod MD, Zuk O, et al. Association of survival and disease progression with chromosomal instability: a genomic exploration of colorectal cancer. Proceedings of the National Academy of Sciences USA. 2009;106(17):7131-7136. doi:10.1073/pnas.0902232106. PMID:19359472.
4. Aran D, Camarda R, Odegaard J, et al. Comprehensive analysis of normal adjacent to tumor transcriptomes. Nature Communications. 2017;8:1077. doi:10.1038/s41467-017-01027-z. PMID:29057876.
5. Irizarry RA, Hobbs B, Collin F, et al. Exploration, normalization, and summaries of high density oligonucleotide array probe level data. Nucleic Acids Research. 2003;31(4):e15. doi:10.1093/nar/gng015. PMID:12582260.
6. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research. 2015;43(7):e47. doi:10.1093/nar/gkv007. PMID:25605792.
7. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B. 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
8. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software. 2010;33(1):1-22. doi:10.18637/jss.v033.i01. PMID:20808728.
9. Pedregosa F, Varoquaux G, Gramfort A, et al. Scikit-learn: machine learning in Python. Journal of Machine Learning Research. 2011;12:2825-2830.
10. Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. BMC Bioinformatics. 2006;7:91. doi:10.1186/1471-2105-7-91. PMID:16504092.
11. Roberts DR, Bahn V, Ciuti S, et al. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography. 2017;40(8):913-929. doi:10.1111/ecog.02881.
12. Van Calster B, Nieboer D, Vergouwe Y, et al. A calibration hierarchy for risk models was defined: from utopia to empirical data. Journal of Clinical Epidemiology. 2016;74:167-176. doi:10.1016/j.jclinepi.2015.12.005. PMID:26772608.
13. Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent Reporting of a multivariable prediction model for Individual Prognosis Or Diagnosis (TRIPOD): the TRIPOD statement. Annals of Internal Medicine. 2015;162(1):55-63. doi:10.7326/M14-0697. PMID:25560714.
14. Kauffmann A, Gentleman R, Huber W. arrayQualityMetrics—a Bioconductor package for quality assessment of microarray data. Bioinformatics. 2009;25(3):415-416. doi:10.1093/bioinformatics/btn647. PMID:19106121.
15. Becht E, Giraldo NA, Lacroix L, et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology. 2016;17:218. doi:10.1186/s13059-016-1070-5. PMID:27765066.
16. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284-287. doi:10.1089/omi.2011.0118. PMID:22455463.
17. Jassal B, Matthews L, Viteri G, et al. The Reactome pathway knowledgebase. Nucleic Acids Research. 2020;48(D1):D498-D503. doi:10.1093/nar/gkz1031. PMID:31691815.
18. Tullai JW, Schaffer ME, Mullenbrock S, Sholder G, Kasif S, Cooper GM. Immediate-early and delayed primary response genes are distinct in function and genomic architecture. Journal of Biological Chemistry. 2007;282(33):23981-23995. doi:10.1074/jbc.M702044200. PMID:17575275.

# Figure legends

**Figure 1. Study objectives and dataset roles.** GSE44076 (246 HG-U219 arrays;
50 healthy, 98 adjacent-normal, 98 tumor) supports Task B field-effect discovery,
supporting Task A, and Task C development; GSE41258 supplies 233
eligible HG-U133A patient-tissue arrays for secondary tumor-normal evaluation.
Arrows distinguish inference, internal validation, and external transfer; no
error bars or statistical tests are used. Abbreviations: CV, cross-validation.

**Figure 2. Sample inclusion flow.** Audited GEO arrays, tissue classes, paired
structure, exclusions, and the 190-patient external estimand are
shown. Counts follow CEL/metadata reconciliation; no inferential test or interval
is displayed.

**Figure 3. Raw-array quality control.** Cohort-wide boxplot, relative log
expression, MA, PCA, correlation, clustering, arrayQualityMetrics, and GPL96 NUSE
summaries are based on independently RMA-preprocessed CEL files. Samples required
two severe independent custom failures for exclusion; PCA separation alone did
not qualify. NUSE, normalized unscaled standard error; RLE, relative log expression.

**Figure 4. Pre- and post-RMA signal distributions.** Semi-transparent density
curves represent every array in GSE44076 (n=246) and GSE41258 (n=390), before and
after cohort-specific RMA. Curves are descriptive; no hypothesis test or error
bar is used. RMA, robust multi-array average.

**Figure 5. PCA, correlation, and clustering.** Unsupervised views use RMA
gene-expression values with samples colored by audited tissue class. PCA was a
diagnostic and not an exclusion criterion; no confidence intervals are shown.

**Figure 6. Adjacent-normal versus healthy differential expression.** Volcano
plots show U0 tissue-only and U1 age/sex-adjusted limma estimates for 18,490 genes
from 98 adjacent-normal and 50 healthy samples. Dashed lines mark BH FDR 0.05 and
absolute log2 fold change 0.5; red points satisfy both. BH, Benjamini-Hochberg.

**Figure 7. Unadjusted-adjusted concordance.** Each point is one gene; axes show
U0 and U1 log2 fold changes from the same 148 GSE44076 samples. The diagonal is
identity. The reported Pearson coefficient is descriptive; no error bars are used.

**Figure 8. Paired tumor versus adjacent-normal differential expression.** limma
patient fixed effects compare 98 complete pairs (196 arrays). Thresholds are BH
FDR 0.05 and absolute log2 fold change 0.5; no longitudinal interpretation is made.

**Figure 9. High-confidence field signature and gene distributions.** Heatmap
and representative distributions use GSE44076 RMA expression across 50 healthy,
98 adjacent-normal, and 98 tumor arrays. Tier 1 status comes from the prespecified
multi-evidence rule, not clustering; group summaries are cross-sectional.

**Figure 10. Tissue composition and evidence tiers.** MCP-counter 1.2.0 scores
for eight immune and two stromal populations are compared between 50 healthy and
98 adjacent-normal samples by the recorded two-group test with BH correction.
Boxplots show medians and interquartile ranges; whiskers follow the standard
1.5-IQR definition. Tier bars count genes after U0/U1/U2, QC, subgroup,
preprocessing, trajectory, and composition checks.

**Figure 11. Pathway enrichment.** Nonredundant GO Biological Process and
Reactome terms are shown for high-confidence and composition-stratified field
sets. Bar length is -log10 BH FDR. Over-representation used the tested mapped
gene universe; ranked panels are explicitly labeled separately.

**Figure 12. Supporting Task A validation.** Confusion counts are aggregated
repeated OOF predictions from three repeats of five-fold donor/patient-grouped
nested CV on 246 GSE44076 samples. All transformations and tuning were training
only. Bootstrap intervals are reported in the text, not as matrix error bars.

**Figure 13. Primary Task B discrimination and calibration.** ROC, precision-
recall, and calibration curves use aggregated repeated OOF predictions from 50
healthy and 98 adjacent-normal samples under grouped nested CV. The positive
class is adjacent-normal; confidence intervals use 1,000 donor/patient bootstrap
replicates. PR, precision-recall; ROC, receiver operating characteristic.

**Figure 14. Task B compact-panel selection.** Points show mean untouched
outer-fold macro-F1 and bars show outer-fold standard deviation for candidate
panel sizes. Panel size was selected only inside grouped inner CV using the
configured performance and Brier-loss tolerances.

**Figure 15. Secondary Task C discrimination, calibration, and compact-panel
selection.** Repeated grouped nested CV compares 98 tumor and 98 paired
adjacent-normal samples. ROC, PR, calibration, and panel-size curves were produced
without external labels; intervals use patient bootstrap where reported.

**Figure 16. Feature stability and threshold lock.** Feature frequency divides
by all 15 repeat-fold opportunities; strict stability requires frequency >=0.65
and sign consistency >=0.80. The threshold curve uses patient-aggregated repeated
GSE44076 OOF rank probabilities and locks 0.41;
GSE41258 labels were unavailable to selection.

**Figure 17. External tumor-normal evaluation.** ROC, PR, calibration, and
confusion matrix use 233 canonical patient-tissue GSE41258 arrays
from 190 patients, within-sample percentile ranks, the
GSE44076-locked threshold, and the two-gene common-platform transport refit.
Intervals are 95% patient-cluster bootstrap intervals from 1,000 replicates.

**Figure 18. External patient-structure and signature-intersection sensitivity.**
The primary all-canonical-tissue estimand is contrasted with one-array-per-patient
and threshold-0.5 sensitivities. The intersection identifies CEMIP, ETV4
as present on GPL96 and FOXQ1 as missing. No external result influenced
gene, representation, panel-size, or threshold selection.

# Supplementary materials

Supplementary data include complete sample metadata; exclusion logs; CEL archive
and platform audits; probe-to-gene mappings; full U0/U1/U2/QC and paired limma
results; evidence tiers; composition and stress sensitivity; GO/Reactome and
ranked enrichment outputs; fold assignments; repeated OOF predictions;
hyperparameters; feature coefficients and stability; nested panel-selection
decisions; task-specific permutation tests; threshold provenance; external
patient-tissue structure, predictions, and patient-cluster bootstrap intervals;
model cards; and a SHA-256 artifact manifest. Large raw and intermediate files
are reconstructed by the documented workflow and are not committed.
