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
