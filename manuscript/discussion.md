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
