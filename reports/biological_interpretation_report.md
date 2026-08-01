# Biological interpretation report

## Evidence boundary

The primary biological contrast is adjacent-normal mucosa from cancer patients
versus colon mucosa from cancer-free donors. It produced 101 Tier 1 and
1,377 provisional candidate field-associated genes after separate
assessment of adjustment, sample consistency, preprocessing, QC, subgroup, and
trajectory evidence. These associations do not establish causation.

## Nonredundant themes

- Extracellular matrix and stromal remodeling: collagen turnover, ECM
  proteoglycans, degradation, cell-substrate adhesion, and IGF transport dominate
  the high-confidence upregulated signal.
- Vascular and angiogenic signaling: adjacent-normal tissue has higher
  endothelial and fibroblast MCP-counter scores and provisional upregulated genes
  include regulation of angiogenesis.
- Immune context: several immune population scores differ between groups, but
  bulk expression cannot distinguish abundance from within-cell activation.
- Metabolism: fatty-acid, carboxylic-acid, and small-molecule catabolism are
  concentrated in the provisional downregulated tier and are less robust than
  the main ECM signal.
- Composition: 74/101 Tier 1 genes retain U3 support;
  27 attenuate, so the field signature should not be described as
  purely epithelial.
- Immediate-early/stress response: 3 Tier 1 genes overlap the
  curated set and 98 remain after removal; the dominant ECM themes
  persist, but preanalytical contribution cannot be excluded.
- Predictive versus biological genes: the Task B and Task C panels optimize
  training-only prediction under nested CV and are not substitutes for the
  evidence-tier gene set. Task C hyaluronan terms are exploratory because the
  panel contains only three genes.

Every enrichment row records whether it is over-representation or ranked,
direction, database, tested universe, mapped count, raw p value, and BH FDR in
`results/tables/supplementary_enrichment_results.csv`.
