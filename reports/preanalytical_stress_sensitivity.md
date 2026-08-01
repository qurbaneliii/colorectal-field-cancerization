# Preanalytical stress and immediate-early sensitivity

The prespecified flag set contains 24 immediate-early, AP-1,
feedback, and acute-response genes documented from the curated immediate-early
response analysis of Tullai et al. (Journal of Biological Chemistry 2007;
282:23981-23995; DOI: 10.1074/jbc.M703225200). The set is versioned at
`config/immediate_early_stress_genes.tsv`. These genes were flagged rather
than automatically removed because acute signaling can reflect tissue handling,
true mucosal biology, or both.

- Curated genes present in the tested matrix: 22
- Flagged genes in the 101-gene high-confidence field set: 3
- Prevalence in the high-confidence field set: 3.0%
- Flagged genes among the top 100 high-confidence candidates: 3
- High-confidence set after flag-only removal sensitivity: 98 genes

Flagged high-confidence genes: JUN, NR4A1, NR4A2.

Tier-specific enrichment is repeated both with and without these flagged genes.
Conclusions are treated as stress-dependent only when the nonredundant pathway
themes materially disappear from the exclusion sensitivity; individual genes
remain reported in the full results.
