# Cross-platform harmonization gate

RNA-seq and U219 array measurements are different assays. The current microarray
scale and threshold cannot simply be applied to read counts. Stranded mRNA
selection may be particularly unsuitable for measuring the small nucleolar RNA
SNORA12. Public protocol alone does not establish its detection or quantification.

Before outcomes are evaluated, lock genome/annotation releases and their hashes,
reference FASTA/GTF, library orientation, quantifier/version, read-QC criteria,
gene identifiers, duplicate/synonym policy, missing-gene behavior and the rank
universe. Distinguish missing annotation, low detection and an actual zero count.
Require the full five-gene panel for an unchanged Task B model claim. Do not impute
an absent SNORA12 from external outcomes or silently drop it.

Any alternative panel or within-sample rank model must first be developed and
locked on development data using fold-local processing. Test assay bridging
without external labels where suitable reference material exists. Numeric model,
mapping and normalization fingerprints then feed `evaluate_frozen_cohort.py`.
Its JSON coefficients require compatible normalized inputs; a string claiming
`frozen_single_sample` is not scientific proof of assay portability.

Refitting on JGA labels constitutes adaptation/new development and requires a
new independent evaluation. The current proposed validation remains blocked by
access, eligibility and assay-transport evidence. No compatible RNA-seq model has
been validated or substituted for the historical model.
