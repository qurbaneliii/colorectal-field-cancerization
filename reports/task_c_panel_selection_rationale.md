# Task C panel selection rationale

Task C panel size was treated as an inner-CV hyperparameter. For each outer
fold, rankings were constructed from inner-training data, candidate sizes
[3, 5, 10, 15, 20] plus the fold-specific full nonzero panel were
evaluated on grouped inner validation folds, and the configured multi-criterion
rule selected a size without using the outer test fold. External labels were
not read by this workflow.

The final full-development panel size was
`3` and the panel is
FOXQ1, CEMIP, ETV4. A gene is called strictly stable
only when the general Elastic Net selection frequency is at least
0.65 and sign consistency is at least
0.8; other retained genes are named compact
consensus or exploratory panel genes in the signature table.
