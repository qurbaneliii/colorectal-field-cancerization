# Task C threshold selection

The cross-platform transport threshold was locked before external labels were
read. Repeated patient-grouped GSE44076 OOF probabilities were generated after
the label-independent GPL13667/GPL96 gene intersection and within-sample rank
transformation. Probabilities were averaged across repeats per sample. The
prespecified rule selected the balanced-accuracy maximum, resolving ties by
proximity to 0.5 and then the lower threshold.

- Locked threshold: 0.41
- Original Task C signature: FOXQ1, CEMIP, ETV4
- Cross-platform signature: CEMIP, ETV4
- Label-independent common-gene universe: 12039 genes
- Aggregated GSE44076 OOF ROC-AUC: 0.9997
- Aggregated GSE44076 OOF balanced accuracy: 0.9949
- External labels used for threshold selection: no
