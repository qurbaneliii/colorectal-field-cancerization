# Method selection and platform compatibility

fRMA freezes probe effects and variance estimates from reference arrays. Published
HG-U133A vectors exist. The searched Bioconductor annotation catalog did not list
an HG-U219 frozen-vector package. This is not proof that valid custom vectors
cannot be built. Cross-platform vectors must not be substituted by chip name
similarity. [fRMA method](https://pubmed.ncbi.nlm.nih.gov/20097884/)
[Reference vectors](https://bioconductor.org/packages/release/data/annotation/html/hgu133afrmavecs.html)

SCAN models probe-sequence effects within each array, and UPC extends the idea
to activity probabilities. SCAN is a credible future comparator but requires
validated platform probe sequences and dependencies beyond the preserved lock.
Neither SCAN nor UPC was installed or evaluated in this experiment. A package's
existence does not demonstrate correct coverage of this project's panels.
[SCAN method](https://pubmed.ncbi.nlm.nih.gov/22959562/)
[Implementation](https://www.bioconductor.org/packages/release/bioc/html/SCAN.UPC.html)

The implemented minimal comparator uses the locked oligo/affy/preprocessCore
environment and actual manufacturer PM indices on both platforms. Each array
undergoes independent background correction. A mean quantile grid is learned
only from training arrays. Each probeset's median within-array PM percentile is
mapped through that grid and log2-transformed; uniquely annotated probesets are
aggregated by gene median. It omits learned probe-affinity effects and is not
fRMA, standard RMA or a substitute primary analysis. This choice tests a concrete
low-memory, single-specimen interface without assuming it improves performance.

The comparison with joint RMA changes both normalization and summarization.
Differences cannot be attributed exclusively to transduction. A subsequent
training-only probe-effect model and established-method comparison are needed.
