# Future research roadmap

1. Obtain approved access to JGAD000785 or equivalent acquisition-matched
   data. Prespecify a cancer-free versus cancer-associated colonic-mucosa
   contrast within common sites and stage strata. Retain patient groups,
   document tissue handling, and verify overlap and expression-panel coverage.
2. Restore the pinned raw environment and confirm archive identities and
   CEL checksums. Reproduce primary QC, limma estimates, composition and
   enrichment in a fresh worktree, then record observed run provenance.
3. Evaluate an inductive preprocessing route on HG-U219. Learn a normalization
   reference and probe summary effects using training arrays. Apply them to
   held-out CEL files individually. Check invariance to the presence of
   other held-out arrays. Compare with the cohort-RMA benchmark using the
   same patient partitions. [McCall et al.](https://doi.org/10.1093/biostatistics/kxp059)
4. Prespecify new grouped null experiments and site/handling negative controls.
   Compare probe aggregation, stable panel selection and measured composition
   where histological estimates are available. Record effect estimates and
   tested gene universes for independent pathway replication.
5. Keep GSE41258 and GSE4107 as named historical sensitivities. Freeze assay
   mapping, model, threshold and primary performance metrics before a new
   matched-acquisition confirmation. Report patient-clustered discrimination,
   calibration and operating characteristics with uncertainty.
6. Extend producer-issued manifests to raw QC, statistical genomics and
   internal validation and verify every manuscript input against its declared
   run.
