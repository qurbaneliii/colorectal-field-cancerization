# Final acceptance code review

The review followed the raw-CEL workflow from the input pins through the
published artifacts. It covered both the requested baseline revision
`9e68c4d70b93e8e4cf48471762fcf11ae98f455a` and the acceptance fixes.
Run IDs, analysis and workflow SHAs, and execution results are recorded in
the PR description and in the hosted evidence bundles.

## Findings addressed

| Finding | Effect | Resolution |
| --- | --- | --- |
| The hosted workflow uploaded logs and input pins without the generated publication files. The figure test also accepted a pre-existing image. | Nine renderings of three internal-validation figures could survive from checkout and appear to pass regeneration. | Generate the Task A confusion matrix and the Task B and C discrimination and calibration panels from the grouped out-of-fold predictions. Record their plotted data and source hashes. Require fresh publication outputs, verify them against the run's inventory, and upload them with both commit identities. |
| Task C verified the digests present in its lock without requiring all four inputs. The external outcomes were opened before the development refits. | An incomplete lock could reach the evaluation stage. The data-access order made the fit/evaluation boundary harder to audit. | Require the four digests and agreement with the unique selected GSE44076 threshold. Check the development sample IDs and refit both representations before opening the external outcomes. External platform gene identifiers are sufficient for the intersection. |
| The manuscript discussed cutoff and calibration transport without quantifying the separately verified expression-route sensitivity. Author contributions and conflicts also needed author confirmation. | Readers could miss the magnitude of the expression-route effect or mistake draft declarations for verified submission metadata. | Report the fixed-cutoff sensitivity with the same panel, probe mapping and external cohort. Identify the author details and declarations as pending, separately from the scientific manuscript gate. |
| The abstract labeled bootstrap medians as Task C point estimates. | It displayed 0.985 and 0.561 where the locked run produced point estimates of 0.983 and 0.558. | Report the actual point estimates with the same patient-cluster intervals. Label internal bootstrap medians as medians. |
| Two sample distribution figures used unseeded point jitter. The checked-in Task C cards predated the complete lock. | Repeat plots moved individual points. The cards lacked the development input fingerprints and the ordered cross-platform rank universe. | Seed horizontal offsets from the project seed, preserving the expression and composition scores on the vertical axis. Refresh the checked-in model cards from the verified current-GEO raw-CEL matrices, retaining their actual generation revision and input hashes. |

## Scientific review

- **Groups and model selection.** Outer and inner splitters check for
  donor/patient overlap. Fold-local selectors and residualizers fit in the
  enclosing training partition. The compact-panel decision is made inside
  each outer split. The final panel and hyperparameters use the entire
  development cohort before transport. The SVM comparator's supervised
  selector sees its enclosing training partition before the grouped
  calibration split. Its calibration scores are therefore not cross-fitted
  through feature selection, although the nested-CV test patients remain
  outside that partition. The primary Elastic Net does not use this
  calibration procedure.
- **Raw preprocessing.** Each cohort undergoes RMA jointly before splitting
  or evaluation. Grouped CV and external metrics are conditional on
  cohort-processed expression. Frozen-reference preprocessing for a single
  incoming array has not been evaluated.
- **Threshold and external outcomes.** The Task C cutoff uses repeated
  GSE44076 grouped predictions averaged across repeats within each sample.
  The development diagnostic is conditional on a panel and hyperparameters
  selected using all development samples. GSE41258's outcomes enter only
  evaluation. Its tumor versus normal-colon endpoint differs from Task B's
  healthy versus adjacent-normal question. GSE4107 is a separate frozen
  association sensitivity and shares the acquisition contrast. These
  public cohorts have been analyzed previously.
- **Statistics.** The primary DE comparisons use prespecified limma
  designs, including patient fixed effects for paired tumor/adjacent
  comparisons, with BH correction. Composition and stress analyses are
  sensitivity analyses. Patient-cluster resampling preserves repeated
  tissue measurements and treats the fitted scores and threshold as fixed.
  The permutation summaries hold the observed-label split assignments and
  modal hyperparameters fixed. They are conditional cohort diagnostics.
  Task B permutations do not identify a cancer-specific effect in the
  presence of the acquisition alias.
- **Artifact provenance.** The published archive pins distinguish the
  current GEO bytes from the historical manifest. Hosted output inventories
  are revision-specific. Cross-run comparisons check scientific values,
  identities and decisions. Run timestamps, source hashes, and rendering
  metadata can legitimately change. The figure evidence additionally binds
  the repaired figures to their saved prediction inputs.
  The checked-in external report carries the verified native execution's
  revision and hashes, including the two matrices and two fitted binaries
  kept outside Git. Its snapshot verifier checks the published files.
  Fresh execution verifies the matrix and model bytes as well.

## Native and hosted output comparison

The pinned Linux replay of `9e68c4d` recovered the same panels, 0.41
development cutoff, external metrics, and patient-cluster interval endpoints.
In one Task B Elastic Net outer fold, the raw-CEL log loss and Brier score
differed from the Windows snapshot by about `5.2e-9` and `1.7e-9`.
The cohort discrimination conclusions were unchanged.

The seeded ranked analysis returned 1,371 adjacent-versus-healthy GO terms
on Linux and 1,370 on Windows. The additional term was `GO:0002604`
(Linux BH-adjusted p-value 0.0472). The two runs selected the same terms
for the main enrichment table. The checked-in enrichment and internal
model comparison snapshots use the hosted values. The Windows run remains
recorded in the raw-CEL audit. This comparison treats individual
near-threshold ranked terms as sensitive to the execution environment.

[Pinned Linux run and artifacts](https://github.com/qurbaneliii/colorectal-field-cancerization/actions/runs/43195342594)

## Open research and submission requirements

| Requirement | Status | Evidence needed |
| --- | --- | --- |
| Historical TAR byte identity and mismatch attribution | BLOCKED | Recover the original TARs, or historical member fingerprints for attribution. Exact byte-for-byte execution requires the original TARs. |
| Acquisition-matched independent Task B validation | BLOCKED | Obtain approved NBDC Type I use of JGAD000785 in an eligible environment. Verify donor, tissue site and sampling distance eligibility before a locked evaluation. |
| Author and submission metadata | BLOCKED | Supply author names, affiliations, correspondence details, approved contributions, funding and conflicts declarations. |
| Clinical performance | NOT RUN | The current study estimates retrospective cohort discrimination and transport. Prospective assay validation and clinical utility are outside the available evidence. |

The two evidence-export inputs are `github-hosted` and `current-geo`.
An optional `analysis_sha` replays a prior full commit with the same
input pins and dependency locks. A replay records missing or stale outputs
explicitly. A run at its own workflow commit fails the publication evidence
gate if any required output is stale.

[Raw-CEL continuation and archive hashes](15_raw_reproduction_continuation.md) ·
[Task C transport evidence](07_cross_platform_transport_report.md) ·
[Publication-readiness gates](../publication_readiness_final.md)
