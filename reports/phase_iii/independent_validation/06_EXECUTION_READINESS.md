# Execution readiness

The implemented `src/data/independent_eligibility.py` validates an explicit schema,
excludes unknown critical fields, preserves raw metadata, applies the site/distance
policy, checks donor consistency and chooses one source-ordered biopsy per donor.
`scripts/prepare_independent_validation.py` writes exclusions, eligible metadata
and an input-hashed manifest. Synthetic identifiers are required in synthetic
mode; authorized-research mode requires an approval record in approved compute.

Required metadata fields are exported as `REQUIRED` and include source identity,
donor/sample IDs, cohort, histology, anatomical site, acquisition, pretreatment,
distance, consent, QC, biopsy order, age, sex, batch, platform and overlap review.
Clinical aliases must be harmonized through a separately reviewed source mapping.
No source values or donor identities are inferred silently.

The integrated Phase II scorer supplies a frozen numeric-model interface and
conditional donor bootstrap. Before real use, produce a manifest covering source
FASTQ checksums, metadata and annotation hashes, pipeline/container or lock digest,
model/reference/mapping/cutoff fingerprints, exclusions, user-approved dataset
permission reference, random seeds, exact revision, environment and score lock.
Create scores before joining endpoint labels and record that artifact's hash.

Run `python -m pytest -q tests/test_independent_eligibility.py tests/test_frozen_cohort.py`
for synthetic contracts. Run the metadata CLI only on synthetic fixtures here.
Real sequencing, quantification and scoring are NOT RUN. Secure execution and
authorized data are external requirements; no background processing is promised.
