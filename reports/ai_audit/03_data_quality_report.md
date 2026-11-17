# Cohort and data quality audit

## Cohort reconciliation

Independent parsing of current GEO text records matched repository accession,
source tissue, title, patient or donor ID, platform, eligibility and available
primary demographics. CEL filenames were matched to the archived inventory.
This check establishes correspondence to the inventory. Archive integrity is
a separate step in the reproduction report.

| Accession | Arrays audited | Eligible | Endpoint tissues | Distinct donors or patients | Complete pairs |
|---|---:|---:|---|---:|---:|
| GSE44076 | 246 | 246 | 50 healthy, 98 distant normal, 98 tumor | 148 | 98 |
| GSE41258 | 390 | 233 | 52 normal colon, 181 primary tumor | 190 | 43 |

The GSE41258 source contains 54 normal-colon and 186 primary-tumor arrays before
eligibility filtering. Of 157 excluded arrays, 150 are outside the endpoint.
There are 77 nonaffirmative author inclusion flags and eight technical replicate
title indicators across all tissues, of which two occur in endpoint tissues.
Reasons overlap. Twelve records have no explicit tissue characteristic. All
eligible patient-tissue strata have one array. The center is described in the
study design as Memorial Sloan Kettering Cancer Center, 1992–2004. The normal
samples are from that neoplasia cohort. [GSE41258 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41258), [Sheffer et al.](https://doi.org/10.1073/pnas.0902232106)

## F01 | Acquisition and identifiability

The `adjacent_normal` model label refers to GEO's **normal distant colon mucosa**.
The authors dissected proximal surgical margins at least 10 cm from tumor.
Cases were retrospectively selected from Bellvitge surgical patients with
stage II MSS disease. Controls were recruited prospectively during colonoscopy
and had no colorectal lesions. GSE44076 has complete sample age, sex and left/right
site information. Sample-level ischemia duration, dates of processing or scanning,
array batches and quantitative histology are unavailable in its GEO metadata.
Stage and MSS are also study-level restrictions. Array metadata contain stage for
tumors, while sample-level MSI fields are empty. [GSE44076 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44076), [Sanz-Pamplona et al.](https://doi.org/10.1186/1476-4598-13-46), [COLONOMICS data descriptor](https://doi.org/10.1038/s41597-022-01697-5)

Case status and sampling procedure induce identical indicator columns for
Task B. A model with both terms is rank deficient. Age/sex/site adjustment and
latent-factor modeling cannot resolve this alias. Published colon ischemia
experiments show that FOS expression can change with handling time, so the
stress signature warrants a preanalytical interpretation as well as a
biological one. [Spruessel et al.](https://doi.org/10.2144/04366RR04)

## F08 | Technical metadata correction

`src/data/metadata.py` now retains the explicit GSE41258 scan-date characteristic.
All 233 eligible records carry a source value. Calendar strings and spreadsheet
serials coexist. No calendar conversion or batch assignment was inferred.
The clinical supplement was independently joined by explicit patient ID.
The remaining differences were source `NA` markers serialized as blanks in
21 MSI and 48 p53 fields across all records. Age, sex, location and stage are missing for nine of 52
normal-colon arrays and present for all 181 eligible tumors. Array-level MSI is
missing for 11 normals and 14 tumors. Missingness is recorded by endpoint tissue
in `results/ai_audit/geo_covariate_missingness.csv`.

## Quality control

The tracked primary route uses HG-U219 with oligo and HG-U133A with affy.
`R/utils.R` implements intensity diagnostics, RLE, PCA, clustering and
multimetric flags. `R/03_preprocess_gse41258.R` additionally runs affyPLM NUSE.
The existing audit reports no primary exclusions under the rule requiring
at least two severe independent failures. `GSM1077736` is a documented sensitivity
exclusion. PCA separation cannot identify a technical outlier by itself when
tissue and acquisition differ. The log's technical eligibility decision does not
override the clinical endpoint filters. See `reports/ai_audit/02_baseline_reproduction.md`
for execution status and `results/ai_audit/geo_metadata_source_hashes.csv` for
the exact GEO metadata inputs.
