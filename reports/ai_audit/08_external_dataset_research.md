# Independent cohort discovery

Eligibility was based on cancer-free and cancer-associated histologically
normal colorectal mucosa, patient identifiers, documented site and acquisition,
assay annotation, and permission to access data. The audit searched GEO and
SRA, ArrayExpress/Expression Atlas, DDBJ/JGA and the original literature.
The CSV registry lists source records, publications, tissue definitions and
decision criteria. GEO permits public download and reuse subject to submitter
rights. The audit uses derived summaries and leaves raw files with their
repositories. [NCBI GEO availability terms](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)

## Acquisition-matched lead

**JGAS000655 / JGAD000785** is the closest match. Ikuta and colleagues used
endoscopic biopsies in both groups and analyzed 388 transcriptomes from
78 participants, including 17 controls and 55 patients with CRC. Six adenoma
participants are a distinct group. Multiple biopsies per participant require
patient clustering and a prespecified colonic site. The study-level
accession is public, while the 2.2 TB FASTQ dataset is Type I controlled
access. No patient-level data were accessed. Clinical site, distance,
gene mapping, overlap and input usability still require verification under
approved access. [Ikuta et al.](https://doi.org/10.1186/s12943-024-02159-9),
[NBDC dataset record](https://humandbs-production.ddbj.nig.ac.jp/dataset/JGAD000785)

## Independent sensitivity

**GSE4107** provides 10 cancer-free controls and 12 early-onset CRC patients
from Singapore on HG-U133 Plus 2.0. Study and specimen provenance differ from
the Barcelona series. Current GEO sample records identify distinct donors,
colonoscopy controls and surgical patient mucosa. The original study sampled
normal-appearing surgical tissue at least 8 cm from tumor. This repeats the
class–acquisition alias. The study is therefore a separate sensitivity for
cross-study cohort association, with limited precision and different age and
ethnicity. [GSE4107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4107),
[Hong et al.](https://pubmed.ncbi.nlm.nih.gov/17317818/)

The sensitivity kept the four Task B genes mapped on GPL570 and refit the
published Elastic Net specification on GSE44076 deposited expression ranks.
GPL570 coverage excluded SNORA12. The model and 0.50 cutoff were fixed before
GSE4107 outcomes were opened. Across 22 donors, ROC-AUC was 0.9917
(patient bootstrap 95% interval 0.9487–1.0000), macro-F1 was 0.3125, balanced
accuracy was 0.5000 and Brier score was 0.5422. Every donor was called healthy.
The ordering signal transfers across cohorts with a shared acquisition
contrast. The probability scale and fixed operating point failed to transfer.
The new cohort was kept separate from model development.

GSE95132 compares a small normal-crypt control group from ACF-bearing
individuals with matched mucosa from CRC cases. GSE90627 offers distances
within rectal cancer patients. FieldEffectCrc pools diverse source studies
whose tissue classes do not share a uniform acquisition protocol. GSE166427
is the COLONOMICS SuperSeries containing GSE44076 itself. A single-cell reanalysis
comparing GSE214695 and GSE132465 is useful for cell-type triangulation,
subject to donor and cross-study effects. The eligibility matrix documents
these decisions. [GSE95132](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95132),
[GSE90627](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90627),
[FieldEffectCrc](https://doi.org/10.18129/B9.bioc.FieldEffectCrc),
[COLONOMICS](https://doi.org/10.1038/s41597-022-01697-5),
[Single-cell reanalysis](https://doi.org/10.1016/j.isci.2026.115671)
