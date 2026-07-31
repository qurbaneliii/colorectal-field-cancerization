# Metadata audit

## GSE44076

Class counts:

| class           |   n |
|:----------------|----:|
| adjacent_normal |  98 |
| tumor           |  98 |
| healthy         |  50 |

- Unique patient/donor IDs: 148
- Complete adjacent-normal/tumor pairs: 98
- Duplicate pair entries: 0
- Missing CEL members: 0
- Unmatched or excluded samples: 0
- Expected-count discrepancies: none

## GSE41258

Preserved canonical category counts:

| class            |   n |
|:-----------------|----:|
| primary_tumor    | 186 |
| normal_colon     |  54 |
| polyp            |  51 |
| liver_metastasis |  47 |
| lung_metastasis  |  20 |
| normal_liver     |  13 |
| other            |  12 |
| normal_lung      |   7 |

- Unique explicit patient IDs: 275
- Technical-replicate candidates: 8
- Main external-validation included samples: 233
- Missing CEL members: 0
- Excluded samples: 157

Technical-replicate candidates:

| geo_accession   | sample_title                 | tissue_class     | patient_id   |
|:----------------|:-----------------------------|:-----------------|:-------------|
| GSM1012280      | Liver Metastasis 00450UR1_ez | liver_metastasis | 450          |
| GSM1012284      | Liver Metastasis 00485UR1_ez | liver_metastasis | 485          |
| GSM1012288      | Liver Metastasis 00823UR1_ez | liver_metastasis | 823          |
| GSM1012291      | Liver Metastasis 00838UR1_ez | liver_metastasis | 838          |
| GSM1012449      | Lung Metastasis 14756VR1_ez  | lung_metastasis  | 14756        |
| GSM1012495      | Primary Tumor A5135AR2_ez    | primary_tumor    | A5135        |
| GSM1012653      | Normal Colon C0667HR1_ez     | normal_colon     | C0667        |
| GSM1012658      | Cell line HCT15rehyb         | other            | Cell line    |

## Parsing rules

1. GSE44076 tissue labels are mapped only from explicit `source_name_ch1`;
   patient IDs come only from the `individual id` characteristic.
2. GSE41258 preserves the explicit `tissue` characteristic before canonical
   mapping. Only Primary Tumor and Normal Colon are eligible for main external
   validation.
3. GSE41258 samples must also carry the author flag `included in analysis: Yes`.
4. Titles matching `(?i)(?:_ez|rehyb)` are
   excluded deterministically as technical-replicate candidates.
5. Missing or unresolved metadata causes exclusion, never imputation.
