# PREPROCESSING AND PLATFORM SHIFT

Joint RMA supports retrospective cohort-level inference. In prediction, held-out arrays help determine the cohort quantile target and probe summaries. The resulting estimates describe cohorts processed jointly.

For a single-array benchmark, fit the platform-specific background, quantile target and probe summarization reference using training CELs in each outer fold. Freeze the annotation and fitted transforms before processing each test array. Compare an array alone and with unrelated test arrays, then estimate paired score changes and patient-bootstrap intervals on common development holdouts. A validated HG-U219 frozen reference remains unavailable to the project. U133A reference support cannot be transferred to HG-U219 by assumption. [fRMA methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC2830579/) · [Training reference vectors](https://pmc.ncbi.nlm.nih.gov/articles/PMC3180392/)

At 0.41, Task C specificity was 0.558 raw versus 0.942 deposited and calibration intercept/slope changed from −2.066/2.620 to 1.372/1.420. Both cohorts change representation. [Transport evidence](../ai_audit/07_cross_platform_transport_report.md)
