# Real-data execution

Analysis commit `cb407ac0069200fce57e9a046ee667069b6ebf8e` completed on verified current-GEO CELs. All 246 GSE44076 and 233 endpoint GSE41258 arrays were independently extracted; 390 external raw-array inventory/QC evidence remains in the primary audit. No raw download was repeated. The local archived joint-RMA matrices are reused baselines, not newly computed RMA.

Source TAR SHA-256 values

- `GSE41258_RAW.tar` — `740aaf0c5ec4c0ba72107b1d4ef50b5a309cb5cab88aa29ea73fc48594f778ee`
- `GSE44076_RAW.tar` — `059abdaa76f731c9d1af116875bb60663fb875c296bbf01e353de48eeed0a030`

There are 60 fold/method rows, 2064 held-out score rows and 1037 invariance comparisons. Maximum gene absolute difference is 0; maximum probability difference is 2.22e-16. All class decisions match.

The independent raw-PM experiment checked ten arrays and found zero differences in background-corrected PM values and ranks. Input/member hashes, manufacturer mapping hashes, train/test IDs, quantile references, numeric models, environments, score locks, independent metric replay and conditional intervals are retained with execution evidence.

The first slow extraction implementation was stopped after feasibility checking; a direct-reader implementation produced identical statistics on the checked sample and completed the extraction. Old partial caches were not mixed into this run. The preserved Windows R environment emitted locale and package-build-version warnings; no dependency version was changed.
