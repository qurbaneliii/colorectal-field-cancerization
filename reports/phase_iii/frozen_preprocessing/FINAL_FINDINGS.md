# Frozen preprocessing findings

The prespecified minimal comparator completed real-CEL evaluation at
`cb407ac0069200fce57e9a046ee667069b6ebf8e`. All 479 eligible arrays were verified
against compressed member hashes and both current-GEO TAR pins. Thirty archived
grouped partitions produced 60 fold/method metric rows and 2,064 predictions.
An independent sklearn implementation verified 744 metrics and recorded hashes.

There were 1,037 transformed-panel/external invariance comparisons. Gene values
matched exactly, maximum probability difference was 2.22e-16 and all decisions
matched. A separate ten-CEL check found zero raw-PM/background/rank differences.

Mean fixed-panel macro-F1 was 1.000 for both Task B methods and 0.994909 for both
internal Task C methods. These are conditional historical-panel results, not
new feature-discovery performance. They do not resolve acquisition confounding.

On 233 previously exposed GSE41258 endpoints, the frozen two-gene L2 refit at
0.41 achieved AUC 0.982894, macro-F1 0.949071, sensitivity 0.988950 and specificity
0.884615. The paired same-classifier joint-RMA comparator had AUC 0.982682,
macro-F1 0.854774, sensitivity 1.000 and specificity 0.615385. Conditional paired
donor-bootstrap 95% intervals for frozen-minus-comparator differences were
0.041886–0.149606 for macro-F1, -0.004051–0.004423 for AUC, and
0.150943–0.385975 for specificity. These estimates were not used to select a model.

The primary archived Elastic Net external result remains separate and unchanged.
The new method changes summarization as well as normalization, and the new
classifier differs from that historical model. Attribution solely to normalization
or an improved clinical claim is unsupported. External calibration intercepts
were -0.271 and -2.114 for the two experimental arms, respectively.

Internal calibration slopes/intercepts are unstable under near-perfect separation.
Small floating-point score changes after baseline alignment changed some fitted
calibration coefficients materially; convergence flags alone do not establish
identifiability. The raw estimates are preserved for audit, not interpreted as
reliable calibrated deployment parameters. Brier/log-loss and classification
results are retained without forcing numerical agreement.

| Gate | Status |
|---|---|
| Real-CEL fixed-panel comparison and metric replay | PASS |
| Reference/alone-versus-batch invariance | PASS |
| Unit and contract suite | PASS — 32 passed, 28 full-data tests deselected |
| Original failed external baseline attempt | FAIL — corrected eligibility alignment and rerun |
| Full nested discovery and stringent Task B sensitivity under frozen processing | NOT RUN |
| Established fRMA/SCAN comparison and clinical deployment | NOT RUN |
| Acquisition-matched independent confirmation | BLOCKED |

Issue #7 remains open for these additional research acceptance requirements.
