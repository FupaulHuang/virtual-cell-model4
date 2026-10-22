# Matched source-only response benchmark

**Labels are exploratory magnitude calls, not DE significance.** No H1 responses or official scoring were used. These source contexts were previously examined during v1 selection; this is a development benchmark.

| Context | Head | Predictor | Signed precision | Signed recovery | Call fraction | Unweighted MSE |
|---|---|---|---:|---:|---:|---:|
| K562 | log2fc | zero | undefined | 0.000000 | 0.000000 | 0.084289 |
| K562 | log2fc | matched_fusion | 0.204238 | 0.047778 | 0.064989 | 0.100389 |
| K562 | log2fc | adapter_v1 | undefined | 0.000000 | 0.000000 | 0.084047 |
| K562 | bulk_delta | zero | undefined | 0.000000 | 0.000000 | 0.009376 |
| K562 | bulk_delta | matched_fusion | 0.170199 | 0.044866 | 0.064009 | 0.011690 |
| K562 | bulk_delta | adapter_v1 | undefined | 0.000000 | 0.000000 | 0.009359 |
| HCT116 | log2fc | zero | undefined | 0.000000 | 0.000000 | 0.091543 |
| HCT116 | log2fc | matched_fusion | 0.214477 | 0.152182 | 0.175848 | 0.138675 |
| HCT116 | log2fc | adapter_v1 | undefined | 0.000000 | 0.000000 | 0.090886 |
| HCT116 | bulk_delta | zero | undefined | 0.000000 | 0.000000 | 0.005723 |
| HCT116 | bulk_delta | matched_fusion | 0.155456 | 0.075155 | 0.070452 | 0.008291 |
| HCT116 | bulk_delta | adapter_v1 | undefined | 0.000000 | 0.000000 | 0.005716 |
| HEK293T | log2fc | zero | undefined | 0.000000 | 0.000000 | 0.062500 |
| HEK293T | log2fc | matched_fusion | 0.178721 | 0.184448 | 0.202624 | 0.123477 |
| HEK293T | log2fc | adapter_v1 | undefined | 0.000000 | 0.000000 | 0.062079 |
| HEK293T | bulk_delta | zero | undefined | 0.000000 | 0.000000 | 0.003920 |
| HEK293T | bulk_delta | matched_fusion | 0.116375 | 0.099705 | 0.089034 | 0.007094 |
| HEK293T | bulk_delta | adapter_v1 | undefined | 0.000000 | 0.000000 | 0.003912 |

Both-head advancement versus matched fusion: **False**. Under-recovery established: **True**.

Matched fusion fixes the common representation to specific coefficient1/shared0 with equal raw-family pooling, amplitude1 and no CD4. It is not the published production fusion. The held-out context is absent from every training input and label; the all-context fitted adapter is never used as a fold predictor.

The JSON retains per-context stable/response errors, coverage, paired target-cluster intervals, undefined precision denominators and the fusion-based gate. Full measured-axis inverse diagnostics are separate normalized mean-composition/bulk errors, not reconstruction of per-cell CPM or official scores.

[Source-label audit](matched_response_source_audit.md) rejected mismatched PIE q-values as raw-population labels. [Frozen protocol](../docs/MATCHED_RESPONSE_BENCHMARK_PROTOCOL_20261021.md) fixes thresholds, folds and criteria.

One separately preregistered response-aware P2 objective may now be designed using clearly qualified magnitude-proxy responsivity labels. No threshold sweep, new H1 scoring cycle, promotion or submission follows automatically.
