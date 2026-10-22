# Frozen H1 Test confirmation

The frozen selected model improves 6/6 raw metrics versus the matched neutral comparator. Both scorer jobs completed successfully using CPU/scanpy and the same reserved Test reference. The model was frozen before these outcomes, and no retuning was performed.

| Metric | Neutral | Selected | Gain in better direction | Direction |
| --- | ---: | ---: | ---: | --- |
| PDS | 0.500000 | 0.841414 | +0.341414 | higher |
| Expression MSE | 1.024846 | 0.973061 | +0.051784 | lower |
| LFC NMAE | 0.999794 | 0.959092 | +0.040701 | lower |
| Direction fidelity | 0.206063 | 0.469806 | +0.263743 | higher |
| Direction reach | 0.049854 | 0.185317 | +0.135463 | higher |
| Significance Jaccard | 0.008694 | 0.069123 | +0.060429 | higher |

The reference contains 19,443 treated cells across 100 targets, 2,000 controls and 18,080 genes. Both predictions contain 200 cells per target. The selected model has source coverage for 100/100 targets; all H1 response sources are excluded.

The JSON report records metric-specific target cohorts, paired target gains, DE call counts, output hashes, frozen-model provenance and terminal process statuses.

- This is one public H1 context and is not an estimate of hidden A/B/C leaderboard performance.
- H1 Validation was used for model selection; this Test confirms reserved target outcomes within H1, not a new cell context never used for selection.
- Both candidates were evaluated once after the selected configuration was frozen; no Test-driven parameter changes are made.
- Only raw metric values are compared on this Test reference; no Validation anchors or cross-panel scaled comparisons are used.
- Reference group sizes vary; predicted groups contain 200 cells each. NMAE may omit targets below its real-side significance gate.
- Usable cohorts come from metric_aggregation.csv n_used: direction reach uses 97 targets and NMAE uses 72; the wide aggregate count column includes the three undefined reach slots.
- Local confirmation emits 200 cells per target, while the challenge requires 400; DE detection power therefore differs from production.
- Coverage100/100 means some source evidence exists for every target, not every target-gene pair is measured in every source.
- Pooling and moment fitting alter variance; the official scorer retains its safeguards against suppressed sampling scatter.
- The selected model calls 206,159 significant target-gene pairs versus 115,411 in the reference; improved metrics do not eliminate excess DE calls.
