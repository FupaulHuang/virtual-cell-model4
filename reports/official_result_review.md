# Independent official result review

Model4 **0.147018685** improves on Model3 **-0.294650307** by **+0.441668992**. Both receipts are terminal `published` results on `vcc2026-val-1` with the identical `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4` anchors. All six raw metrics and all six scaled components improved. The saved Model4 rank is **565**.

| Component | Model3 scaled | Model4 scaled | Change | Contribution to overall gain |
| --- | ---: | ---: | ---: | ---: |
| PDS | 0.018404 | 0.645588 | +0.627184 | +0.104531 |
| Expression MSE | 0.000000 | 0.037668 | +0.037668 | +0.006278 |
| LFC NMAE | -0.000361 | 0.138194 | +0.138555 | +0.023092 |
| Direction fidelity | -1.697095 | -0.018683 | +1.678412 | +0.279735 |
| Direction reach | -0.007310 | 0.079313 | +0.086623 | +0.014437 |
| Significance Jaccard | -0.081540 | 0.000032 | +0.081573 | +0.013595 |

The largest gain is the recovery of direction fidelity from -1.697095 to -0.018683; it explains 63.3% of the overall improvement. PDS explains another 23.7%. The new model is materially better than the prior entry, rather than merely a successful upload.

## Remaining weakness

**Gene-level response calibration remains the main unresolved limitation.** Direction fidelity is still below the context-mean baseline (`-0.018683`), and significance overlap is effectively at that baseline (`0.00003234`). Raw Jaccard is 0.030320. Strong PDS (`0.645588`) shows useful perturbation discrimination, but it does not establish accurate identification of the affected downstream genes.

Expression MSE also improves only slightly on the baseline (`0.037668`). Against the saved October 20 leader snapshot (iCell, 0.473365), the remaining overall gap is 0.326347; the largest component gap is expression MSE (0.754508 on the scaled component). Model4 has improved substantially but has not reached leading performance. This is a dated snapshot comparison, not a fresh live-leaderboard claim.

The reserved H1 and K562 confirmations each improved all six raw metrics over neutral controls. They support transfer, but do not turn local scores into official score estimates. On H1 Test, 206,159 predicted significant target-gene pairs versus 115,411 reference pairs expose a calibration concern. That is consistent with pursuing better effect magnitudes, uncertainty and DE-call calibration; it does not prove the hidden-panel mechanism, and the already-tested agreement gate was not the selected solution.

## Evidence and limits

- Current receipt: `submissions/m4_multisource_dual_v1_status.json`.
- Historical receipt: `../legacy_model3/submissions/m3_e13_full_prediction/published_status.json`.
- The JSON review records receipt hashes, exact raw/scaled changes and component contributions. Both published averages were checked against the mean of their six components.
- Rank 740 and rank 565 were observed on different dates; do not treat their difference as a controlled ranking comparison.
- Multiple model changes were bundled, so this comparison cannot assign a causal gain to any one implementation change.
- No per-context official metrics were supplied; no A/B/C-specific diagnosis is claimed.
- No model files, submissions or evaluation outcomes were changed by this review.
