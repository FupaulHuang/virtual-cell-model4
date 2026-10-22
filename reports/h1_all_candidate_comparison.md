The numerical leader on this frozen H1 Validation bundle is **raw_dual_full_cpu** (0.242756), ahead of raw_dual_cd4_cpu by **0.017787**. This report does not make the model-selection decision.

All 12 arms use the same CPU/scanpy definitions, measured controls, and fresh five-split competition-enrolled reference bundle. Every overall score was checked against the mean of exactly six official `from_replicate` components. These are local validation scores.

**Reference-scaled metrics (higher is better)**

| Arm | Covered | Overall | PDS | MSE | NMAE | Fidelity | Reach | Jaccard |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| raw_dual_full_cpu | 50/50 | 0.242756 | 0.852530 | 0.000000 | -0.034997 | 0.447099 | 0.025097 | 0.166808 |
| raw_dual_cd4_cpu | 50/50 | 0.224969 | 0.792107 | 0.000000 | -0.045133 | 0.409464 | 0.051113 | 0.142262 |
| raw_dual_hct2_cpu | 49/50 | 0.220519 | 0.772281 | 0.000000 | -0.051577 | 0.438208 | 0.003062 | 0.161143 |
| raw_dual_cpu | 49/50 | 0.219109 | 0.752455 | 0.000000 | -0.049404 | 0.421403 | 0.041092 | 0.149111 |
| raw_dual_agreement_cpu | 49/50 | 0.210776 | 0.811934 | 0.000000 | -0.060058 | 0.381076 | 0.015341 | 0.116361 |
| hct_raw_dual_cpu | 47/50 | 0.207272 | 0.671261 | 0.000000 | -0.075057 | 0.465835 | -0.002257 | 0.183849 |
| raw_dual_half_cpu | 49/50 | 0.156044 | 0.759063 | 0.000000 | -0.079814 | 0.209961 | -0.000463 | 0.047514 |
| raw_single_cpu | 49/50 | 0.142330 | 0.751511 | 0.000000 | -0.088963 | 0.199822 | -0.046848 | 0.038458 |
| pie_dual_cpu | 22/50 | 0.105982 | 0.379532 | 0.000000 | -0.076969 | 0.236118 | 0.014256 | 0.082957 |
| pie_single_cpu | 22/50 | 0.069775 | 0.377644 | 0.000000 | -0.100417 | 0.193415 | -0.098619 | 0.046626 |
| pie_dual_pool1_cpu | 22/50 | 0.047166 | 0.265295 | 0.000000 | -0.076909 | 0.046195 | 0.009413 | 0.039000 |
| neutral_cpu | 0/50 | -0.023070 | -0.012273 | 0.000000 | -0.120416 | 0.086212 | -0.097389 | 0.005447 |

**Raw metrics**

PDS, fidelity, reach and Jaccard are higher-is-better; raw MSE and NMAE are lower-is-better.

| Arm | PDS | MSE | NMAE | Fidelity | Reach | Jaccard |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| raw_dual_full_cpu | 0.873878 | 0.967012 | 0.952608 | 0.458353 | 0.191048 | 0.073540 |
| raw_dual_cd4_cpu | 0.847755 | 0.971724 | 0.958297 | 0.432365 | 0.212028 | 0.063011 |
| raw_dual_hct2_cpu | 0.839184 | 0.972358 | 0.961915 | 0.452213 | 0.173279 | 0.071110 |
| raw_dual_cpu | 0.830612 | 0.975089 | 0.960695 | 0.440609 | 0.203947 | 0.065949 |
| raw_dual_agreement_cpu | 0.856327 | 0.984206 | 0.966675 | 0.412763 | 0.183181 | 0.051901 |
| hct_raw_dual_cpu | 0.795510 | 0.969248 | 0.975095 | 0.471290 | 0.168990 | 0.080850 |
| raw_dual_half_cpu | 0.833469 | 0.997256 | 0.977765 | 0.294606 | 0.170437 | 0.022369 |
| raw_single_cpu | 0.830204 | 0.975215 | 0.982901 | 0.287605 | 0.133032 | 0.018485 |
| pie_dual_cpu | 0.669388 | 0.975623 | 0.976168 | 0.312668 | 0.182306 | 0.037573 |
| pie_single_cpu | 0.668571 | 0.975691 | 0.989330 | 0.283181 | 0.091285 | 0.021988 |
| pie_dual_pool1_cpu | 0.620000 | 0.975709 | 0.976134 | 0.181524 | 0.178401 | 0.018717 |
| neutral_cpu | 0.500000 | 1.020046 | 1.000556 | 0.209156 | 0.092276 | 0.004325 |

**Source and generation configuration**

| Arm | Sources and weights | Mode | CPM/bulk amplitudes | Pool | Agreement filter |
| --- | --- | --- | --- | ---: | --- |
| raw_dual_full_cpu | K562 raw ×2; HCT116 raw ×1; HEK293T raw ×1; CD4 DE ×0.5 | dual | 0.6/0.3 | 4 | False |
| raw_dual_cd4_cpu | HCT116 raw ×1; HEK293T raw ×1; CD4 DE ×0.5 | dual | 0.6/0.3 | 4 | False |
| raw_dual_hct2_cpu | HCT116 raw ×2; HEK293T raw ×1 | dual | 0.6/0.3 | 4 | False |
| raw_dual_cpu | HCT116 raw ×1; HEK293T raw ×1 | dual | 0.6/0.3 | 4 | False |
| raw_dual_agreement_cpu | HCT116 raw ×1; HEK293T raw ×1 | dual | 0.6/0.3 | 4 | True |
| hct_raw_dual_cpu | HCT116 raw ×1 | dual | 0.6/0.3 | 4 | False |
| raw_dual_half_cpu | HCT116 raw ×1; HEK293T raw ×1 | dual | 0.3/0.15 | 4 | False |
| raw_single_cpu | HCT116 raw ×1; HEK293T raw ×1 | single | 0.6/0.3 | 4 | False |
| pie_dual_cpu | HCT116 PIE ×1; HEK293T PIE ×1 | dual | 0.6/0.3 | 4 | False |
| pie_single_cpu | HCT116 PIE ×1; HEK293T PIE ×1 | single | 0.6/0.3 | 4 | False |
| pie_dual_pool1_cpu | HCT116 PIE ×1; HEK293T PIE ×1 | dual | 0.6/0.3 | 1 | False |
| neutral_cpu | Measured controls only | neutral | —/— | 4 | False |

**Interpretation limits**

- The numerical lead is a point-estimate ranking on 50 validation targets, not proof of superior generalization.
- Doubling the HCT116 weight changes the Orion-only score by +0.001410; adding CD4 changes it by +0.005860. These small margins do not establish stable ordering beyond this panel.
- Small absolute score margins are uncertain: candidates are correlated and no independent candidate-generation replicates or performance confidence intervals were computed.
- Five replicate-anchor splits calibrate the score scale; they are not five independent candidate-performance trials.
- NMAE can use fewer finite targets than other components; effective counts are recorded separately. The normalized expression metric is aggregate-derived, not a per-target mean.
- Multiple candidates were compared on this Validation set; selection can overfit. No Test or K562 confirmation outcomes were read for this report.
- These are local scores on public H1 data and locally estimated anchors, not hidden 2026 A/B/C leaderboard estimates.

Earlier one-split pilot values remain in `reports/h1_initial_pilot1_comparison.*` as explicitly provisional diagnostics. They were not used here. No old GPU anchors were used.

The JSON records all raw/scaled values, finite-target counts, source arguments and hashes, coverage, paired variant contrasts, reference hashes, and the receipt chain from each `local_scoring_input.receipt.json` to its prediction and source configuration.
