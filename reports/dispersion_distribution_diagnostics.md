# Frozen dispersion arms: distribution diagnostics

These are post-freeze descriptive comparisons with the measured H1 Training development cells. They are not new model-selection results or challenge scores. No settings were changed and no scorer was rerun.

All 50 targets have 400 measured and 400 predicted cells. One shared gate uses real-control arithmetic mean CPM >5 across 4,000 controls, with the perturbed gene excluded. Variances use per-cell CPM on the 18,080-gene native axis and sample ddof1. Log-variance summaries use the same positive-variance gene cohort for every arm within each target; zero-variance exceptions are reported separately.

| Arm | Macro mean absolute log2 variance error | Mean target-median log2 variance ratio | Underdispersed fraction | Overdispersed fraction | Detection MAE | Detection bias | CV² log-error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| incumbent | 1.427432 | -1.507147 | 0.9858 | 0.0142 | 0.080026 | 0.079166 | 1.417053 |
| sd15 | 0.624503 | -0.525018 | 0.8601 | 0.1399 | 0.058797 | 0.047519 | 0.591101 |
| sd20 | 0.547441 | 0.284126 | 0.3490 | 0.6510 | 0.068205 | -0.028471 | 0.511560 |

Negative log-variance ratio means less variable than real cells. Detection MAE is the mean absolute difference in positive-cell fractions; zero-detection MAE has the same value. Negative detection bias means the prediction has more zeros. Detection errors are fractions: 0.08 corresponds to eight percentage points. Each displayed statistic is first calculated per target and then averaged equally over targets.

Both amplified arms move the gene-wise variance summaries closer to the measured development cells. **The 1.5× arm is more consistent on detection:** it improves variance and detection-rate errors for all 50 targets and brings the variance-trace ratio close to one. **The 2× arm overshoots total variance:** it has the smaller average log-variance error, but 65.10% of gated genes are overdispersed and detection error worsens on eight targets. These findings describe distribution proximity; they do not select an arm or establish better challenge scores.

The common gate contains **10,777 genes per target**, and every gated gene has positive variance in the real cells and all three predictions. Thus the main log-variance comparison loses no genes through its positive-variance requirement. The fixed mean-matched subset averages **7,789 genes per target**.

## Paired changes relative to the incumbent

| Arm | Variance error change | Targets improved / worsened | Detection error change | Targets improved / worsened |
| --- | ---: | ---: | ---: | ---: |
| sd15 | -0.802929 | 50 / 0 | -0.021230 | 50 / 0 |
| sd20 | -0.879991 | 50 / 0 | -0.011821 | 42 / 8 |

Lower error is closer to the measured development distribution. These comparisons describe the already-frozen outputs; they do not alter which candidate is preferred by the six-metric protocol.

## Mean differences and diagnostic limits

| Arm | Mean-matched variance log-error | Mean-matched detection MAE | Predicted-zero / real-positive variance fraction |
| --- | ---: | ---: | ---: |
| incumbent | 1.516162 | 0.072542 | 0.000000 |
| sd15 | 0.568389 | 0.051139 | 0.000000 |
| sd20 | 0.521998 | 0.065437 | 0.000000 |

The mean-matched subset contains only genes whose incumbent/reference mean CPMs are within factor 1.1. It is identical across arms and is descriptive, not a tuning criterion. CV² partially accounts for scale, but neither it nor this subset proves a pure dispersion effect. Mean response errors, variable biological states and depth dependence remain relevant.

## Independent preservation audit

- sd15: exact raw gene totals and row depths in all 50 groups; zero per-gene CPM-bound violations; maximum used fraction of the declared bound 0.999977924. Every saved gene-audit axis and before/after/change/tolerance field agrees with independently recomputed counts within floating-point tolerance.
- sd20: exact raw gene totals and row depths in all 50 groups; zero per-gene CPM-bound violations; maximum used fraction of the declared bound 0.999996893. Every saved gene-audit axis and before/after/change/tolerance field agrees with independently recomputed counts within floating-point tolerance.

The full JSON retains per-target values, paired changes, zero-variance cohorts, worst 10 distribution-error changes, complete preservation checks and input hashes. The reference hash is the frozen count-identical categorical-metadata repair recorded in the development reference receipt. Closed H1 Test and K562 confirmation sets were not accessed.

Detailed diagnostics and input hashes (local artifact omitted: `dispersion_distribution_diagnostics.json`).
