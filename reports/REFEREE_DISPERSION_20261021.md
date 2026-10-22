# Referee assessment and completed 400-cell dispersion experiment

**Neither tested dispersion arm meets the requested objective.** Both preserve
the other four metric point estimates or improve them, but both reduce DE
direction fidelity and significant-set overlap. Keep the submitted Model4 as
the incumbent. No new model was promoted or submitted.

## Why the top 10 score higher

The official October 21 snapshot at 10:18 HKT shows a top-10 mean score of
**0.341648**, versus Model4's **0.147019**, on the same validation panel and
anchors. The current displayed leader is **0.471613** and the tenth entry is
**0.317528**. These are timestamped displayed entries, not every team's
historical best.

| Component | Model4 scaled | Top-10 mean scaled | Share of the overall gap |
| --- | ---: | ---: | ---: |
| Expression MSE | 0.037668 | 0.524676 | 41.7% |
| Direction reach | 0.079313 | 0.319617 | 20.6% |
| PDS | 0.645588 | 0.827769 | 15.6% |
| LFC NMAE | 0.138194 | 0.273039 | 11.5% |
| Direction fidelity | −0.018683 | 0.074129 | 7.9% |
| DE significance overlap | 0.000032 | 0.030657 | 2.6% |

**MSE and reach account for 62.3% of the observed gap.** Every top-10 MSE
component is at least 0.4544. Some high-ranking teams still have weak overlap;
the second-ranked entry even has a slightly negative scaled fidelity score.
Optimizing the two weakest-looking Model4 columns alone would therefore miss
most of the competitive deficit.

This decomposition is arithmetic, not a causal explanation of private methods.
All ten current API descriptions are empty. Related primary publications discuss
large perturbation corpora, biological priors and compact linear predictors,
but no reproducible recipe tied to these ten current submission IDs was verified.
Neither model names nor submission counts establish architecture, overfitting
or misconduct. Independent contexts, complete data provenance and controlled
ablations remain necessary for stronger claims.

Evidence: [dated referee audit](../research/current_top10_referee_20261021.md),
machine-readable comparison (local artifact omitted: `top10_referee_20261021.json`),
[official leaderboard source](https://virtualcellchallenge.org/api/leaderboard?get_final=false).

*Component comparison and numerical score-gap decomposition (figure excluded from the source-only export).*

## Experiment and engineering checks

The preregistered experiment uses 50 metadata-selected public H1 Training
targets, **400 genuine reference cells and 400 predicted cells per target**, and
the same 4,000 measured controls for every arm. All H1 response sources are
excluded. A disjoint 25-target reserve remains unextracted. This is new
development evidence in a previously studied context, not a new hidden-context
validation result.

Both transforms act on the same frozen 20,000 × 18,080 incumbent matrix. The
1.5× and 2× requests refer to aggregate CPM residual standard deviation. Across
all 50 groups, the achieved ranges were **1.4850–1.4883×** and **1.9801–1.9892×**.
Individual genes do not uniformly attain those multipliers.

Independent checks on the actual final integer matrices established:

- Every cell depth and every gene's group count total is preserved exactly.
- Every gene's arithmetic mean CPM stays within `max(1e-8 CPM, 0.01% of baseline)`.
- Gene/cell/target axes, nonnegative integer counts and complete 400-cell groups pass.
- No duplicate cells were introduced. Saved per-gene audits match recomputation.
- The frozen module passed 13 tests plus 24 independent exact-rational cases.
- All three official-definition CPU/scanpy metric jobs finished successfully.
  Reference, controls, configuration, cached real-DE content and scoring-copy
  identities were verified; no predicted p-values were edited.

The initial reference-category error and an early implementation-hash
synchronization stop are preserved in the audit records. They produced no
scored experimental result. The repaired reference was checked count-for-count
against the original selected source cells; the final algorithm was frozen
before the two clean arm runs.

Evidence: [protocol](../docs/DISPERSION_EXPERIMENT_20261021.md),
[implementation evidence](constrained_dispersion_implementation.md),
[wrapper review](dispersion_wrapper_review.md),
terminal process records (local artifact omitted: `dispersion_terminal_statuses.json`).

## Local raw metrics: both arms fail the primary objective

These are raw metrics on the fixed H1 development panel, not leaderboard scores.
Higher is better except for MSE and NMAE. No cross-panel anchors or average of
heterogeneous raw metrics is used.

| Metric | Incumbent | 1.5× dispersion | 2× dispersion |
| --- | ---: | ---: | ---: |
| PDS | 0.910204082 | 0.910204082 | 0.910204082 |
| Expression MSE ↓ | 0.970348453 | 0.970348453 | 0.970348453 |
| LFC NMAE ↓ | 0.962901741 | 0.962897978 | 0.962899227 |
| Direction fidelity | **0.510809091** | **0.500658067** | **0.340200986** |
| Direction reach | 0.115049906 | 0.129836861 | 0.152301126 |
| DE significance overlap | **0.092176112** | **0.070284180** | **0.030675343** |

PDS is invariant as intended. MSE is also unchanged; these results show no
expression-error gain from the altered sampling correction. NMAE changes are
small and favorable, and reach improves. However, neither candidate improves
either primary point estimate, so **both fail the preregistered rule**. A gain
in reach cannot substitute for the requested gains in fidelity and overlap.

The [paired results and uncertainty report](dispersion_experiment_results.md)
records finite target cohorts, confidence intervals, DE true/false positives
and negatives, and the independent final-matrix checks. The
[referee review](dispersion_referee_review.md) separates engineering correctness
from performance and biological evidence.

For 1.5×, the paired 95% interval for the Jaccard change is
**[−0.03189, −0.01252]**, supporting deterioration on this development panel.
Its fidelity interval **[−0.02250, +0.00076]** crosses zero: the point estimate
is lower, but this interval does not establish a fidelity decrease. At 2×,
both fidelity and Jaccard intervals are entirely negative. All metric-specific
finite target cohorts are identical across arms. These are paired target
bootstrap intervals conditional on this panel, without multiplicity adjustment.

The call diagnostics explain why fewer false positives were insufficient:

| Pooled target-gene pairs | Incumbent | 1.5× | 2× |
| --- | ---: | ---: | ---: |
| Predicted significant | 282,674 | 145,309 | 23,817 |
| True positives | 33,512 | 18,786 | 3,922 |
| False positives | 249,162 | 126,523 | 19,895 |
| False negatives | 26,485 | 41,211 | 56,075 |

The reference has 59,997 significant target-gene pairs. The 1.5× transform
removes about half the false positives, but also loses about 44% of the
incumbent's true positives. The stronger transform loses still more genuine
hits. These descriptive pooled counts use strict adjusted p < 0.05 and exclude
the perturbed gene; the official Jaccard is averaged per target and is not the
Jaccard calculated from pooled counts.

## What this teaches us

The baseline was underdispersed relative to measured H1 cells: 98.6% of eligible
genes had lower variance. The 1.5× arm brought the variance-trace ratio from
0.438 to 0.973 and improved variance and detection-rate errors on all 50 targets.
The 2× arm overshot total variance to 1.737 and worsened detection error on
eight targets, although its average per-gene log-variance error decreased.
[Distribution diagnostics](dispersion_distribution_diagnostics.md).

Thus a closer marginal cell distribution did not imply a better DE prediction.
Significance depends on the joint consequences of effect size, variability,
ranking and the reference response. Uniformly increasing dispersion cannot
repair wrong mean-response directions. This is a rejection of the two tested
global amplification policies for the stated objective; it does not prove that
all selective or context-specific calibration is ineffective.

The strongest next priority from the score audit is a controlled test of bulk
mean-response construction. No clear unit or log-base bug was found. A
[separate proposed ablation plan](../docs/BULK_RESPONSE_REFEREE_PLAN.md) distinguishes
common-response subtraction from bulk amplitude while retaining DE/PDS guards.
Its four arms are proposals only. It requires improvement in uncorrected mean
error as well as the scored MSE, so a sampling-correction gain alone cannot be
misrepresented as better mean prediction.
