# Constrained dispersion development results

All arms use the frozen 50-target, 400-cell H1 Training development cohort. Values below are raw official-definition metrics; no heterogeneous overall score or cross-panel anchors are used.

| Metric | Incumbent | SD1.5 | SD2.0 | Better |
| --- | ---: | ---: | ---: | --- |
| PDS | 0.910204082 | 0.910204082 | 0.910204082 | higher |
| Expression MSE | 0.970348453 | 0.970348453 | 0.970348453 | lower |
| LFC NMAE | 0.962901741 | 0.962897978 | 0.962899227 | lower |
| Direction fidelity | 0.510809091 | 0.500658067 | 0.340200986 | higher |
| Direction reach | 0.115049906 | 0.129836861 | 0.152301126 | higher |
| Significance Jaccard | 0.092176112 | 0.070284180 | 0.030675343 | higher |

## sd15: FAIL STRICT

Independent matrix invariants: True; PDS invariant: True; finite cohorts identical: True; all protected point estimates nonworse: True. Both primary gains positive: False; both paired 95% intervals above zero: False.

Evidence status: strict criterion not met. No model promotion.

| Metric | Gain in better direction | Paired 95% interval |
| --- | ---: | --- |
| PDS | +0.000000000 | [+0.000000000, +0.000000000] |
| Expression MSE | -0.000000000 | Global ratio-of-sums with whole-panel correction budget; no invalid per-target-normalized bootstrap is performed. |
| LFC NMAE | +0.000003763 | [-0.000001551, +0.000009704] |
| Direction fidelity | -0.010151024 | [-0.022504281, +0.000759886] |
| Direction reach | +0.014786956 | [+0.003368527, +0.027598704] |
| Significance Jaccard | -0.021891932 | [-0.031890120, -0.012516834] |

Achieved aggregate SD multiplier, min/median/max: [1.4850490163972543, 1.4863152287709132, 1.4883497135421642]. Target reached in 50/50 groups; duplicate-cell increase flags: 0. Maximum realized LFC drift: 0.000144273533.

With verified identical whole-panel pseudobulks and fixed reference, MSE movement comes from the capped sampling-noise correction (apart from numerical effects), not improved predicted expression means or better mean-response biology.

## sd20: FAIL STRICT

Independent matrix invariants: True; PDS invariant: True; finite cohorts identical: True; all protected point estimates nonworse: True. Both primary gains positive: False; both paired 95% intervals above zero: False.

Evidence status: strict criterion not met. No model promotion.

| Metric | Gain in better direction | Paired 95% interval |
| --- | ---: | --- |
| PDS | +0.000000000 | [+0.000000000, +0.000000000] |
| Expression MSE | -0.000000000 | Global ratio-of-sums with whole-panel correction budget; no invalid per-target-normalized bootstrap is performed. |
| LFC NMAE | +0.000002514 | [-0.000004225, +0.000008827] |
| Direction fidelity | -0.170608105 | [-0.226366473, -0.116076718] |
| Direction reach | +0.037251220 | [+0.009107568, +0.068343106] |
| Significance Jaccard | -0.061500769 | [-0.089241960, -0.034811204] |

Achieved aggregate SD multiplier, min/median/max: [1.980148351532032, 1.9834639011348387, 1.9892299762175365]. Target reached in 50/50 groups; duplicate-cell increase flags: 0. Maximum realized LFC drift: 0.000144275717.

With verified identical whole-panel pseudobulks and fixed reference, MSE movement comes from the capped sampling-noise correction (apart from numerical effects), not improved predicted expression means or better mean-response biology.

The JSON includes source/config identity checks, actual finite target IDs, paired win/loss/tie counts, strict p_adj<0.05 reference-relative TP/FP/FN, full independent count audits and terminal statuses.

- Development evidence in one previously used H1 context; not an untouched biological validation.
- The 25-target reserve remains unextracted; closed H1 Test/K562 confirmations are not reused.
- Requested aggregate residual SD is not each gene's multiplier; inspect achieved distributions and variance decreases.
- Strict numerical point nonregression is not statistical proof of absence of harm.
- Bootstrap intervals are paired target percentiles for the fixed finite cohorts, without familywise multiplicity adjustment; separate confirmation is required.
- Four thousand measured controls differ from the official context-control count; no leaderboard score is inferred.
