# Initial H1 held-out-context comparison

All arms exclude H1 response banks and use identical CPU/scanpy metric definitions, real controls and fresh five-split reference anchors. Scores are local to this public H1 bundle.

| Arm | Covered | Overall | PDS | MSE | NMAE | Fidelity | Reach | Jaccard |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| neutral_pool4 | 0/50 | -0.023070 | -0.012273 | 0.000000 | -0.120416 | 0.086212 | -0.097389 | 0.005447 |
| pie_dual_pool4 | 22/50 | 0.105982 | 0.379532 | 0.000000 | -0.076969 | 0.236118 | 0.014256 | 0.082957 |
| pie_single_pool4 | 22/50 | 0.069775 | 0.377644 | 0.000000 | -0.100417 | 0.193415 | -0.098619 | 0.046626 |
| pie_dual_pool1 | 22/50 | 0.047166 | 0.265295 | 0.000000 | -0.076909 | 0.046195 | 0.009413 | 0.039000 |

The table reports reference-scaled components. Full raw values, coverage-stratified raw summaries and post-rounding moment errors are in `h1_initial_comparison.json`.

- Public H1 reference is not the hidden A/B/C panel; local scores do not estimate leaderboard scores.
- PIE sources cover only 22 of 50 H1 validation targets.
- Pool1 versus pool4 changes donor selection and depths as well as variance.
- Single-profile mode replaces the requested arithmetic mean with the bulk-derived profile, so its comparison also changes effective arithmetic-mean effect amplitudes.
- Initial runs did not request DE-table export; significant-set call counts are unavailable.
- No H1 Test result was used in this comparison.
- The expression scorer flagged suppressed across-perturbation pseudobulk scatter from moment-constrained counts; its bounded correction remains active.
- The initial evaluator processes completed and wrote all metric artifacts, but their shell wrappers subsequently exited 127 because the live wrapper file was edited. Final scoring validates persisted run metadata; metrics were not restarted.
