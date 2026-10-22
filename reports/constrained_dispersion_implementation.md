# Constrained dispersion implementation evidence

Implemented only the new [module](../scripts/constrained_dispersion.py) and [tests](../tests/test_constrained_dispersion.py). Frozen prediction/scoring code and source datasets were not modified. Detailed API, hashes and observed results are in the JSON receipt (local artifact omitted: `constrained_dispersion_implementation.json`).

The transform makes integer reallocations between nearby-depth cell pairs. Each move preserves both cells' library depths and every gene's raw total exactly. Integer bounds enforce the remaining accumulated CPM budget before balanced rounding; a final recomputation from the resulting counts checks every gene against **max(1e−8 CPM, 0.01% of its incumbent mean)**. It injects no expression jitter. Identity is bit-exact, including the input dtype.

The multiplier denotes aggregate residual **standard deviation**, computed as the square root of summed sample CPM variances. Thus nominal 1.5 and 2.0 correspond to variance multipliers 2.25 and 4. The independent target-achievement tolerance is 1%. `constraints_passed` and `target_reached` are separate; `partial` and `unattainable` results retain their limitations. Unattainable refers to this bounded search, not a global infeasibility theorem.

## Tests and small real-data benchmark

All **13 tests passed**. They cover identity, exact margins, cumulative CPM bounds, equal/unequal depths, zero/single-count genes, one-cell and tied inputs, deterministic seeds, integer overflow protection, duplicate/density diagnostics and adverse rounding boundaries. The independent reviewer also reported 12 tests plus 24 rational-arithmetic CPM cases passing on the numerically identical earlier revision.

One real incumbent group, **ACVR1B, 400 cells × 18,080 genes**, was read without changing its dataset or computing response metrics:

| Requested SD | Achieved aggregate SD | Seconds | Maximum relative CPM drift | NNZ after | New duplicate cells |
| ---: | ---: | ---: | ---: | ---: | --- |
| 1.0 | 1.000000 | 0.327 | 0 | 3,850,946 | No |
| 1.5 | 1.486834 | 1.247 | 0.0098754% | 3,743,483 | No |
| 2.0 | 1.985519 | 1.190 | 0.0099978% | 3,462,616 | No |

Every row depth and gene total matched exactly. Both changed arms reached their aggregate target within the declared 1% tolerance in two passes. Per-gene median SD multipliers were only **1.1025** and **1.4331**: the aggregate increase is not a uniform per-gene increase. The largest absolute mean changes, 0.1905 and 0.1958 CPM, occur within the gene-dependent relative bound; maximum absolute change is not the applicable universal tolerance.

The final module hash is `945d102322e9e2688f64bbfae98b1f815622b354086d3f08f6e3805bc79000ea`. The benchmark and initial independent review used `6d7eebc5c0e4e0ed0418c204aa9e391f680020ab835fa717fb45c647fa160981`; subsequent additions only clarified bounded-search feasibility and added two diagnostic fields. No numerical transform changed. The final 13-test result includes the additional zero-budget test. The hash change was immediately disclosed to root and reviewer; execution records must pin the actual final hash.

## Integration and limits

Call `redistribute_dispersion(counts, multiplier, seed=20261021)`. The result contains `counts`, JSON-serializable `diagnostics`, and a `gene_stats` dictionary of arrays aligned to the gene axis. The latter includes before/after CPM means and variances, per-gene changes/tolerances, variance ratios and SD multipliers. Transformed counts use int64 to avoid overflow from a narrow input dtype.

Exact arithmetic CPM preservation is possible for equal-depth moves. General near-depth moves are approximate under the explicit per-gene bound, and tight or zero budgets can prevent a change. Do not relax that budget silently. Exact raw column totals protect the pseudobulk profile; the bounded mean drift alone is **not** a pass for the stricter measured NMAE requirement. Scored MSE's sampling correction and DE reach may also change. Independent full-row validation and six-metric comparisons determine those outcomes. This small benchmark establishes engineering feasibility and approximate runtime, not biological realism or an expected score gain.
