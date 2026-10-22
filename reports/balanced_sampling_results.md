# Balanced raw-donor precision and null study

| Context | Arithmetic error ratio (95% interval) | Bulk error ratio (95% interval) | Precision | Null median / p90 | Null |
|---|---:|---:|---|---:|---|
| H1 | 0.9845 (0.9796–0.9893) | 0.9587 (0.9518–0.9652) | FAIL | 0 / 0 | PASS |
| A | 0.9876 (0.9821–0.9934) | 0.9624 (0.9558–0.9687) | PASS | 0 / 0 | PASS |
| B | 0.9885 (0.9777–0.9994) | 0.9605 (0.9518–0.9690) | PASS | 0 / 0 | PASS |
| C | 0.9860 (0.9802–0.9915) | 0.9386 (0.9283–0.9487) | PASS | 0 / 0.1 | PASS |

Combined engineering gate: **FAIL**. Ratios are balanced/uniform; lower is better. All200 frozen draws are used in each context.

Precision and null engineering only; nonlinear bulk ratio and downstream metrics need separate evaluation.

Full evidence (local artifact omitted: `balanced_sampling_results.json`). No biological model selection or submission follows from this diagnostic.
