# Control-only donor sampling precision

| Context | Arithmetic log-error ratio (95% interval) | Bulk log-error ratio (95% interval) | Exact arithmetic mean variance ratio | Precision gate |
|---|---:|---:|---:|---|
| H1 | 0.9964 (0.9914–1.0014) | 0.9889 (0.9813–0.9963) | 0.9790 | FAIL |
| A | 0.9988 (0.9934–1.0041) | 0.9894 (0.9808–0.9977) | 0.9820 | PASS |
| B | 0.9911 (0.9804–1.0020) | 0.9991 (0.9887–1.0096) | 0.9748 | PASS |
| C | 0.9999 (0.9941–1.0057) | 0.9988 (0.9859–1.0117) | 0.9836 | PASS |

Ratios are stratified/uniform; lower is better. There are200 predeclared draws per policy/context. Arithmetic and bulk log-errors use the complete gene axis. Exact variance applies only to the linear arithmetic mean, not log error or bulk ratios.

All-context precision gate: **FAIL**. A separate exact-policy null pass is also required before a nonzero H1 experiment.

Monte Carlo uncertainty conditional on fixed controls; not biological prediction error or a score.

Full results (local artifact omitted: `sampling_precision_results.json`); all draws (local artifact omitted: `sampling_precision_draws.csv`).
