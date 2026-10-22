# Sampling-policy null audit

The actual uniform and global depth-stratified samplers each selected400of the same4000raw donors per technical group. No extra guide/batch quotas were imposed. Every count was preserved by the zero-effect emitter, and comparison controls were disjoint. Populations, selected indices, named seeds, code and thresholds were frozen before new DE outcomes.

| Context | Uniform median / p90 | Stratified median / p90 | Allowed median / p90 | Decision |
| --- | ---: | ---: | ---: | --- |
| H1 | 0 / 0 | 0 / 0 | 2 / 5 | pass |
| A | 0 / 1 | 0 / 0 | 2 / 5 | pass |
| B | 0 / 1 | 0 / 0 | 2 / 5 | pass |
| C | 0 / 0 | 0 / 0 | 2 / 5 | pass |

Overall engineering gate: **PASS**. Uniform median≤2/p90≤5 was separately required in every context; a failed uniform baseline is inconclusive. Quantiles use linear interpolation. These are pragmatic triage limits, not empirical FDR guarantees or promotion rules for perturbed-cell performance.

Counts above use q<0.05 after own-alias removal **after** the pinned helper’s reference-CPM>5 gate and per-alias BH. Raw official diagnostic counts are also retained. H1 cached measured-null results receive the same post-BH removal; A/B/C aliases match no genes.

The JSON report (local artifact omitted: `stratified_sampling_null_results.json`) contains group distribution summaries, cached measured-null results and provenance; the CSV (local artifact omitted: `stratified_sampling_null_counts.csv`) distinguishes raw and own-excluded call counts. Per-draw guide/batch counts and total-variation imbalances are in the frozen sampling-margins artifact. The earlier within-guide replay experiment is a different sampling policy. Shared populations make these ten aliases correlated technical diagnostics.
