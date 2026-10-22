# Post hoc directional-fidelity decomposition

Exact pinned scorer helpers reproduce every per-target fidelity and every arm aggregate (maximum absolute difference 0). This is an explanatory diagnostic, not a new selection criterion.

| Arm | Seed | Official fidelity | Bare sign agreement¹ | Called genes (sum) | Targets under budget | Zero-call targets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| current_current | 20261021 | 0.510809 | 0.514743 | 282,674 | 2/50 | 0/50 |
| current_current | 20261022 | 0.511310 | 0.515407 | 277,168 | 2/50 | 0/50 |
| current_preserving | 20261021 | 0.390880 | 0.528137 | 41,221 | 19/50 | 0/50 |
| current_preserving | 20261022 | 0.383271 | 0.532548 | 40,287 | 21/50 | 0/50 |
| adapter_current | 20261021 | 0.487432 | 0.511453 | 137,066 | 6/50 | 0/50 |
| adapter_current | 20261022 | 0.478558 | 0.514938 | 107,206 | 10/50 | 0/50 |
| adapter_preserving | 20261021 | 0.007701 | 0.600000 | 8 | 50/50 | 45/50 |
| adapter_preserving | 20261022 | 0.002536 | 0.500000 | 8 | 50/50 | 42/50 |

¹ Bare sign agreement is the mean k/n_pred over targets with nonempty scored calls; it is not DE-call precision. All current/current and current/preserving arms have50 nonempty-call targets, so their target cohorts match. Adapter/preserving has only5 and8 eligible targets. Gene-call sets differ between arms; these conditional percentages do not establish better all-gene or true-effect sign accuracy. The reference directional budget totals59,997 for every arm.

**The current response with the preserving emitter loses fidelity through reduced call coverage, despite slightly better sign agreement among its calls.** Across seeds, bare agreement rises from51.47–51.54% to52.81–53.25%, while call totals fall from277,168–282,674 to40,287–41,221. Under-budget targets rise from2/50 to19–21/50. The macro budget penalty grows from0.00393–0.00410 to0.13726–0.14928. Remaining sign agreement around53% is still weak; eliminating a null artifact does not create correct response directions.

**Adapter plus preserving emitter is almost entirely silent:** only8 scored calls per seed, with6 or4 agreeing signs, and45/50 or42/50 targets making no calls. All50 targets are below their directional budget. Its near-zero fidelity primarily reflects severe under-calling; a precision percentage based on8 pairs is not evidence of a useful response predictor.

The helper defines n_pred as predicted-significant genes with an adjudicable real direction, after own-target exclusion. k counts agreeing committed signs within that set; the reference gene need not be significant. N_conf counts reference-significant genes independently and supplies the fixed budget. Coverage n_pred/N_conf may exceed1 and is not reference-call recall. Official fidelity is k/max(n_pred,N_conf), using the helper’s degenerate conventions. Thus neither k nor N_conf should be called the number of correct DE overlaps.

Rules were not reimplemented. Per-target components, exact settings and code/input hashes are in JSON (local artifact omitted: `redesign_direction_components.json`) and CSV (local artifact omitted: `redesign_direction_components.csv`). No model, significance threshold, scorer or reserve was changed.

Raw significance counts are obtained separately with the official `de_nsig_counts` and `exclude_on_target` helpers. Here all own-excluded significant calls are direction-adjudicable, so their totals happen to equal n_pred; the definitions remain distinct. Including own genes, raw predicted call totals are retained in JSON/CSV.
