# Preregistered bulk-response development experiment, 2026-10-21

This protocol fixes the four authorized contrasts before full generation or their metric outcomes. It follows [the bulk-response proposal](BULK_RESPONSE_REFEREE_PLAN.md). It is an experiment protocol, not a model freeze or promotion. The earlier dispersion outcomes have already been observed; this Training cohort is development data.

## Fixed cohort and inputs

Use exactly the 50 targets and 400 genuine reference cells per target selected in cohort_freeze.json (local artifact omitted: `../evaluations/dispersion_h1_training_400/cohort_freeze.json`), with the same 4,000 controls and 18,080 genes. Predict 400 cells per target. The 25 metadata-selected reserve targets are not extracted or scored. Previously closed H1 Test and K562 confirmations remain closed. Exclude every H1 response source. Keep frozen source files, source weights and their order, full source target populations, destination controls, donors, depths, pool size 4, seed 20261021, count emitter, integer repair and existing model code unchanged. Do not apply the rejected dispersion transform.

The reference is the metadata-repaired reference.h5ad (local artifact omitted: `../evaluations/dispersion_h1_training_400/reference.h5ad`), SHA-256 `6c2f8f2becaa5f6a4fbddd884c62027fe604bc4d171064f6d917c76ea7adbdfd`. The prior categorical-level repair changed no counts, cells, genes or ordering; its provenance remains in reference_category_fix.json (local artifact omitted: `../evaluations/dispersion_h1_training_400/reference_category_fix.json`). Cohort/source-row freezes remain unchanged.

## Four arms and isolation

| Arm | Raw-source bulk common subtraction | Final fused bulk amplitude |
| --- | ---: | ---: |
| incumbent | 1.0 | 0.30 |
| common05 | 0.5 | 0.30 |
| bulk015 | 1.0 | 0.15 |
| bulk060 | 1.0 | 0.60 |

All four arithmetic-CPM heads use common subtraction 1.0 and amplitude 0.6. Center raw sources on their full source target populations with the original on-target exclusions. In `common05`, keep centered CD4 LFC and its destination-control bulk conversion fixed: only raw-source bulk centering changes. The two amplitude contrasts apply to the same final fused bulk effect, including CD4. Do not add combinations, fit target-specific amplitudes, change sources or replace a failed arm.

The new incumbent must reproduce the complete 50-target original incumbent's count matrix, row/gene axes and group order before reusing existing incumbent metrics (local artifact omitted: `../evaluations/dispersion_h1_training_400/incumbent_cpu`). File-byte identity is not required if serialization metadata differs; exact integer matrix and semantic-axis identity is required and recorded. All arms must retain exact per-row depths, finite nonnegative integral counts and validator limits. Gene totals may change intentionally, so PDS invariance is not assumed. Requested arithmetic profiles must match across arms; realized arithmetic CPM and DE behavior must be measured after fitting and integer emission, not assumed constant.

Use the original emitter's unchanged projection limit (L1 at most 0.03) and final fitting limit (maximum mean/bulk L1 at most 0.001), including its other feasibility checks. Record requested, projected, fitted and realized profiles and actual rounding drift. If any target encounters a genuine unchanged-emitter feasibility rejection, label the whole arm `GENERATION_INFEASIBLE`, preserve the failure JSON and nonzero terminal status, and do not score its partial output. This means infeasible for the pinned emitter/settings, not mathematically impossible under every generator. Input, code or provenance failures are errors, not scientific infeasibility. Continue the remaining prespecified valid arms without tolerance relaxation. Every scored arm must contain the same complete 50-target panel.

## Endpoints and descriptive diagnostics

Use the pinned official six raw metrics and CPU scanpy settings, with the exact reference, real controls and verified cached real DE. Do not compute a heterogeneous raw overall score or transfer validation anchors to this panel. Compare each arm with the existing incumbent only after full identity and provenance checks.

Primary mean biology error is defined from pooled raw count sums: first form `b[g] = log1p(50000 * counts[g] / sum(counts))` on the complete gene axis, then remove the perturbation's own gene. For each target average squared predicted-minus-reference differences over the remaining genes, then take the equal-weight mean across the 50 targets. This matches the official MSE numerator gene normalization, rather than assuming legacy `expr_mse` has the same mask. Report this uncorrected error for requested, projected, fitted and realized bulk profiles, where available. The realized profile is the primary biology endpoint.

Separate uncorrected mean residual changes from official sampling correction changes. Unlike the fixed-bulk dispersion experiment, this experiment changes the global predicted-spread bound as well as individual sampling corrections. Under identical reference/masks, the change in effective prediction correction equals the change in uncorrected error minus the change in capped MSE numerator; do not attribute a scored gain entirely to either mechanism without this decomposition.

For each target, let `b0` be the same control log-pseudobulk, `d = b_pred - b0`, and `e = b_real - b_pred`, using the same own-target exclusion. Report their mean dot product and cosine where norms are nonzero, before/after each contrast. Positive alignment indicates that extending the current log-profile direction would locally reduce mean error; it does not identify the optimum generator amplitude, which also renormalizes and projects. Do not fit a new amplitude from these diagnostics.

Predeclare descriptive target strata by the number of eligible raw source banks (0, 1, 2, or 3) under the frozen builder's target-availability rule. Record CD4 availability separately. Empty strata are reported, not merged in response to outcomes. These are diagnostic strata, not new target exclusions, selection objectives or source-weight adjustments. Own-target coordinates are excluded from all primary mean-error/alignment diagnostics. No response-derived source-disagreement cutoffs or post hoc subgroup selection are authorized.

## Decision and uncertainty

A candidate passes the prespecified development point gate only if raw normalized capped MSE is strictly lower **and** realized uncorrected log-bulk mean error is strictly lower than incumbent, while each of PDS, NMAE, fidelity, reach and Jaccard is nonworse within absolute numerical tolerance `1e-8` in its favorable direction. All finite target cohorts must be identical across compared arms for each endpoint; do not silently compare intersections. Do not redefine small protected-metric losses as success. Infeasible arms cannot pass.

Report paired target bootstrap 95% intervals using 2,000 draws and seed 20261021 for uncorrected per-target mean error and the five ordinary per-target metrics, where valid. Show effective finite target counts. Do not bootstrap the global capped MSE as if its per-target components were independent: its budget is shared across targets. The intervals are descriptive conditional on this cohort and its shared controls; targets can be correlated, and no multiplicity correction is claimed. Report strict `p_adj < 0.05` DE call counts/overlap when useful, without adapting the gate. Point dominance is the predeclared development gate; uncertainty must accompany it and cannot turn a failed point gate into a pass.

No arm is automatically promoted from this already-observed single-context development study. The authorized conditional next step is one reserved confirmation only if a candidate passes every development point gate and its paired 95% mean-error improvement interval is strictly positive. Among eligible candidates select the lowest raw normalized MSE, then the lowest realized uncorrected mean error, then registered arm order (`common05`, `bulk015`, `bulk060`) for exact ties. Freeze that one recipe before extracting any reserved responses. Compare incumbent and selected recipe once on the 25 preselected reserve targets with 400 real/predicted cells each and the same 4,000 controls, all H1 sources excluded, all six raw metrics, and no retuning. This authorization is local only and includes no official submission. The previously closed confirmations remain closed. If no arm qualifies, retain the incumbent, leave the reserve unextracted and stop without adding arms or recasting a tradeoff as improvement without affecting other metrics.

## Artifact layout

Predictions: `outputs/bulk_h1_400_{incumbent,common05,bulk015,bulk060}/prediction.h5ad`. Baseline metric source: `evaluations/dispersion_h1_training_400/incumbent_cpu`. New metrics: `evaluations/bulk_h1_training_400/{common05,bulk015,bulk060}_cpu`. Input/cache checks: `scripts/run_bulk_cpu_metrics.py` and dispersion_real_de_cache.json (local artifact omitted: `../reports/dispersion_real_de_cache.json`). Root records `reports/bulk_implementation_freeze.json`, complete incumbent identity evidence and `reports/bulk_terminal_statuses.json`; the new summary writes `reports/bulk_ablation_results.json` and `.md`. Independent engineering and final scientific verdicts go in `reports/bulk_referee_review.md`. Existing frozen files and old metric artifacts remain unchanged.
