# Prospective structural redesign experiment

Authorized by the user's request to proceed with the next work on 21 October
2026. Implements NEXT_MODEL_REDESIGN_20261021.md in isolated files. No official
submission or change to the published model is part of this experiment.

## Emitter hypothesis and fixed algorithm

Fit requested arithmetic probabilities and depth-weighted bulk probabilities
on the full inference control population (4,000 H1 cells), keeping each cell's
original library size. Start from each actual control's raw count profile,
never four-cell averages. The moment iteration is attributed to the pinned MIT
AtlasShift kernel. Use at most 150 iterations, L1 target 0.0002, maximum accepted
fit error 0.001 and maximum bulk-feasibility projection 0.03. Exact depth support
bounds replace the old 1% inward bounds, so the original control moments are
feasible. The support prior is alpha=min(0.01,0.01*response_L1), approaching zero
with the requested response. There is no constant smoothing floor.

For each sampled actual donor, let lambda be its fitted expected counts. Thin
genes where lambda<x with independent binomial removals; allocate the total
removed counts among lambda>x genes by a multinomial distribution proportional
to their surplus. This preserves every donor library exactly and E[Y|x]=lambda.
It also permits expression at donor-zero genes. Zero response returns each
donor's counts unchanged. Empirical group gene totals are not repaired.

Sample 400 donors without replacement per target, with target-specific seeded
sampling from all inference controls. Reuse donor rows between the two response
arms at a given seed. Frozen seed identifiers are 20261021 and 20261022; candidate
ranking may not choose the better seed. Report both or their paired summary.
Calibration is independent of emission seed. Report requested, fitted population,
sampled conditional expectation, and actual integer mean/bulk profiles separately.
Bulk of a sampled group is a ratio estimate; full-population expected bulk is not
claimed to be the exact expectation of that random ratio.

## Response adapter and comparisons

Use one raw bank each for K562, HCT116, and HEK293T. Every H1 response bank is
excluded from both inputs and labels. Context-level outer holdout excludes that
context from every training task. A small regularized adapter uses source-specific
and shared effects, source cell counts, availability, agreement and destination
control expression. Explicit heads are log2(CPM+1) difference and natural-log
log1p(50,000*bulk_probability) difference. Record native gene-panel and bulk
denominator limitations. The initial adapter omits CD4 inferred bulk; therefore
its factorial response intervention includes this source-policy change, which must
be disclosed, rather than attributed solely to learning. Ridge selection uses only
non-H1 context folds; no H1 outcome-driven parameter changes.

Four arms use the existing 50-target H1 development reference at 400 cells/target:
current response/current emitter, current response/new emitter, adapter/current
emitter, adapter/new emitter. The first is the verified existing incumbent at
seed20261021. New emissions use two predeclared seeds; a second current/current
seed is required for seed-paired comparisons. Retain all six official raw metrics,
uncorrected log-bulk mean error, DE diagnostics, source coverage and realized
profile drift. A second held-out source context provides adapter validation;
no general claim of improved cell-level transfer follows from source-bank error
alone. Full second-context count evaluation is required before such a claim.
It is a gated next stage: proceed to that additional count panel only if a recipe
passes the H1 no-regression gate at both seeds. If every recipe fails, retain the
source-context fold results and stop promotion testing without opening reserves.

## Gates and interpretation

First run meaningful identity/activation/conditional-expectation/cap tests and
control-only H1 and A/B/C null checks. A/B/C have a separately frozen matched
ntc_id split. Diagnostic aliases are not perturbations and must never be submitted.
Existing H1 donor-replay DE can be reused only after exact output-count identity
is proved with the new emitter. Null results are not predictive performance.

Candidate scoring requires counts, input hashes, reference/control identities,
code hashes and scorer settings to be verified. Preserve the historical strict
no-regression point gate (1e-8 tolerance) for all six raw metrics; additionally
require a gain greater than 1e-8 in direction fidelity or significance Jaccard,
identical finite target cohorts, and passage at both seeds. Report seed stability.
This gate is operational, not proof that tiny declines are statistically real.
No raw-metric average, no retrospective rescue threshold, no reserve access before
a complete recipe is selected and frozen. The 25 reserved targets remain closed.
If a structural candidate fails, retain the result and explain which metrics
improved or regressed; do not silently select a different parameter or seed.
