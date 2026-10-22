# P0 sampling-policy null protocol

This prospective engineering check compares two ways of selecting 400 measured donors from a fixed 4,000-cell population. It does not fit perturbation effects, access hidden responses or reserved targets, rerun a performance score, or authorize submission. Root and independent referee reviewed the design before freezing or observing new DE outcomes.

## Fixed populations and comparisons

For A/B/C, reuse the exact donor-population and comparison-control artifacts from `evaluations/control_preserving_null_abc/{A,B,C}`. Each contains 4,000 cells. The two roles are disjoint and their population `ntc_id` margins match exactly. The earlier measured-null DE tables are reused only after checking their input, helper, parameter and result hashes.

For H1, use the 4,000 existing inference controls in `evaluations/dispersion_h1_training_400/context_H1.h5ad`. The comparison consists only of the 4,000 `non-targeting` rows in the frozen `evaluations/emission_null_400/reference.h5ad`. Confirm disjoint original cell IDs. Retain the original matched batch and guide population margins, including the previously documented joint-stratum reallocation. No actual perturbed-response file or reserve is read. Reuse the original measured-null DE table after checking its cache and parameter hashes.

Use the existing ten technical aliases in each context: `NULL_TECH_01` through `NULL_TECH_10` for A/B/C, and the previously frozen ten aliases for H1. H1 aliases are labels attached to measured controls, not new observations of those knockdowns.

## The only intervention

Call the actual `scripts/stratified_donor_sampling.py` API `select_donors(depths, seed, policy, cells=400)` for both policies.

- `uniform`: select 400 distinct rows without replacement from all 4,000 donors.
- `depth_stratified`: stable-sort the original donor rows by integer library depth, partition into 400 consecutive blocks of ten, and select one row uniformly from each block. Depth ties retain the original frozen row order.

There is no additional guide or batch quota in either draw. Report each draw's guide/batch counts, total-variation distance from the donor population and comparison population, and depth-tie frequencies; do not modify draws to improve those diagnostics. This differs from the earlier within-guide pool/replay experiment.

For each context/alias, the shared sampling seed is `(int.from_bytes(SHA256(('sampling-null-v1:' + context + ':' + alias).encode()).digest()[:4], 'little') + 20261025) % 2**32`. It is identical between policies. Freeze every selected donor index before computing new DE results. Aliases may share donors and share comparison controls, so the ten draws are correlated technical diagnostics, not independent biological replicates.

At zero effect, fit the existing control-preserving emitter to the full donor population's own arithmetic and raw-bulk profiles, select the corresponding fitted rows, and emit counts. Require each expected row probability and every emitted count to match the selected raw donor exactly. Also require 400 unique selected donors per group, valid integer counts, correct gene axes and disjointness from comparison cells. The new sampler and emitter hashes are pinned with the protocol, script and existing helper hashes.

## DE and descriptive diagnostics

Use the same pinned `cell_eval2.de_compute.compute_de` scanpy path as the earlier null test: integer counts, per-cell CPM normalization to 1e6, arithmetic LFC, epsilon 1e-9, reference-control mean CPM strictly greater than five, tie-corrected Wilcoxon and per-alias Benjamini–Hochberg. The gate and BH calculation are unchanged.

Report both (a) the official diagnostic count of `p_adj < 0.05` and (b) that count after removing the alias's own gene **after BH**, without recomputing BH. The latter is used for the engineering gate below. Apply the same post-BH exclusion to the cached measured-null table. A/B/C aliases match no genes; in H1, record whether each own gene was present and removed. The official `de_wilcoxon_nsig_counts` diagnostic itself does not remove the own gene, so these two counts are explicitly distinguished.

Report depth quantiles, detected genes per cell, gene detection error and sample CPM variance relative to the same comparison controls. These distribution and guide/batch summaries are descriptive, not additional pass/fail criteria. Keep raw-count variance, if reported, distinct from CPM variance.

## Preregistered engineering gate

Use the ten post-BH own-alias-excluded call counts per context. Compute the median and 90th percentile with NumPy's `method='linear'`. The approved thresholds are:

1. All count-identity, unique-donor, axis and disjointness checks must pass in all four contexts.
2. The uniform baseline must have median ≤2 and 90th percentile ≤5 in every context. If it fails, label the null comparison **inconclusive**, rather than certifying a policy against a poorly calibrated baseline.
3. In every context, the stratified median must be ≤`max(2, 2 × uniform median)` and its 90th percentile ≤`max(5, 2 × uniform 90th percentile)`.

All four contexts must pass. These are pragmatic engineering triage thresholds, not empirical FDR guarantees, confidence bounds or promotion criteria for perturbed-cell performance. Do not change thresholds, seeds, populations or sampling policy after seeing results. Failure does not authorize an extra policy or parameter sweep.

## Execution and retained evidence

Create only new artifacts under `evaluations/stratified_sampling_null` and `reports/stratified_sampling_null_*`; preserve existing artifacts. Freeze protocol/code/input hashes and selected indices before DE. Run eight DE jobs (four contexts × two policies) with at most four single-threaded jobs concurrently. Writable NumPy/Numba/Matplotlib runtime caches stay in the workspace. Retain the actual count files, validation receipts, sampled indices, guide/batch diagnostics, complete DE tables, source-cache checks and explicit context-wise gate decisions. No performance scorer is run.
