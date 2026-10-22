# Matched source-only response benchmark, 2026-10-21

This is a source-context **development** benchmark. K562, HCT116 and HEK293T were already inspected during adapter-v1 selection; these are not untouched validation contexts. No H1 responses, H1 evaluation, reserved responses, CD4 inputs or official scoring are used.

## Inputs and folds

Use exactly the three raw banks in the frozen adapter-v1 provenance: `data/K562_GWPS_CPM_full_statistics.npz`, `data/processed/HCT116_full_statistics.npz` and `data/processed/HEK293T_full_statistics.npz`. Require their SHA-256 values and the immutable adapter module/model hashes to match. Keep each raw/DE representation within one family; no DE representation is an extra predictor.

Each outer fold holds out one whole context. That context is absent from every training input and label, including source-common feature calculations. Train on the two remaining destination tasks, each using the other remaining context as its input. Evaluate the held-out destination from both other contexts. This one-input-training/two-input-testing distribution shift is a limitation. Destination matched controls are allowed; destination responses enter labels only.

Before outcomes, freeze source hashes, context/target/gene axes, and every sampled pair in `reports/matched_response_manifest.json` and `reports/matched_response_pairs.npz`. Evaluate every destination target with at least 20 cells. Include uncovered target/gene coordinates rather than restricting evaluation to transferred pairs. Sample 100,000 target-balanced, uniformly selected native-axis gene pairs per context without replacement within a target; exclude each target's own gene. Use seed 20261023 plus the fixed context index. Training uses 30,000 source-covered pairs per task and seed 20261021, with the frozen v1 metadata sampler. The training and evaluation pair identities are saved, not resampled after outcomes.

## Matched predictors and units

Compare zero, **matched fusion**, and adapter-v1 refitted within each outer fold. Do not use the all-context fitted adapter on a held-out context. Fix ridge regularization at 0.03 for both heads, as already selected in v1; there is no new regularization or amplitude search.

All predictors use the same raw input prior mass 100,000, minimum 20 cells, availability masks, native panels, matched controls and explicit response spaces. Labels are observed unshrunk source responses:

- Arithmetic response: `log2(target arithmetic CPM + 1) - log2(matched control arithmetic CPM + 1)`.
- Bulk response: `ln(1 + 50000 * target native-panel count probability) - ln(1 + 50000 * matched control probability)`.

The shared representation is the frozen v1 basis: pooled target-specific response plus the full-source common component. Here “specific” does not mean a separate learned coefficient for each source family. Own-target coordinates are omitted when computing each source common component. Matched fusion fixes the specific coefficient to 1, the shared coefficient to 0, and all interactions to 0. The adapter learns coefficients on that same basis and its fixed interactions. Zero fixes all coefficients to 0. Use equal valid-family pooling, no CD4 and amplitude 1. Clip every predictor identically to [-3, 3]; observed evaluation truth and training labels remain unshrunk and **unclipped**, preserving error on extreme measured responses. Matched fusion is therefore **not** the published 2/1/1 fusion with 0.6/0.3 amplitudes; it is the controlled comparator for this benchmark. The comparison tests the learned coefficient policy, including learned shared effects, conditional on a common representation.

Primary call and magnitude metrics operate in the explicit native response spaces with the common prediction clip and unclipped observed truth. Also run all predictors through the same frozen exact `+1` arithmetic inverse and natural-log bulk inverse on the **complete measured destination gene axis**, not the sampled pairs. Unmeasured union-axis placeholders are excluded and must carry no bulk mass. This separate diagnostic compares normalized arithmetic-mean compositions and bulk probabilities with their observed full-axis counterparts, excluding the own gene after normalization. Renormalized arithmetic-mean composition is not asserted to reconstruct per-cell CPM on another panel. The inverse diagnostic is not an official score and is not an advancement endpoint.

## Labels and endpoints

The independent source-provenance audit must finish before freezing/running outcomes. Available PIE Orion q-values describe stronger-KD filtered cell populations and cannot label the broader raw banks; K562's selected raw subset has no matched DE table. Unless that audit establishes exact matching, use the following fixed **exploratory magnitude calls**, never significance, official DE, or biological ground truth:

| Head | Reference or predicted magnitude call |
|---|---|
| Arithmetic | absolute response >= 0.25 log2(CPM+1) units |
| Bulk | absolute response >= 0.10 natural-log bulk units |

Stable coordinates are the complementary magnitude class; this does not establish biological absence of an effect. No p-values are manufactured from means. The same threshold applies to labels and all predictors; no threshold sweep is authorized.

Per context/head/predictor report measured source coverage, prediction call fraction, reference call fraction, signed precision (correct-sign responsive predictions / predicted calls), signed recovery (correct-sign responsive predictions / reference calls), sign accuracy among jointly called coordinates, stable false-call fraction, and class-unweighted, stable and responsive squared errors. Zero predicted calls make precision/sign quality undefined, not perfect. Retain sample counts and target-balanced weighting explicitly.

Each target has equal weight; genes sampled uniformly within that target have equal weight. “Unweighted error” means no response-class reweighting. Average context summaries equally. Uncertainty uses 2,000 paired target-cluster bootstrap draws with seed 20261023, resampling the same targets for all compared predictors. Recompute ratio metrics from target primitives in each draw. If a ratio denominator vanishes in any draw, report the undefined draws and withhold its interval rather than silently filtering them. These intervals are conditional/descriptive; they do not represent cell-sampling uncertainty or independent source contexts.

## Advancement and the conditional P2 trigger

The comparator is matched fusion, not zero. For **each head**, advancement requires higher equal-context signed recovery, nonworse equal-context signed precision and class-unweighted error, and no context-level signed precision or recovery regression. Precision must be defined in every compared context. Use absolute tolerance 1e-8 for nonregression and improvement greater than 1e-8. Both heads must pass. Report paired uncertainty even when the point gate fails; intervals do not rescue a failed gate.

P1 establishes under-recovery if adapter signed recovery is below matched fusion by more than 1e-8 in at least two of three contexts for either head. Only then may ONE small response-aware objective be preregistered in a separate P2 protocol before its outcomes. Any probability head must be called a magnitude-proxy responsivity probability unless credible matched significance labels become available. Freeze class weights, regularization, calibration, decision threshold and the fusion-based gate before P2 evaluation. Do not add methods in response to P1 results, change this benchmark, or revisit H1/reserves.

If P1/P2 cannot beat matched fusion under the registered source-only rule, stop the next H1 scoring cycle and report the missing evidence. No promotion or submission follows from this benchmark alone.
