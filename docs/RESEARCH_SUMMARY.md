# Model 4: evidence and modeling decisions

Started and submitted 2026-10-20 (Asia/Hong_Kong). Status: published. Model4 scored
**0.1470186853**, versus Model3 **−0.2946503069**, on the same validation panel and
reference anchors. Rank was **565** at approximately 23:54 HKT. The full result,
reproducibility instructions and limitations are in [the final report](../reports/FINAL_REPORT.md).

## What the earlier experiments established

The published results progress from model1 pseudobulk transfer (-0.990279), through conservative/control replay (~-0.302), model2 delta transfer (-0.297379), and model3 H1/GRN hybrid (-0.2946503069). The model3 score was rechecked using the live production CLI on 2026-10-20: entry [PRIVATE_IDENTIFIER_REMOVED] is published. The dominant scaled failure was DE direction fidelity (-1.697095), with little set overlap (raw Jaccard .000393823). Raising perturbation discrimination alone did not improve the overall score.

Evidence: ../legacy_model1/outputs/submission_results.md; ../legacy_model2/submissions/*_submit.json; ../legacy_model3/docs/IMPLEMENTATION_STATUS.md; the production vcc status response.

The previous H1-only supervised data directly covers only25/300 current targets; available low-efficacy 10x K562 covers about11 usable targets. The model3 GRN exact-definition local gain is only .00001098 and does not establish useful transfer to new contexts. Old delta-space proxies disagreed with the actual leaderboard. Same-context heldout-target validation does not test the main2026 generalization requirement.

A concrete implementation mismatch was found independently in model2/model3: learned differences of mean(log1p(CP10K)) were exponentiated as raw-count fold changes. Moreover, multiplying nonzero counts cannot introduce expression in an originally zero cell/gene. These issues must be corrected, although they do not by themselves prove the cause of the low official score.

## Current challenge and competitive context

Official website https://virtualcellchallenge.org/ confirms no new training responses this year; validation has300 CRISPRi targets across3 anonymized contexts and final release is October22. Current production format:360000 cells,18533 genes,400 cells per context-target pair, raw nonnegative integers, no controls. Six reference-scaled metrics are averaged, and own-target expression is excluded. Two scoring submissions per UTC day.

A public production leaderboard snapshot retrieved at 23:13 HKT is `reports/leaderboard_20261020_late.json` (1601 entries): leader Illumina AI/iCell .4733651969; twentieth .3099502990; hundredth .2651830284. These are validation scores and do not forecast final-phase performance. Ranks and entry counts are retrieval-time snapshots.

## Model4 direction

1. Broaden measured same-target responses using public genome-scale K562, public CD4 DE summaries and official Arc PIE X-Atlas HCT116/HEK293T summaries. Keep full source-panel centering, assay provenance and missing measurements explicit.
2. Rebuild local public2025 H1 arithmetic CPM and raw-count summaries with batch-matched controls. Repeated controls across public split files must be deduplicated and verified.
3. Transfer regularized measured response effects to the destination's control means. Retain separate arithmetic-normalized and library-weighted profiles where source data supports them; identify approximations explicitly for DE-only sources.
4. Use AtlasShift's MIT count emitter, with depth-matched control pooling, fitted expression moments and integer count conservation. Record realized moments after rounding, rather than assuming requested effects were realized.
5. Evaluate rendered counts in a heldout source cell context, excluding that entire source from response fitting. Preserve heldout test targets for confirmation after validation choices.
6. Package only after independent structure/count checks and production vcc prep. Submit today's selected candidate and wait for the actual official outcome.

The completed work improved the official score by **0.4416689922**. Its strongest
scaled component is discrimination (0.645588); direction fidelity (−0.018683)
and significant-set Jaccard (0.000032) remain near or below their reference
baseline. The result supports broader measured transfer and corrected count
generation as a useful combined change; it does not isolate their causal
contributions or approach the saved 0.473365 leader score.

Twelve local candidates were compared on H1 Validation. The chosen model scored
0.242756 there, and both reserved context/target confirmations improved all six
raw aggregate metrics over a neutral comparator. Those panels have different
gene axes, eligible-target counts and reference scales; their values must stay
separate from the official result. Excess local DE calls and source-dependent
normalization remain concrete limitations. See the [next iteration plan](NEXT_ITERATION_PLAN.md)
for controlled tests addressing them.

## Public methods and attribution

AtlasShift: https://github.com/kaipengm2/Virtual-Cell-Challenge-2026.git (pinned commit in reports/upstream_sources.json), MIT. Its authors report .1545618019 on September10 using multi-source transfer and dual-moment counts. Sidechain reports .1684 on October16; rpinho/vcc2026 reports .1528 on October18. Their reported scores belong to their full configurations and data, not to this workspace.

Official source: https://huggingface.co/datasets/arcinstitute/PIE_x_atlas_orion . DE-derived arithmetic-mean response heads do not supply measured library-weighted raw pseudobulk; any approximation will be labeled.

Numerical source kernels remain in pinned external checkouts. Original model4 orchestration, receipts and tests live in scripts/, reports/ and tests/.
