# Experiment queue — 2026-10-20

Selection observations are restricted to public H1 Validation. The H1 Test reference has been prepared using a fixed seed and the same measured controls, but has not been scored or used to choose parameters.

The initial PIE-only runs E2/E3/E4 compare dual profile, single profile and pooling policy on the same 22 covered targets of the 50-target H1 Validation panel. They can diagnose emission behavior but cannot decide the value of a broad source ensemble.

The next raw-source runs, once complete banks pass their audits, are fixed before viewing their validation scores:

| Name | Sources (H1 completely excluded) | Arithmetic amplitude | Bulk amplitude | Pool | Purpose |
| --- | --- | ---: | ---: | ---: | --- |
| raw_dual | HCT116 + HEK293T; weight 1 each | 0.6 | 0.3 | 4 | Broad-coverage base candidate |
| raw_dual_half | Same | 0.3 | 0.15 | 4 | Conservative effect-size ablation |
| raw_single | Same | tied to bulk | 0.3 | 4 | One-profile expression/DE tradeoff |

All use prior100000, source-common subtraction1, seed20261020,200 cells/target and no promoter prior during this first source/emission comparison. Add a promoter ablation only if its known source-independent prior is selected for the production model; do not silently add it after local selection.

If K562 raw subset acquisition completes in time, compare one additional ensemble with K562 weight2 and both Orion sources weight1, holding selected emission settings fixed. This is a source-coverage addition and must be labeled as a different candidate. CD4 remains optional until its file is complete and verified. Do not block a valid, broadly covered same-day candidate solely on an unfinished optional download.

A validation table must include all six metrics. Any selection change made after these initial comparisons requires a specific observed reason recorded in MODEL_SELECTION.md. After selecting, run the untouched Test reference once with the same external-source configuration and frozen parameters. Final A/B/C inference can then include the rebuilt H1 response source (weight2), since H1 is external to those destination contexts.

## Added before scoring, 21:43 HKT

`hct_raw_dual`: complete HCT116 raw bank alone, default dual/pool4 settings. This isolates whether broad unfiltered measured response transfer improves on the restricted PIE-only starting point. It is also the explicit single-source comparator for the upcoming HCT+HEK ensemble. Motivation: completed HCT bank has290/300 challenge targets at20+cells; HEK raw is still processing. The final source-combination decision will use the completed comparison, not acquisition order.

## Added before scoring, 22:04 HKT

`raw_dual_agreement`: same complete HCT+HEK inputs and default dual/pool4 settings, with source-direction confidence gating. Where at least two independent sources measured a gene-target pair, require at least two effects of magnitude >=0.02 log2 units agreeing with the fused sign, and at least75% of measured sources agreeing. Where only one source measured a pair, preserve the effect. Suppressed genes receive zero transferred effect in both heads. Motivation: the HCT-only count predictions call160453 DEpairs versus68223 realpairs, with only55.4% direction agreement on the real significant set. This is one targeted ablation; do not assume gating improves the six-metric outcome before scoring.

## Added before scoring, 22:09 HKT

`raw_dual_hct2`: HCT weight2,HEK weight1; all other default dual/pool4 settings unchanged. Motivation: HCT-only transfer has betterfidelity/Jaccard, while the equal blend has betterPDS/reach/NMAE; source-only on-target efficacy also differs (remaining expression .319 vs .523). Squared knockdown magnitudes give an approximately2:1 reliability starting ratio. This is a tested hypothesis, not a fitted or validated quality formula. Record allsix metrics and retain the pre-existing alternatives.

## Source-family additions, fixed before scoring at22:11HKT

`raw_dual_cd4`: default HCT1+HEK1 dual/pool4 plus publisher-quality-filtered CD4weight0.5. Once K562 is verified, `raw_dual_full` adds K562weight2 to exactly this configuration. These are source additions with fixed prior/default amplitudes; compare complete six-metric outcomes. H1 remains completely excluded for H1 Validation/Test. No source-agreement gate is silently transferred to a different number of source families.

## Selection boundary

The source-family additions, one agreement-gating comparison and one HCT-weight comparison above close the current selection queue. No further parameter grid is planned unless a concrete validation defect requires correction. Once full-source evidence and provisional reference scaling are available, freeze the chosen configuration, run the reserved H1 Test and K562 context confirmation without retuning on either, and generate/submit the fullofficial panel. Continue five-split scaling as an independent cross-check.

Selection safeguard, recorded before the complete-source score: the existence of a prebuilt full-panel artifact gives that candidate no preference. Select using held-out evidence, with fullfive-split scaling when available and clearly labeledpilot scaling otherwise. Small differences are uncertain, not proof of superiority. If another candidate wins, generate and package that model; retain the provisional complete-source artifact as unpromoted.
