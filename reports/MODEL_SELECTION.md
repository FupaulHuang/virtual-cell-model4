# Model4 selection — frozen before confirmation

Selected on 2026-10-20 at approximately 22:43 HKT: **m4_multisource_dual_v1**.
The sections below preserve the original pre-confirmation decision. Current
state: both confirmations passed, the artifact was submitted and its official
score is **0.1470186853**. See the completed outcome appended below.

## Evidence for selection

The complete-source candidate has the highest completed provisional H1 validation
score, **0.244315**, compared with **0.226368** for the next best candidate using
Orion+CD4. These are explicit one-split local diagnostic scores, not five-split
results or leaderboard estimates. The independent five-split calculation is
still running and will cross-check the comparison.

Adding K562 improves five of six raw metrics over the Orion+CD4 candidate:

| Metric | Orion+CD4 | Complete source ensemble | Better direction |
| --- | ---: | ---: | --- |
| Perturbation discrimination | 0.847755 | 0.873878 | Higher |
| Expression MSE | 0.971724 | 0.967012 | Lower |
| LFC NMAE | 0.958297 | 0.952608 | Lower |
| Direction fidelity | 0.432365 | 0.458353 | Higher |
| Direction reach | 0.212028 | 0.191048 | Higher |
| Significant-set Jaccard | 0.063011 | 0.073540 | Higher |

The reach tradeoff remains visible. Expression accuracy is still no better than
the mean-response baseline under reference scaling, and scaled LFC accuracy is
negative. The overall local gain comes primarily from discrimination, fidelity
and significance overlap. This is a stronger measured-transfer model, not a
claim of complete biological accuracy.

Other completed alternatives include neutral counts, filtered PIE transfer,
one-profile emission, one-donor pooling, half-amplitude transfer, HCT alone,
HCT-weight2, and a source-agreement filter. None exceeds the selected candidate's
provisional six-component score. The agreement filter raises discrimination but
loses other DE metrics; it is not promoted. The source-weight change gives a
small gain over the equal Orion blend but falls below the complete ensemble.

## Frozen method

- Validation sources: K562 weight2, HCT116 weight1, HEK293T weight1, CD4 weight0.5.
  All H1 response data is excluded for H1 evaluation.
- Production additionally uses public H1 responses at weight2, since H1 is an
  external context relative to A/B/C.
- Arithmetic-CPM effect amplitude0.6; pseudobulk effect amplitude0.3; clipping3;
  source-common subtraction1; prior100000 counts.
- Four pooled destination-control donors per predicted cell; seeded dual-moment
  integer count emission; seed20261020.
- No agreement filter and no promoter prior.
- Full source availability covers all300 challenge targets with at least20
  measured cells in one eligible source. Gene-level missingness remains masked.

The exact configuration, bank hashes and generated-artifact identity are in
`models/m4_multisource_dual/model.json`. The full prebuilt prediction matches this
candidate and has passed the independent all-row/group/hash validation and
strict production VCC packaging. Its existence did not determine selection.

## Reserved confirmation and promotion

After this freeze, evaluate once on:

1. H1 Test:100 target identities, no H1 response source, using its frozen reference.
2. K562:50 fixed-seed targets in a second context, excluding all K562/GWPS
   response sources; use H1+Orion+CD4 and the same frozen method settings.

Compare each with a fixed zero-response control-moment model using the same
gene axis, cell counts, controls and seeds. Do not tune on either confirmation
result or reuse H1 Validation reference scaling for another panel/context.
Confirmation outcomes and the five-split cross-check will be appended before
final promotion/submission. Negative or mixed confirmation evidence must be
reported rather than hidden.

## Completed five-split cross-check

The independent five-split CPU/scanpy evaluation confirms the same complete-source winner at **0.2427562224**, versus **0.2249690328** for Orion+CD4 (margin **0.0177871896**). All twelve candidates used the same reference bundle; each overall score equals the mean of its six reference-scaled components. No H1Test or K562confirmation outcome was available for selection. Full comparison: `reports/h1_all_candidate_comparison.json/.md`; immutable freeze snapshot cross-reference: `models/m4_multisource_dual/selection_crosscheck.json`.

## Completed confirmation and official outcome

The frozen method was evaluated once on H1 Test and K562, excluding each
respective response family. Both confirmations improved all six aggregate raw
metrics over the fixed neutral emitter. Finite target counts and limitations
are recorded in [H1 Test](h1_test_confirmation.md) and [K562](k562_confirmation.md).
No settings were changed in response to those results.

The same prebuilt artifact was promoted after confirmation. Official entry
`[PRIVATE_IDENTIFIER_REMOVED]` is published with score **0.1470186853**, rank **565** at
approximately 23:54 HKT on October 20. See [the final report](FINAL_REPORT.md).
The immutable model/config files retain their freeze-time status and hashes;
release.json (local artifact omitted: `../models/m4_multisource_dual/release.json`) records current state.
