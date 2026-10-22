# Model 4 implementation plan

Updated: 2026-10-21. This records the completed October 20 work in `virtual-cell-model4`.

Outcome: the October 20 implementation and submission are complete. The official
score is **0.1470186853**, up **0.4416689922** from Model3 on the same panel and
anchors. See [the final report](reports/FINAL_REPORT.md). Follow-up work is
specified in [NEXT_ITERATION_PLAN.md](docs/NEXT_ITERATION_PLAN.md); its execution
status is recorded below.

October 21 update: the user requested a referee assessment and authorized the
400-cell constrained dispersion experiment. The [top-10 audit](research/current_top10_referee_20261021.md)
is complete; MSE and reach explain 62.3% of the observed gap to the top-10 mean.
The [experiment protocol](docs/DISPERSION_EXPERIMENT_20261021.md) is frozen,
and both transformed matrices passed complete count/profile validation. All
three metric runs completed: both transformed arms fail the primary objective
(fidelity and overlap fall), despite preserved or favorable protected point
estimates. See [the referee result](reports/REFEREE_DISPERSION_20261021.md).
Neither arm is promoted. Other experiments remain proposed.

## Objective and intended deliverables

Use the lessons from Models 1–3 to build a model suited to the 2026 Virtual Cell
Challenge: predict gene-knockdown responses in a new cell context using its
unperturbed controls and response data measured elsewhere. Produce a reproducible,
validated prediction and submit a selected candidate on October 20. Report the official
result honestly, including any failure to improve.

Deliverables:

1. An evidence-based summary of earlier models and the current competition.
2. Reproducible source-response banks, model configuration, generation code and
   input/output provenance.
3. A comparison on a held-out cell context using actual generated counts and
   the official metric definitions.
4. A complete, independently validated `.h5ad` and packaged `.vcc`.
5. A production submission receipt and terminal scoring result, compared with
   the previous best score of **−0.2946503069**.

A higher score is a research objective, not an assumed result. The live leader
was **0.473365** when checked on October 20; reproducing a public method scoring about
0.15 does not establish that we are competitive with that leader.

## Why this model direction

| Finding | Implementation decision | Why it matters |
| --- | --- | --- |
| Earlier H1 responses directly cover only 25 of the 300 current targets. | Acquire broader, measured same-target responses from public K562, HCT116, HEK293T and CD4 screens. | The task is primarily transfer to a new context; a measured response to the same intervention is stronger evidence than a distant target neighbor. |
| Earlier code exponentiates differences of mean-log expression as raw-count fold changes. | Rebuild arithmetic CPM and raw pseudobulk statistics; generate counts from explicitly defined expression profiles. | The old quantities are mathematically different, especially for sparse genes. Correct units are required for meaningful DE predictions. |
| Old generators cannot add expression at donor-zero entries. | Use a tested count generator that can populate missing donor support and conserves integer cell depths. | Upregulation and detection changes must be representable. |
| Proxy scores disagreed with official submission scores; fidelity was especially weak. | Evaluate the rendered cell matrix with all six official metric definitions. | Better response-vector similarity alone does not establish better predictions. |
| Model3's GRN gain was extremely small and did not generalize reliably. | Start with a regularized multi-source response ensemble; require a measured gain before adding more model complexity. | This focuses effort on demonstrated coverage and count-generation failures. |
| The final challenge uses different cell contexts. | Exclude the entire destination response source during local validation; use its controls only. | This tests the central zero-shot context-transfer requirement. |

Evidence and source links: [research summary](docs/RESEARCH_SUMMARY.md),
[current public methods](research/current_methods_20261020.md),
[count-emission review](docs/EMISSION_REVIEW.md).

## Work packages and dependencies

Status is based on artifacts inspected on 2026-10-20. Running jobs are not
counted as completed results.

| ID | Work and concrete output | Status | Completion evidence |
| --- | --- | --- | --- |
| P0 | Audit earlier models, refresh old production score and current leaderboard. | Complete | Research summary, saved leaderboard JSON and verified production status. |
| P1 | Rebuild H1 arithmetic CPM/raw-count bank with batch-matched controls. | Complete | `data/H1_2025_full_statistics.npz` and provenance: 300 targets, 38,176 unique controls, 76,352 repeated rows checked and omitted. |
| P2 | Build strongly filtered official Arc HCT116/HEK293T DE banks. | Complete | `data/processed/PIE_*_DE_statistics.npz`, downloaded-file hashes, measured-pair masks and source reports. Together these cover only 148 official targets; they are not full coverage. |
| P3 | Build broader raw X-Atlas HCT116/HEK293T banks with matched controls; acquire K562/CD4 evidence. | Complete | All raw/DE banks audited. Combined eligible sources cover all 300 challenge targets; source-specific gene missingness remains explicit. |
| P4 | Implement configurable transfer, count generation and independent validation. | Complete | Numerical/axis tests plus a complete 360,000 × 18,533 matrix validated across every count and group. |
| P5 | Compare count-generation variants on held-out H1; recompute reference anchors on a common backend. | Complete | Twelve candidates scored against the same five-split CPU/scanpy reference bundle. Complete-source model leads at 0.242756 versus 0.224969. |
| P6 | Select and freeze the model using validation; run reserved confirmations. | Complete | Immutable model/config records; five-split ranking confirms selection. H1 Test and K562 each improve all six aggregate raw metrics with their destination response source excluded. See their confirmation reports. |
| P7 | Generate and verify the full 360,000-cell prediction; package it. | Complete for frozen candidate | All-row/group/hash and 64-bit-pointer audits passed; strict production prep accepted all 900 groups, counts preserved, no warnings. Package 3.97 GB. |
| P8 | Submit the selected candidate on October 20 and verify scoring. | Complete | Entry `[PRIVATE_IDENTIFIER_REMOVED]`; upload size/MD5 verified, official status published, score 0.1470186853. |
| P9 | Write final results and reproducibility instructions. | Complete | `reports/FINAL_REPORT.md` links the exact configuration, frozen hashes, validation, confirmations and official receipt. |

Main dependency chain: **P0 → P1/P2/P3 → P5 → P6 → P7 → P8 → P9**.
P4 proceeds alongside source acquisition. Source coverage and integrity are
checked before source banks enter a selected model.

## Model to implement and evaluate

The initial candidate transfers regularized source responses to each supplied
context's controls. Each source contributes only where it measured the requested
target and gene. Source-specific average responses are subtracted using a recorded
centering panel. Missing measurements never become zero-effect votes.

Two expression profiles are modeled separately:

- **Arithmetic per-cell CPM**, which supports differential expression effects.
- **Library-weighted pseudobulk**, which supports expression-distance and
  discrimination metrics.

Raw source data supports both profiles. DE-only source tables support the first;
their pseudobulk response is an explicitly labeled approximation. Do not describe
DE-derived approximate pseudobulk as a measured raw-count quantity.

The initial source weights and effect amplitudes came from the released
AtlasShift configuration. The completed comparison supports retaining the selected
settings among the alternatives tested; it does not establish a global optimum.
Count generation uses context-matched
controls, a configurable donor-pooling level, fitted expression moments, seeded
integer rounding and exact cell-depth conservation. The final matrix's realized
moments are recorded after rounding.

The strongly filtered and unfiltered versions of one dataset are correlated
evidence. Do not combine them as independent source families without an explicit
within-family mixing decision and an ablation.

## Execution rules and decision points

- Earlier workspaces remain read-only. New scripts, data, logs and outputs belong here.
- Use source-response metadata and public studies only. Challenge perturbation
  outcomes remain unavailable and are not used for fitting or local selection.
- Keep the 100-target H1 Test response panel out of tuning. Including its target
  names in external-source acquisition is allowed; its H1 outcomes are reserved
  for the final confirmation evaluation.
- Track source completion by verified files and process completion, not by a
  download's filename or age. Resume existing live jobs rather than launching duplicates.
- The October 20 submission is a validation submission. The final D/E/F contexts are
  released on October 22 and are a subsequent adaptation step.
- If an optional source remains unavailable near today's submission window,
  assess the best fully verified configuration available. Record the omitted
  source and coverage limitation explicitly; do not label partial acquisition
  as the planned full-source model.
- Keep the objective active until the requested model, summary and submission
  outcome have been verified. A passing smoke test is only an implementation check.

Detailed next actions are in [EXPERIMENT_PLAN.md](EXPERIMENT_PLAN.md) and
[SUBMISSION_PLAN.md](SUBMISSION_PLAN.md).

## Authorized bulk-response follow-up

The user authorized the [bulk-response plan](docs/BULK_RESPONSE_REFEREE_PLAN.md).
Its [separate protocol](docs/BULK_ABLATION_PROTOCOL_20261021.md) fixes four arms,
full-incumbent identity, unchanged emitter feasibility limits, lower genuine
mean error plus raw MSE, and nonregression in the five other metrics. Four full
predictions are validated; three new metric runs completed. All three variants
fail the strict no-regression criterion, so no reserved confirmation is triggered.
The incumbent is retained and the reserve remains unextracted.
[Final result](reports/BULK_ABLATION_FINAL_20261021.md). No official submission
was made in this stage.

## Reassessment and completed null diagnostic

After the unsuccessful scalar ablations, the user authorized proceeding with a
structural diagnostic. [Matched real controls and donor replay](reports/EMISSION_NULL_FINAL_20261021.md)
produced zero DE calls, while zero-response pooled/dual synthetic cells produced
thousands. All new reference cells are genuine non-targeting cells and are
disjoint from inference controls. No reserved perturbation outcomes were read.
The [next-model plan](docs/NEXT_MODEL_REDESIGN_20261021.md) defines the necessary
zero-effect identity contract and a subsequent response-adapter/factorial study.
No new prediction recipe is promoted and no new challenge score is claimed.
