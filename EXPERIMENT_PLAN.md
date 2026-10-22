# Model 4 experiment and selection plan

Updated: 2026-10-21. This file defines what is being tested and how each result
affected the selected model. The experiments and two reserved confirmations are
complete; official score: **0.1470186853**. Local scores are evidence about
transfer, not estimates of the official score. See [FINAL_REPORT.md](reports/FINAL_REPORT.md).

## October 21 follow-up

The user authorized a separate [400-cell constrained dispersion comparison](docs/DISPERSION_EXPERIMENT_20261021.md).
It uses 50 metadata-selected H1 Training targets, 400 genuine reference cells and
400 predicted cells per target, plus the same 4,000 controls across all arms.
Twenty-five disjoint targets remain reserved without response extraction.
The submitted method generates the incumbent with all H1 response sources excluded.
Two transforms request 1.5× and 2× aggregate CPM residual standard deviation;
both complete matrices preserve exact cell/gene count totals and bounded CPM means.
All six raw metrics are complete: both arms preserve PDS/MSE and improve NMAE/reach
point estimates, but reduce fidelity and overlap. Both **FAIL STRICT**; neither is
promoted. [Paired results and uncertainty](reports/dispersion_experiment_results.md)
and [referee conclusion](reports/REFEREE_DISPERSION_20261021.md) retain the negative
result. No old reference anchors or closed confirmation outcomes were reused.

## Questions to resolve

1. Does broader measured same-target coverage improve transfer relative to the
   limited H1/10x data used earlier?
2. Does fitting arithmetic means separately from pseudobulk improve DE while
   maintaining expression accuracy?
3. How much does donor pooling change significance calls and direction fidelity?
4. Do source-quality filters improve transfer enough to justify their lost coverage?
5. Does a model selected on one held-out context retain useful performance on
   targets not used for selection?

These questions matter because four of the six challenge metrics concern DE,
while the preceding models were often selected using vector-space proxies.

## Data roles and leakage controls

| Data | Permitted role | Exclusion/check |
| --- | --- | --- |
| Official A/B/C controls and 300 target identities | Condition the final prediction; check submission axes. | No target-context perturbation outcomes are available. |
| Public H1 2025 responses | A source for the final A/B/C model; a destination reference for local holdouts. | Exclude **all H1 response banks** when evaluating transfer to H1. |
| Public H1 controls | Define the H1 destination's unperturbed state. | These remain available during H1 holdout evaluation. |
| H1 Validation, 50 targets | Select implementation/model settings. | Frozen reference and consistent sampling across candidates. |
| H1 Test, 100 targets | One confirmation after model selection. | Do not inspect its outcome scores to tune candidate settings repeatedly. |
| K562/HCT116/HEK293T/CD4 public responses | External response sources. | Verify provenance, per-pair availability and control matching. Exclude an entire source family if it becomes an evaluation destination. |

The current frozen H1 validation bundle is `evaluations/h1_validation/`: 18,080
genes, 2,000 controls and 200 reference cells per target except MAX, which has
161. Predictions use 200 cells per target. The manifest records the exact
reference hash and the original sampling source.

The historical Model3 H1 prediction used H1 training responses. It can be
reported as a historical comparator, but it is not a model under the same
leave-context-out constraint.

## Ordered experiment matrix

Keep source responses, target order, seeds and destination controls fixed within
each ablation. Change one modeling decision at a time where practical.

| ID | Experiment | Purpose | State / expected output |
| --- | --- | --- | --- |
| E0 | Recompute mean-response baseline and five replicate anchors on CPU/scanpy. | Put all local results on a common scale/backend. | Complete; all 12 candidates use the same fresh five-split reference. |
| E1 | Neutral control-moment emitter, pool 4. | Measure generation without perturbation information. | Complete; local score −0.023070. |
| E2 | Strongly filtered PIE sources, dual profiles, pool 4. | Test measured transfer with corrected units and emission. | Complete; score 0.105982, only 22/50 validation targets covered. |
| E3 | Same PIE sources, single profile, pool 4. | Compare the complete single-profile policy. | Complete; score 0.069775. |
| E4 | Same PIE sources, dual profiles, pool 1. | Test the alternative pooling policy. | Complete; score 0.047166. |
| E5 | Raw Orion ensemble, followed by CD4 and K562 additions; H1 excluded. | Test broader measured source coverage. | Complete; Orion 0.219109, +CD4 0.224969, +K562 0.242756 (selected). |
| E6 | Raw single profile, half amplitude, HCT-only, HCT weight 2, agreement filter. | Check whether simpler or more selective transfer improves all-metric performance. | Complete; none exceeds the selected model. Full 12-arm table linked below. |
| E7 | Freeze settings; evaluate reserved H1 Test once. | Check generalization beyond selection targets. | Complete; 6/6 aggregate raw metrics improve over neutral; H1 response sources excluded. |
| E8 | Frozen K562 confirmation with fixed target sampling. | Check transfer in a second context. | Complete; 6/6 aggregate raw metrics improve over neutral; K562/GWPS sources excluded. |

Full raw/scaled results and actual source configurations are in
[the twelve-candidate comparison](reports/h1_all_candidate_comparison.md).
The original queue and one-split pilot remain historical evidence; the completed
five-split result confirms the same selection. Scores above are local H1
Validation values, not hidden A/B/C performance.

The second confirmation context, K562, was prepared from the verified source
subset: 50 targets selected by a fixed seed and source cell count (at least 100),
their genuine cells and 2,000 measured controls. The entire K562/GWPS response
family was excluded. The globally selected settings were applied once, using
H1/Orion/CD4 as response sources.
This context is not used to tune settings; it tests whether the H1-selected policy
transfers to a different context. Its 8,246 unique-gene panel (8,248 original
columns, with two duplicate symbols combined) and sampled source cells must be
reported as limitations, and its scores cannot be compared directly
with the official 18,533-gene hidden contexts.

Initial hypotheses: mean-effect amplitude 0.6, bulk-effect amplitude 0.3,
source-common subtraction 1.0, prior strength 100,000 counts, pool 4. These were
public-method starting settings. They survived the tested comparisons; they are
not established optima. Half amplitude, a single profile and an agreement filter
were tested and not selected. New experiments must address a specific observed
failure; see [the next iteration plan](docs/NEXT_ITERATION_PLAN.md).

The current pool-1 versus pool-4 implementation changes the number of selected
donors and the depth distribution as well as averaging. Treat E4 as a comparison
of complete pooling policies; it cannot isolate variance as the sole cause of a
score difference. A variance-only claim would require matched donors and depths.

Likewise, the current single-profile arm sets both desired profiles to the
bulk-derived profile. Its effect transformation and amplitude therefore differ
from the arithmetic-profile head as well as the number of fitted profiles.
A win for the dual arm supports that complete configuration; it does not isolate
two-profile fitting as the sole causal explanation.

## Evaluation procedure

1. Confirm source banks are complete, numeric values are finite and measured
   masks/target identifiers are correct. Record target coverage before scoring.
2. Generate raw integer-count predictions from the exact candidate configuration.
3. Run `scripts/validate_prediction_receipt.py`, with H1 source exclusion for H1
   context holdouts. Inspect source manifests as well as filenames.
4. For the local scorer only, append the exact measured reference controls using
   `scripts/prepare_metric_input.py` (the wrapper in step 5 performs this automatically).
   The installed scorer requires matching
   perturbation-label sets before it applies `control_source=real`. Keep those
   local scoring files separate from submission artifacts, which contain no controls.
5. Run `scripts/run_h1_cpu_metrics.sh` using the same official `vcc2026` profile,
   CPU scanpy backend and reference controls for every candidate.
6. Scale results against the newly computed baseline and five-split anchors.
   Do not mix stored GPU/gpudge scores with these new CPU results.
7. Save the entire six-component result table and per-target diagnostics before
   selecting a candidate.

The one-split pilot supported the original freeze while the full anchor job
ran. The completed five-split reference subsequently confirmed the same winner
at 0.242756, ahead of 0.224969. Pilot outputs retain their provisional labels.
Neither reference scale was reused for H1 Test or K562; those confirmations
report paired raw metrics and the eligible target counts for each metric.

## What to record and how to select

Record overall scaled score and each of PDS, expression MSE, LFC NMAE,
direction fidelity, direction reach and significant-set Jaccard. Also record:

- Covered and uncovered target strata, and target-level response magnitudes.
- Significant-gene call counts, false-direction behavior and source agreement.
- Realized arithmetic-mean/pseudobulk errors after integer rounding.
- Cell-depth distribution, stored entries, library caps, generation time and hashes.
- Source family exclusions, data hashes and exact model arguments.

Select using the complete six-component score and the biological diagnostics.
A small aggregate gain with a severe fidelity failure is not sufficient evidence
for promotion. Prefer a stable result across many targets over a gain attributable
to a few large-effect targets. Report uncertainty from paired target-level
differences when it changes the decision; do not treat a tiny difference as proven.

The 22/50 coverage of the initial PIE-only H1 experiment limits what it can show.
It can test emission and provide a baseline, but it cannot settle the performance
of the later broad raw-source model.

## Completed selection record

`reports/MODEL_SELECTION.md` and `models/m4_multisource_dual/release.json` record:

- Exact selected source banks and weights, including any within-family mixing.
- Mean and bulk amplitudes, prior strength, centering panel, pooling, promoter
  prior decision and random seed.
- A table of completed experiments and all six scores.
- The reason for selection and observed weaknesses.
- Confirmation-test outcome, shown separately from validation.
- A statement that A/B/C hidden response data was not used.

The selected submission must match this record. To preserve the October 20
upload window, a clearly labeled provisional full-panel candidate may be
generated and packaged while its held-out evaluations run. This reversible
preparation does not select or promote the candidate. If it loses selection,
retain it as unpromoted and generate the chosen model instead. Submission still
requires the frozen decision, confirmation results and all artifact checks.

Both reserved confirmations completed without retuning. The full artifact was
promoted only after confirmation, submitted and published. Immutable pre-promotion
config/model snapshots retain their original status fields; `release.json` is
the authoritative current state. No local metric is substituted for the official
0.1470186853 score.

## Completed bulk-response ablation

The [registered four-arm bulk study](docs/BULK_ABLATION_PROTOCOL_20261021.md)
is complete. All full predictions pass engineering checks and the isolated
incumbent reproduces the previous artifact exactly. The three variants each
fail the no-regression objective; a smaller amplitude lowers uncorrected mean
error but loses more sampling-correction allowance and worsens scored MSE.
[Final report](reports/BULK_ABLATION_FINAL_20261021.md). No recipe is promoted
and the conditional25-target reserve remains unextracted.

## Completed zero-response emission diagnostic

The [null protocol](docs/EMISSION_NULL_PROTOCOL_20261021.md) compares ten400-cell
technical groups against4,000 independent measured controls, matched on batch
and guide margins. Only official DE-count diagnostics are computed; no overall
score is assigned to pseudo-perturbations. All five runs finished successfully.
The median calls are0 for measured controls,0 for donor replay,4,672 for pooling
without fitting and4,323.5 for neutral dual emission. The registered excess flag
is positive. [Full report](reports/EMISSION_NULL_FINAL_20261021.md).
