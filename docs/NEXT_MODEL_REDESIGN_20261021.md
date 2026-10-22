# Next model: fix null behavior, then learn response specificity

The two scalar-ablation studies produced no accepted replacement for Model4.
The subsequent [control-only diagnostic](../reports/emission_null_results.md)
identifies a structural weakness worth addressing before another parameter grid.
The redesign below has now been implemented and scored at two seeds. The emitter
passes exact identity and null checks in H1 and A/B/C, but all three new recipes
fail the complete prediction gate. See the [completed referee report](../reports/REDESIGN_FINAL_20261021.md)
and [next sampling/response-recovery plan](RESPONSE_RECOVERY_PLAN_20261021.md).
The original plan is retained below as the rationale and scope of that study.

## What has actually advanced

The original Model4 improved the official score from −0.294650 to 0.147019.
The later dispersion and bulk studies did not improve the submitted model under
the specified no-regression rule. They showed that matching variance globally,
uniformly changing bulk amplitude, and restoring a shared response component
are insufficient solutions. This is useful negative evidence, not a new score gain.

The null test provides a more specific finding: measured controls and raw donor
replay produce zero significant calls in ten matched null groups, while pooled
controls and the zero-response dual emitter produce about 4,300–4,700 calls per
group. Pooling alone is sufficient to create a large apparent perturbation signal
in this H1 experiment. The DE test detects a distributional difference introduced
by synthetic-cell construction; this is not evidence that the test is malfunctioning.

The null-emitter call set overlaps 53.7% of the existing response probe's calls
under the same comparison controls. This is evidence of a substantial
response-independent background in model calls, not a classification of those
individual genes as false against actual perturbation outcomes.

## R1 — Require control-preserving behavior at zero effect

Implement a new emitter in a new module, leaving the published model frozen.
Its defining contract is: with zero perturbation effect, each output cell is the
selected measured control cell without expression modification. Use actual control
cells as the base states. Do not average four donors into one synthetic cell at
this identity point, or force each generated group to have exact population gene
totals. Match intended moments in expectation and audit realized sampling error.

For nonzero effects, use count changes that vanish continuously as the effect
vanishes. A candidate design is conditional thinning for decreases and a count-birth
process for increases, with a control-derived prior allowing activation at donor-zero
entries. This is a hypothesis requiring a separate freeze, not an already validated
algorithm. Avoid reinstating Model3's mean-log/raw-count mismatch or its inability to
activate donor-zero genes. Source effect units must remain explicit.

Required artifacts and gates:

- `scripts/control_preserving_emitter.py` and meaningful tests for exact identity,
  activation, integer counts, reproducibility, support and depth/density caps.
- Zero-response comparison against measured controls on H1 and the supplied A/B/C
  controls, using disjoint comparison cells. This needs no hidden perturbed outcomes.
- Compare cell-depth distributions, gene detection, gene variance and DE null calls.
  Require the identity emitter to match donor replay for the same sampled cells;
  do not merely raise variance until significance disappears.
- Report requested versus realized arithmetic and bulk profiles across emission
  seeds. A group mean is a sample estimate; exact per-group moment enforcement is
  not automatically a biological virtue.

The present null result is from one H1 panel. Repeat its essential check in the
challenge controls before claiming that the same failure magnitude occurs in A/B/C.

## R2 — Learn which responses transfer, rather than suppressing every response

Correcting null behavior is necessary but may initially lower recall or apparent
overlap. The completed dispersion study already showed that discarding false calls
can discard genuine hits too. The response model must improve gene selection,
direction and magnitude alongside the emitter.

Build a small regularized response adapter using the verified public source banks.
Predict response confidence and magnitude from measured source effects, source
uncertainty/cell count, availability, source agreement and destination control
expression. Keep correlated raw/DE representations within one source family.
Include shared-response components as explicit features, not an assumed universally
beneficial offset. A trained, validated soft-confidence model is different from the
hard agreement gate that previously failed.

Train and select with leave-one-source-context-out folds. If H1 is the evaluation
destination, exclude every H1 response representation from adapter training inputs
and labels, not just from final inference. Target identities and destination controls
remain allowed. Record native-panel normalization and CD4's inferred bulk quantities;
do not treat unlike assays as identical raw-count measurements. Start with a small
model appropriate to the few available source contexts, rather than presuming a
large neural network will solve the problem.

## R3 — Separate emission and response contributions

Freeze a four-arm experiment before scoring:

| Response model | Emitter | Purpose |
| --- | --- | --- |
| Current | Current | Existing incumbent |
| Current | Control-preserving | Isolate emission change |
| Adapter | Current | Isolate response prediction change |
| Adapter | Control-preserving | Test their combined effect |

Use 400 cells per target and real perturbation references for performance tests.
Report all six metrics, uncorrected mean-profile error, direction/DE call diagnostics,
source coverage and multiple emission seeds. Include a second held-out context
before making a general transfer claim. The null diagnostic is not a replacement
for perturbed-cell validation, and its pseudo-labels must never be submitted.

## Selection and honest interpretation

Do not claim performance progress from passing engineering or null checks alone.
Do not retroactively change the recorded failures of previous experiments. Retain
the requested no-regression promotion rule while reporting effect sizes, uncertainty
and tradeoff frontiers explicitly. The earlier 1e−8 point gate is an especially
strict operational rule, not statistical proof that every small observed decline
is real. Any practical replacement rule must be specified prospectively rather
than chosen to rescue a favored result.

If an overall local score is needed to compare tradeoffs, build fresh official
reference scaling on the same fixed panel/backend; do not average heterogeneous raw
metrics or reuse anchors from a different panel. Keep the 25 reserved targets closed
until a recipe is selected and frozen. No new official submission is authorized or
made by this plan.
