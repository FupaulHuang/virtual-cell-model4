# Next gates: reduce donor-sampling error and recover real responses

## Decision from the structural experiment

The control-preserving emitter passes exact identity and control-only DE checks
in H1 and A/B/C. Nevertheless, all three new response/emitter recipes fail the
registered six-metric gate at both emission seeds. The published incumbent remains
the accepted model. The new emitter is an experimentally validated building block,
not an accepted replacement prediction recipe.

The exact scorer diagnostic identifies a coverage problem. With the current
response model, replacing the emitter sharply reduces scored significant direction
calls. With the small adapter and new emitter, only eight adjudicable significant
direction calls remain across all 50 targets at each seed. The scorer's fixed
reference directional budget is 59,997. Its fidelity formula, k/max(n_pred,N_conf),
requires enough correct calls as well as correct signs. A correct zero-response
generator cannot supply missing perturbation-response information.

This does not establish that one particular training choice caused the failure.
The adapter intervention also omitted CD4, used equal source pooling, included
learned shared effects and changed the arithmetic inverse. These changes were
declared in advance, but they prevent attributing the whole loss to learning alone.

There is a separate, measurable sampling problem. For the current response/new
emitter, population-fitted mean error is 0.00222193. Selecting 400 donors raises
the conditional expected-profile error to 0.00255527–0.00255637, before stochastic
count transport. Final error is 0.00257580–0.00258006. About 93–94% of this
population-to-output error increase is already present in donor selection.

## P0 — Test sampling precision using controls only

Before changing the response model, compare uniform sampling with one predeclared
depth-stratified policy: sort the 4,000 controls by library depth, divide them into
400 disjoint blocks of ten, and randomly choose one actual cell per block. Every
cell then has inclusion probability 1/10, the sample size remains 400, and no cell
is averaged or modified at zero effect. This preserves the population expectation
while testing whether depth stratification reduces finite-sample profile error.
It is a hypothesis, not a claimed improvement.

Freeze this policy and its seeds before testing. Measure arithmetic and bulk
sampling error across repeated control-only draws, plus detection, variance and
DE against disjoint matched controls. Require exact count identity and acceptable
null behavior under this exact policy. Avoid choosing a lucky donor set or seed.
Only then compare it with uniform sampling using the already frozen response
profiles. The current response model must remain fixed to isolate sampling.

P0 and the source-only response benchmark can proceed independently. Neither
justifies a submission until complete prediction metrics support a replacement.

## P1 — Establish a matched response-only benchmark

Do this before another full cell-generation/scoring cycle. Use only the verified
non-H1 source families for training and model selection, with each destination
context absent from every training input and label. Freeze targets and folds first.

Compare zero response, frozen source fusion and the learned adapter under the SAME
source inclusion, normalization, common-effect representation and response inverse.
Keep CD4 either absent from all compared predictors or added identically after each
predictor; do not let source omission masquerade as an isolated learning effect.
Report signal recovery against fusion, not merely a small MSE advantage over zero.

Audit source DE/control metadata before making response-confidence labels. A DE
table must match its raw source family, target, control and assay conditions.
Correlated raw/DE representations remain one family. Where valid DE labels cannot
be established, report that limitation rather than inventing significance from
means alone. Keep every H1 response representation and the 25 reserves closed.

Required evidence per held-out context: response sign among predicted calls,
reference response recovery, call coverage, magnitude error, and uncertainty.
Uniform gene-level squared error alone is not a sufficient advancement gate when
most genes have little response. Separate stable genes from responsive genes in
the report; retain an unweighted error so a new objective cannot hide overall harm.

## P2 — Test one response-aware objective on source folds

If P1 confirms response under-recovery, predeclare a small model with a response
probability head and conditional sign/magnitude prediction. A diagnostic estimate
of source-mean sign agreement is not a DE-response probability. Use source-only
folds to check whether the proposed probability is calibrated and whether it
recovers more real responses without losing sign quality.

Predeclare class balancing, regularization choices, calibration procedure and
selection criteria before reading new fold outcomes. Keep model capacity small
given the few available contexts. Require improvement over the matched frozen
fusion predictor, not only over zero response. If that fails, stop another H1
scoring cycle and identify the missing source evidence; a larger model or a new
amplitude sweep is not justified by the current results alone.

## P3 — Reconnect response and emission only after P1/P2 pass

Freeze one recipe and its source-fold evidence before revisiting count generation.
Use the control-preserving emitter's requested/population-fitted/sample-expected/
realized profile audit to separate response error from donor-sampling error. Any
new donor-balancing policy must independently pass a matched null test under that
exact sampling policy. Do not force gene totals or average donors just to recover
the old call volume.

Retain all six metrics, the same-seed baseline comparison and explicit call-count
diagnostics. A second held-out context and the untouched reserve become relevant
only after a complete candidate passes the predeclared development gate. The user
subsequently authorized one official submission if an adequately enhanced model
is found. The current candidates do not meet that condition; the authorization
has not been used.

## Scope and reason to stop here

P0–P3 are the next proposed experiments, not work already performed in the completed
structural study. The current study has supplied a concrete negative result: fixing
synthetic null behavior plus the first small adapter is insufficient. The next
decision should depend on stronger transferred-response evidence, rather than
repeating emission or global dispersion tuning on the same development outcomes.
