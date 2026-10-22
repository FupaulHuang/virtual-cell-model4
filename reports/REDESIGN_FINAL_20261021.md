# Structural redesign: null behavior improved, no replacement accepted

**Decision: retain the published Model4. Do not submit these candidates.** All
three new response/emitter recipes fail the predeclared no-regression gate at
both seeds. The user subsequently authorized one submission if a sufficiently
enhanced model is found; that condition is not met here, and the authorization
has not been used. The published official score remains **0.1470186853**.

The completed work implements a control-preserving count emitter, a small
source-context response adapter, and the full response-by-emitter experiment.
There are useful advances in diagnosis and one aspect of mean prediction, but
no accepted improvement to the complete prediction recipe.

## 1. The null defect is fixed and confirmed beyond H1

The new emitter retains every selected donor count at zero effect. Nonzero
effects use population moment calibration followed by coupled binomial thinning
and multinomial count redistribution. Each cell retains its original library
size; donor-zero genes can activate; empirical group gene totals are not forced.

Exact identity passes in H1 and A/B/C. In each challenge context, ten 400-cell
technical groups were compared with disjoint, guide-matched real controls.

| Context | Measured-null median DE calls | New zero-effect emitter median | Four-donor pooling median |
|---|---:|---:|---:|
| A | 0 | 0 | 6,004.5 |
| B | 0 | 0 | 5,687 |
| C | 0 | 0 | 5,191 |

Across A/B/C, measured-null groups have 0–1 calls, replay/new-emitter groups
0–2, and pooled groups 4,715–6,431. This identifies a distributional artifact
in our generator. It does not establish how top teams built their submissions.
Pooling changes detection, depth and variance together; these technical groups
share controls and are not independent biological replicates.
The null checks use their frozen matched replay draws; the perturbation experiment
uses uniform 400-donor draws, whose sampling error is audited separately below.

See the [complete null audit](control_preserving_null_results.md) and
exact identity gate (local artifact omitted: `control_preserving_identity_gate.json`).

## 2. The full prediction experiment fails the promotion gate

All eight prediction files passed integer-count, axis, source, donor and profile
checks. The current/current seed-one artifact is byte-identical to the incumbent,
so its successful score was reused. Seven new official local CPU/scanpy runs
finished successfully with identical reference cells, scorer settings and verified
reference-DE tables. The comparison uses 50 previously frozen H1 development
targets, 400 cells per target and two fixed seeds. Every H1 response representation
was excluded from adapter training inputs and labels.

The table averages the two seeds **within each metric**. Decisions use each seed
separately; no composite average across unlike metrics was calculated.

| Response / emitter | PDS ↑ | MSE ↓ | NMAE ↓ | Direction fidelity ↑ | Reach ↑ | Significance overlap ↑ |
|---|---:|---:|---:|---:|---:|---:|
| Current / current | 0.910612 | 0.970357 | 0.962576 | 0.511060 | 0.113587 | 0.092085 |
| Current / preserving | 0.776531 | 1.066609 | 0.970225 | 0.387076 | 0.120392 | 0.045010 |
| Adapter / current | 0.793673 | 1.006793 | 0.992349 | 0.482995 | 0.073505 | 0.069050 |
| Adapter / preserving | 0.613265 | 1.013689 | 0.993840 | 0.005118 | 0.063260 | 0.000016 |

These are raw local metrics, not leaderboard-scaled values. NMAE uses the same
42 eligible targets in every arm; the other ordinary metrics use all 50. The
normalized MSE is a full-panel ratio. The gate requires all six metrics nonworse
within 1e-8, an improvement greater than 1e-8 in fidelity or overlap, identical
finite cohorts and passage at both seeds. No candidate passes.

*Two-seed metric comparison (figure excluded from the source-only export).*

The [full paired analysis](redesign_experiment_results.md),
machine-readable evidence (local artifact omitted: `redesign_experiment_results.json`) and
metric table (local artifact omitted: `redesign_experiment_metrics.csv`) contain individual seeds,
descriptive paired intervals, coverage and per-gene profile-audit links.

## 3. Why the direction score falls despite better null behavior

The exact pinned direction metric is k/max(n_pred,N_conf), where k counts
agreeing committed directions among adjudicable predicted-significant calls,
n_pred counts those calls, and N_conf is the reference-significant budget.
It penalizes both insufficient call coverage and incorrect directions. These
quantities are distinct from true-positive DE overlap.

With current responses, the preserving emitter slightly improves conditional
sign agreement among its own calls, from 51.47–51.54% to 52.81–53.25%. Both arms
have 50 nonempty-call targets, but their called gene sets differ, so this is not
proof of improved all-gene sign prediction. Meanwhile, call totals fall from
277,168–282,674 to 40,287–41,221; under-budget targets increase from 2 to 19–21.
The resulting coverage penalty lowers the official fidelity metric.

DE overlap also falls: the incumbent recovers 33,083–33,512 of 59,997 reference
calls, while current/preserving recovers 6,664–7,246. Fewer unsupported calls were
accompanied by many fewer recovered reference calls. Adapter/preserving makes
only eight own-excluded significant calls per seed, recovering two or three
reference calls. Its near-zero fidelity mainly reflects severe under-calling.

The [direction decomposition](redesign_direction_components.md) uses the exact
pinned scorer helpers and reproduces all 400 per-target fidelity values. It is
a post hoc explanation, not a replacement selection criterion. Reference DE
calls are an experimental comparison set, not infallible biological truth.

## 4. Two distinct mean-profile findings

**The adapter improves uncorrected mean error with the current emitter.**
Own-target-excluded whole-panel log-bulk error falls by about 3.64% at both seeds,
from approximately 0.00222139 to 0.00214058. This is a real local mean-profile
improvement. Corrected MSE nevertheless worsens: the effective prediction-noise
correction decreases by about 0.0001458, outweighing the roughly 0.0000808
uncorrected-error improvement. Its other metrics also regress. Therefore it is
incorrect to describe every aspect of this arm as worse, or to promote it on
mean error alone.

**Uniform donor selection adds substantial finite-sample error.** For current
responses with the new emitter, the full-population fitted error is 0.00222193.
After selecting 400 donors, the conditional expected-profile error is
0.00255527–0.00255637; final integer error is 0.00257580–0.00258006. About 93–94%
of that increase from population fit to output is already present before count
transport. Mean realized arithmetic L1 drift from the requested profile rises
from roughly 0.0016 with the old emitter to 0.0197 with the new one. Matching
population moments does not guarantee precise profiles in a random 400-cell draw.

These conclusions use whole-axis normalization before own-target exclusion and
actual final integer counts. Requested profiles are exactly identical across
emitters for a fixed response model. The [distribution audit](redesign_distribution_audit.md)
also reports depth, detection and CPM variance across all arms and seeds.

The preserving emitter improves these distribution diagnostics even on the real
perturbation panel: with current responses, mean gene-detection error falls from
about 0.0792 to 0.0236, and the mean predicted/real CPM variance-trace ratio moves
from 0.438 to 1.180. Adapter/preserving has a variance ratio near 0.988 yet almost
no recovered significant responses. Matching variability is therefore insufficient
without correct target-specific response information and precise sampled means.

## 5. What this says about the response adapter

The adapter is a small regularized model, with explicit response units and
source-context holdout protection. Source-fold squared-error gains over zero
response were only about 0.1–0.9%; this was weak evidence for response recovery.
The experiment confirms that this prototype is insufficient for the scored
direction and overlap tasks. It does not establish that every learned adapter
would fail.

The adapter intervention also omits CD4, changes the arithmetic inverse and
learns shared-response contributions. It pools available raw families equally.
These are declared parts of the response change, so its failure cannot be
attributed to learning alone. See the source-fold report (local artifact omitted: `response_adapter_loso.json`)
and receipt metadata clarification (local artifact omitted: `redesign_receipt_metadata_notes.json`).

## 6. Next work and completed scope

The [next plan](../docs/RESPONSE_RECOVERY_PLAN_20261021.md) separates two problems:

1. Test one depth-stratified raw-donor sampling policy using controls only. It
   keeps each cell's inclusion probability at 1/10 while testing whether profile
   noise can be reduced without pooling cells or altering their counts at zero.
2. Benchmark transferred responses under matched source and normalization
   policies. Require sign, coverage and magnitude improvement over frozen fusion,
   rather than relying only on a small MSE advantage over zero response.

Neither proposal has been run in this completed experiment. The second held-out
count panel and the 25 reserved targets remain unopened because no candidate
passed the development gate. No new official submission was made.

Validation comprises 21 emitter tests, 11 adapter tests and six analysis guard
tests, plus independent code/metric reviews, exact H1/A/B/C identity checks,
nine successful A/B/C DE diagnostics and all eight successful scoring records.
See the implementation freeze (local artifact omitted: `redesign_implementation_freeze.json`),
analysis freeze (local artifact omitted: `redesign_analysis_freeze.json`),
terminal records (local artifact omitted: `redesign_terminal_statuses.json`) and
completion audit (local artifact omitted: `redesign_completion_audit.json`).
