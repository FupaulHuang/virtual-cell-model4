# P0: frozen control-only sampling precision experiment

This implements P0 in RESPONSE_RECOVERY_PLAN_20261021.md. The experiment compares
uniform sampling without replacement against stable depth stratification. No
perturbed-cell outcome, source-response label or reserve is used by this stage.

## Inputs and algorithms

Use the existing 4,000 H1 inference controls and the independently frozen
4,000-cell donor populations from the A/B/C null study. The input rows, genes,
file hashes and code are frozen before sampling outcomes. For each population,
stable-sort raw library depths, retaining frozen input row order at ties. Divide
into 400 consecutive blocks of ten and independently choose one donor per block.
Uniform sampling chooses 400 distinct cells from the same population. Both have
cell inclusion probability 1/10; arithmetic profile means and raw count means are
unbiased linear estimators. Finite-sample bulk ratios are not claimed unbiased.

Use exactly 200 seed IDs, 20261025 through 20261224 inclusive. The integer seed
actually passed to the RNG is SHA256('sampling-precision-v1:'+context+':'+seed_id)
first four bytes little-endian. Both policies use that seed. Retain every draw;
no favorable seed selection or alternative stratification is allowed.

## Precision outcomes and tests

Compare each 400-cell draw with its own complete 4,000-cell population. Normalize
on the complete gene axis. Primary errors are mean squared difference of
log1p(arithmetic CPM) and log1p(50,000*bulk probability). Also retain arithmetic/bulk
L1 error, library-depth quantiles/mean, detection error, CPM variance trace and
guide/batch composition differences. These are sampling errors, not biological
prediction errors or competition scores. Per-gene support is never selected by
observed performance.

Compute exact finite-population per-gene variance of the arithmetic mean under
both designs. Enumerated toy-population tests independently verify unbiasedness,
equal inclusion probabilities and these variance formulas. This separates
theoretical sampling precision from Monte Carlo chance. The exact formula is not
applied to nonlinear bulk ratios or log-transformed errors.

Use paired seed bootstraps (5,000 resamples, seed 20261025) for mean gains and
ratios; intervals describe Monte Carlo uncertainty conditional on these four
fixed control populations, not biological generalization.

## Prospective engineering gate

Proceed to a nonzero H1 sampling comparison only if BOTH precision and the
separately frozen exact-policy null study pass:

- H1 mean bulk log-error ratio stratified/uniform <=0.95, and the paired 95%
  interval for uniform-minus-stratified mean error is strictly positive.
- H1 arithmetic log-error ratio <=1.02, including the bootstrap ratio upper95%
  endpoint. For A, B and C, both error ratios and their upper95% endpoints <=1.05.
- Null checks require exact donor/count identity in all four contexts. Each
  uniform baseline must have median<=2 and p90<=5 calls or the null comparison
  is inconclusive. Stratified median<=max(2,2*uniform median) and
  p90<=max(5,2*uniform p90), with linear quantiles. Call counts use post-BH own-alias
  exclusion consistently; raw official nsig counts are also reported. Report
  guide/batch imbalance without changing the registered sampler.

These thresholds are engineering triage choices, not FDR guarantees or a relaxed
six-metric promotion rule. If any required gate fails, stop before another
H1 perturbation-scoring cycle and retain all results. The response-only benchmark
can proceed independently. Passing permits a fixed current-response comparison
against the existing uniform-emission arms at seeds 20261021 and 20261022; no
response amplitudes, profiles, emitter parameters or thresholds may be changed.
Any model promotion still requires the historical six-metric no-regression gate.

The user's conditional authorization for one sufficiently supported official
submission remains unused. No control-only diagnostic can satisfy that condition.
