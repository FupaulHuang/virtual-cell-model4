# Multivariate balanced raw-donor experiment

This new experiment implements the first branch of CROSS_DOMAIN_VALIDATION_PLAN.
The earlier depth-only P0 failure and every published-model artifact are retained.
No perturbed H1 response or reserve is used. The comparison populations and 200
precision seed IDs are exactly those frozen in sampling_precision_freeze.json.

## Single fixed design

Use 4,000 measured control cells and select 400 distinct unweighted donors, each
with prescribed inclusion probability 0.1. Balance an intercept, standardized
library depth, and four standardized PC scores. PCs come from population-standardized
arithmetic CPM of control-mean-CPM>5 genes, using controls only. The feature module
records the solver, seed, singular values, rank tolerances, signs and gene masks.
No feature-count search or response-dependent feature choice is permitted.

The sampler uses cube-style martingale flight steps followed by randomized
fixed-size pivotal landing. Flight preserves the auxiliary constraints in exact
arithmetic; landing preserves marginal probabilities and sample size but can
relax auxiliary balance. Record numerical drift and final feature residuals.
Never replace landing with a nearest-subset search, retry a sample because its
profile looks poor, modify counts, or force gene totals.

Algorithm sources: [Deville and Tillé (2004)](https://doi.org/10.1093/biomet/91.4.893),
[Tillé's author overview](https://www150.statcan.gc.ca/n1/pub/12-001-x/2011002/article/11609-eng.pdf),
and [ordered pivotal sampling](https://arxiv.org/abs/1211.5442). This is a new
implementation with a specified pivotal landing, not a claim that it reproduces
an external package's particular landing optimization.

## Frozen precision comparison

Use the same P0 named RNG seeds, 200 draws per context, and identical complete-axis
arithmetic log1p(CPM) and log1p(50,000*bulk) errors. Reuse prior uniform/depth results
only after code/input/output hashes and a fresh uniform draw's numerical statistics
agree. Each new draw records actual raw profiles, detection, variance, depths,
guide/batch imbalance, selected identities and feature-balance residuals. Report
all draws and marginal inclusion frequencies. No seed is selected after outcomes.

Apply the unchanged P0 engineering thresholds against uniform: H1 bulk error
ratio<=0.95 with positive paired 95% mean-gain interval; H1 arithmetic ratio and
upper95% ratio endpoint<=1.02; A/B/C both ratios and upper endpoints<=1.05. Reuse
the fixed 5,000 paired-seed bootstrap definition. Compare against depth-only
descriptively. These intervals quantify Monte Carlo variation on fixed populations.

## Exact-policy null check

Use the same ten technical aliases, 4,000 donor and disjoint 4,000 comparison
controls per context as the completed P0 null experiment. Use the same named null
seeds. Reuse the exact frozen feature matrix from precision preparation. At zero
effect, fit each population's own moments and require every emitted count to equal
the selected donor. Test with the pinned official compute_de helper and its previous
reference-only CPM>5 gate, scanpy Wilcoxon and per-alias BH. Exclude any own-alias
row only after BH; also retain raw call counts. No guide quotas are added.

Reuse prior uniform and depth-only DE tables only after full hashes/parameters
agree. Require all four uniform baselines median<=2 and p90<=5; balanced median
<=max(2,2*uniform median), p90<=max(5,2*uniform p90), linear quantiles, and all count
identities. Report technical-group dependence and guide/batch imbalance. This is
an engineering guard, not a calibrated FDR guarantee.

## Decision and separate oracle branch

Both all-context precision and null gates must pass. Failure stops this candidate;
there is no alternate PC count, feature weight, landing rule or threshold rescue.
Passing permits inclusion as an experimental arm in the separately frozen
independent injection–recovery diagnostic. It does not select a biological model,
open reserved responses, or justify the user's conditional official submission.
The injection study's current/uniform arms remain independently authorized even
if this balanced candidate fails.
