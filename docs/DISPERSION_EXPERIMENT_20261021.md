# Constrained dispersion development experiment — 2026-10-21

This protocol is recorded before the new development metric outcomes. It defines
an experiment and referee criteria; it does not change or replace the frozen,
submitted `m4_multisource_dual_v1` model. Root owns scoring and model decisions.

## Question and fixed arms

Can increasing cellular dispersion improve official-definition DE direction
fidelity and significant-set Jaccard while preserving the other four metrics?

Use one incumbent prediction and two transformations of that same integer
matrix. The requested multipliers **1.5 and 2.0 refer to aggregate normalized-CPM
residual standard deviation**, corresponding to variance factors 2.25 and 4 if
attained. Identity multiplier 1 must return the input bit-for-bit. No additional
arms or tolerance changes are permitted after viewing this comparison's metrics.

The transform is an experimental distributional model. Seeded tie-breaking or
new heterogeneity is not evidence of biological realism. Record its rule, seed,
requested multiplier and actually attained dispersion. An unattainable requested
multiplier is a partial or infeasible transformation, not permission to weaken
the invariants. The implementation may preserve only aggregate spread rather
than each gene's spread; report that distinction and per-gene change summaries.

## Development data and information boundaries

Authoritative frozen selection:
`evaluations/dispersion_h1_training_400/cohort_freeze.json`.

- Source: public H1 2025 `adata_Training.h5ad`; seed `20261021`.
- 114 targets satisfy the metadata eligibility rule: at least 400 genuine cells
  and absence from the prior H1 Validation, H1 Test and K562 confirmation target
  lists. Select 50 development targets uniformly using the recorded seed.
- Sample 400 genuine reference cells per selected target, without replacement,
  plus 4,000 unique measured controls. Reference shape is 24,000 × 18,080.
- Every arm predicts 400 cells per target on the identical ordered 18,080-gene
  axis, giving 20,000 treated prediction rows. The exact same controls and
  reference are used for every arm.
- A disjoint 25-target reserve is selected now using metadata only. Its responses
  are not extracted or scored in this development experiment.
- `source_rows.csv` records source row/cell identities and has its hash in the
  freeze. The final manifest records reference/control file hashes.
- Exclude all H1 response sources from incumbent generation. Keep the frozen
  source weights, mean-effect settings and other modeling decisions; the new
  control draw and development seed are recorded, not selected by outcomes.
- Closed H1 Test and K562 confirmation response panels are not reopened for
  tuning, selection or additional scoring.

This is new Model4 development evidence, not a claim of previously untouched
biological data: H1 Training contributed to earlier source-bank work. H1 is also
the context used for previous model selection. Four thousand controls are
manageable and identical across arms, but differ from the official 18,400;
significance power can therefore differ from the challenge even though the
400-cell predicted group size now matches it.

## Exact and bounded matrix invariants

For each target independently, let integer counts be `X[i,g]`, fixed cell depth
`L[i] = sum_g X[i,g]`, raw gene total `C[g] = sum_i X[i,g]`, and
`mu[g] = mean_i(1e6 * X[i,g] / L[i])`.

Hard requirements, checked on the **final emitted integer matrix** against the
incumbent, not merely on an intermediate floating-point solution:

1. Same cell order, target labels, gene order and matrix dimensions.
2. Finite, nonnegative integral counts; no overflow; valid cell-depth and sparse
   storage limits. No control cells enter prediction files.
3. Every `L[i]` is preserved **exactly**, using integer arithmetic.
4. Every `C[g]` is preserved **exactly**, using integer arithmetic.
5. For every gene with positive incumbent mean,
   `abs(mu_new[g] - mu[g]) <= max(1e-8 CPM, 1e-4 * mu[g])`.
   This is a 0.01% relative mean-CPM bound with a numerical floor, not a claim of
   exact weighted-mean preservation. Zero-total genes remain exactly zero.
6. Identity is bit-identical; a fixed seed is reproducible; no reference
   perturbation outcome is accessed by the transformation.

A two-row, two-gene integer cycle adds `d` to `(i,a)` and `(j,b)` and subtracts
`d` from `(i,b)` and `(j,a)`. It exactly preserves row and raw column totals.
Its mean-CPM change for gene `a` is

`1e6 * d / n_cells * (1/L[i] - 1/L[j])`,

with the opposite change for gene `b`. Thus equal-depth cycles preserve CPM
means exactly, while nearby depths only bound the change. Track the **cumulative
per-gene drift** across every pass and independently recompute it at the end.
Per-swap bounds alone do not establish the final invariant. A low-count gene can
already be maximally dispersed at its fixed total, so universal per-gene
amplification is not guaranteed.

Under fixed real controls, the relative CPM bound limits an individual positive
gene's LFC change to approximately `1.443e-4` log2 units (with the scorer's fixed
epsilon). This is only an LFC bound. It does not itself prove that the aggregate
NMAE meets the strict metric nonregression criterion below.

## Which metrics are protected, and which must be checked

The pinned scorer is `cell-eval2 0.16.0`, commit
`5e64833518a6603a0301cbe28185d49c30f4a986`, using the `vcc2026` profile.

- Exact raw column totals and the fixed full panel preserve count pseudobulks
  and hence PDS. Verify actual PDS equality rather than relying on this argument
  alone; retain the same gene-exclusion and tie conventions.
- Arithmetic mean CPM determines DE LFC under fixed real controls. It is bounded,
  not exact here; inspect realized LFC drift and actual NMAE.
- Scored expression MSE includes a delete-one-cell sampling correction. Fixed
  pseudobulks alone do not fix this correction. With fixed means and cohort its
  numerator subtracts `min(B,Q)`, where `B` is the fixed across-perturbation bulk
  spread budget and `Q` is the weighted sum of prediction corrections capped by
  real corrections. MSE stays fixed if `min(B,Q)` stays fixed, for example if the
  budget remains binding. This condition must be measured, not assumed.
- Wilcoxon significance depends on cellular ranks/ties and distributional shape,
  not variance alone. Adjusted p-values use per-target BH on the fixed
  real-control CPM gate. The scorer computes p-values from submitted counts;
  post-hoc editing of DE tables is not an experiment arm.
- Reach ranks by predicted significance/confidence. It can change even when
  means, PDS and LFC do not, and is a protected metric in this test.

## Scoring and preregistered referee decision

All arms use the same reference, controls, complete 50-target panel, CPU/scanpy
backend, real-control convention and official preset. Keep the implementation
and wrappers fixed while jobs run. Scoring-only files may append the exact
measured controls to satisfy the installed scorer's input gate; never modify
the original prediction files. Export DE tables and metric aggregation cohorts.

Report the **six paired raw metrics only**. No H1 Validation anchors, new
cross-panel scaling or arithmetic average of heterogeneous raw metrics is used.

Strict development decision, relative to the incumbent:

1. Every hard matrix invariant passes.
2. PDS aggregate and available per-target values are unchanged within absolute
   tolerance `1e-10`.
3. All four protected metrics—PDS, normalized expression MSE, LFC NMAE and
   direction reach—are nonworse at the point estimate up to numerical tolerance
   `1e-8` in their own raw units. Lower is better for MSE/NMAE; higher is better
   for PDS/reach. Even a small genuine tradeoff gives **FAIL STRICT**.
4. Both fidelity-yield raw and Jaccard have positive point-estimate gains.
5. Before any future promotion, both gains require paired target-bootstrap
   confidence intervals above zero. Use 2,000 paired target resamples, fixed seed
   `20261021`, and preserve each metric's actual finite real-side cohort. Report
   intervals and win/loss/tie counts, not only aggregate gains. If either interval
   crosses zero, the result remains inconclusive for promotion even if strict
   point-estimate dominance holds.

Report uncertainty for protected metrics separately; a nonworse point estimate
is not statistical proof of no harm. Do not turn the mean-CPM engineering bound
into an extra metric degradation allowance. Do not increase tolerances or remove
targets after observing results. MSE is a ratio of sums: do not bootstrap it as a
mean of per-target normalized MSE values. Its correction budget also uses the
whole panel, so a naive subset rescore changes the metric.

This first-stage single-context development experiment does not promote a model
or authorize a submission. A promising arm can motivate a separately recorded
confirmation on the untouched 25-target reserve; it does not reopen closed
H1 Test or K562 confirmation panels.

## Required artifacts and implementation review

Preserve incumbent and transformed hashes, source/row provenance, exact
transformation arguments, seed, input/output margins, final CPM/LFC drift,
requested/achieved dispersion, per-gene variance-change summaries, runtime,
sparse density, all six raw metrics, actual finite cohorts, DE-call diagnostics,
paired uncertainty and authoritative process statuses.

The independent referee inspects `scripts/constrained_dispersion.py` and its
tests without editing the implementation. Required cases include identity,
equal/unequal depths, empty or infeasible movement, zero/low-count genes,
cumulative budgets, tight rounding bounds, deterministic seeds and final
recomputation independent of incremental drift bookkeeping. The signed review
goes in `reports/dispersion_referee_review.md`; it must distinguish engineering
correctness, feasibility, observed metric dominance and biological evidence.
