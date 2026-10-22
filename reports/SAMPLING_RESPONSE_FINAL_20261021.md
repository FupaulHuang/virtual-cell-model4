# Sampling precision and response recovery: completed P0–P2 study

**Decision: retain the incumbent. No new H1 perturbation scoring, reserve extraction
or official submission is justified by these results.** The user's conditional
authorization for one sufficiently supported submission remains unused.

The next-work plan was implemented through its two independent branches and the
single conditional response-aware experiment. All registered runs completed;
the failures below are advancement failures, not software execution failures.

## P0 — Depth stratification is valid but too weak

Both policies select 400 actual donors from the same 4,000-cell population.
Uniform sampling draws without replacement. Stratified sampling stable-sorts
library depths into 400 blocks of ten and chooses one cell per block. Both have
inclusion probability 0.1. Counts are never averaged or repaired.

There were 200 fixed draws per policy in each of H1, A, B and C: 1,600 draws total.
Ratios below are stratified/uniform error; values below one favor stratification.

| Context | Arithmetic log-profile MSE ratio | Bulk log-profile MSE ratio | Exact arithmetic-mean variance ratio |
|---|---:|---:|---:|
| H1 | 0.9964 | 0.9889 | 0.9790 |
| A | 0.9988 | 0.9894 | 0.9820 |
| B | 0.9911 | 0.9991 | 0.9748 |
| C | 0.9999 | 0.9988 | 0.9836 |

H1 bulk error improves about **1.11%**, with ratio interval 0.9813–0.9963.
This is a modest positive effect, but it misses the predeclared **5%** requirement.
The protected precision checks pass. Exact depth-mean sampling variance falls by
at least 99.5%, yet arithmetic-profile mean variance falls only 1.6–2.5%.
Thus library depth captures little of the normalized gene-profile sampling error.

The exact-policy null study passes: all 40 stratified groups have zero significant
calls; uniform groups have zero or one. All 80 groups reproduce selected raw donors
exactly. The four contexts have disjoint comparison controls and the same frozen
population margins; per-draw guide/batch imbalance is reported without changing
the sampler. These are correlated technical groups, not independent biological
replicates or an empirical FDR guarantee.

See [precision results and intervals](sampling_precision_results.md),
[null results](stratified_sampling_null_results.md), and the
independent DE-table audit (local artifact omitted: `stratified_sampling_null_table_audit.json`).
Passing the null check does not override the failed precision gate. No new
nonzero-response sampling candidate was generated or scored.

## P1 — The matched comparison confirms response under-recovery

The audit found no coherent matched DE-label set for all three raw source banks.
K562 has no corresponding q-value table; PIE HCT116/HEK293T DE uses different
perturbation-cell populations and control policies from the raw banks. This is a
population/estimand mismatch, not proof that the published DE data are wrong.
See the [source-label audit](matched_response_source_audit.md).

Consequently, this experiment uses explicitly declared **magnitude proxies**:
absolute response at least 0.25 in log2(CPM+1) units or 0.10 in natural-log bulk
units. They are not DE significance, biological truth or official challenge scores.

Each source context was excluded from every training input and label in its fold.
Zero response, matched fusion and adapter-v1 used the same source families, priors,
units, prediction clipping and full-axis inverse. CD4 was omitted from all three.
Evaluation sampled 100,000 target-balanced coordinates per context, retaining
uncovered coordinates and excluding each target's own gene. No H1 response was used.

Adapter-v1 made **zero proxy calls in every context and both heads**, despite
source coverage of about 89.6%, 99.8% and 99.0%. It reduced squared error but failed
to recover responses found by matched fusion. The preregistered under-recovery
trigger fired, permitting one new response-aware procedure. The comparator here
is equal-family matched fusion, not the published production configuration.

See the [P1 results](matched_response_results.md) and
[frozen protocol](../docs/MATCHED_RESPONSE_BENCHMARK_PROTOCOL_20261021.md).

## P2 — More calls did not mean better information

The single frozen hurdle procedure fits response occurrence, conditional sign and
magnitude separately, with within-training-source probability recalibration.
The response decision uses the declared class-balanced prevalence boundary,
not posterior probability 0.5 or a selected threshold. No call quota or rescue
variant was used. The same original fusion-based advancement gate was retained.

Equal-context summaries below remain separate for each response head:

| Head | Predictor | Signed precision ↑ | Signed recovery ↑ | Class-unweighted MSE ↓ |
|---|---|---:|---:|---:|
| Arithmetic | Matched fusion | 0.1991 | 0.1281 | 0.12085 |
| Arithmetic | Hurdle P2 | 0.1639 | 0.3611 | 0.20916 |
| Bulk | Matched fusion | 0.1473 | 0.0732 | 0.00903 |
| Bulk | Hurdle P2 | 0.0826 | 0.3240 | 0.02138 |

P2 recovers more proxy responses but loses precision and increases error in both
heads. The harm also appears on genes whose matched control mean CPM exceeds five;
it is not confined to low-expression coordinates. For example, HCT116 arithmetic
signed precision on that subset falls from 0.2058 to 0.1798, and bulk precision
falls from 0.1913 to 0.1466.

Probability calibration also fails to transfer reliably: HCT116 arithmetic mean
predicted responsivity is 0.4745 against observed proxy prevalence 0.2478; HEK293T
is 0.4457 against 0.1963. Training-source recalibration cannot be treated as evidence
of calibration in a new context. The calibration split is fit data with full-source
common preprocessing, not untouched target/context validation.

**Both heads fail.** No further objective, cutoff adjustment, or H1 scoring cycle
was run. This result applies to the whole registered hurdle policy, which changes
model structure, training losses, calibration and decisions together; it does not
isolate a single cause or establish that all response-aware models must fail.

See [P2 results](matched_response_p2_results.md),
full diagnostics (local artifact omitted: `matched_response_p2_results.json`),
[frozen protocol](../docs/MATCHED_RESPONSE_P2_PROTOCOL_20261021.md), and
successful execution receipts (local artifact omitted: `matched_response_execution_receipts.json`).

## What the cross-domain analogies contribute

Our assessment is that **weather forecasting is the closest overall benchmark
analogy**. [WeatherBench 2](https://research.google/pubs/weatherbench-2-a-benchmark-for-the-next-generation-of-data-driven-global-weather-models/)
and its [evaluation code](https://github.com/google-research/weatherbench2/blob/main/scripts/evaluate.py)
separate ensemble-mean accuracy, finite-ensemble corrections and distributional
spread. These are useful parallels to the mean/variability tradeoff observed here.

Three concrete ideas warrant separately frozen follow-ups:

1. **Survey/geographic balanced sampling:** use multiple control-derived expression
   features to select representative actual cells with known inclusion probabilities.
   The [cube method](https://doc.rero.ch/record/296198/files/914893.pdf) supplies a
   mathematical foundation. P0 tested depth stratification, not this full method.
2. **Astronomical injection–recovery:** feed known independent signal distributions
   through the emitter and detector to measure where recovery fails. The
   [Kepler study](https://www.ipac.caltech.edu/publication/2015ApJ...810...95C) illustrates
   this principle. The emitter under test must not generate its own reference truth.
3. **Weather ensemble calibration:** [ensemble copula coupling](https://arxiv.org/abs/1302.7149)
   offers a way to preserve dependence while calibrating marginals. Integer counts,
   row-depth limits and exact zero-effect identity would still need to hold.

These methods address sampling, emission and detection. They do not supply missing
biological response information. The matched-source evidence here calls for better
matched response labels and credible context-transfer calibration before another
candidate is considered for submission.

The [cross-domain research note](../research/CROSS_DOMAIN_METHODS_20261021.md) and
[concrete validation plan](../docs/CROSS_DOMAIN_VALIDATION_PLAN_20261021.md) describe
those proposals. No cube/ECC or injection–recovery experiment has been executed
in this study. The future plan now records P2's failure and separates those new
diagnostics from another source-response fitting attempt.

## Verification and limits

There were 31 passing tests: 12 sampling tests (including exhaustive finite-population
proof checks), nine P1 guard tests and ten P2 tests. Independent reviews covered the
sampling formulas, held-out-context exclusion, source provenance, probability weights,
fixed gates and failure handling. All four precision runs, eight null-DE jobs, P1
and P2 completed successfully. Frozen code/data identities and report links are
checked in the completion audit (local artifact omitted: `sampling_response_completion_audit.json`).

The three source contexts were already used during adapter-v1 development; they
are not untouched validation. Training folds have one input source while held-out
evaluation has two, and native source panels differ. Magnitude labels and descriptive
bootstrap intervals do not replace matched DE labels or biological validation.
The 25 reserve targets remain closed. The published model and its score of
0.1470186853 are unchanged by this work.
