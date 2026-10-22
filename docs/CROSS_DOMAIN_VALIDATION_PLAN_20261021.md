# Proposed cross-domain validation after P2

This is a future proposal, not an executed experiment or an amendment to the frozen P2 benchmark. P2 has now finished: it recovers more magnitude-proxy responses but loses precision and increases error, so it fails advancement. None of the steps below by itself qualifies a model for new H1 response scoring, access to the 25-target reserve, or the user's conditionally authorized official submission.

## 1. Balance expression features as well as depth

P0 is informative about what to balance. Its exact depth-mean sampling variance fell by at least **99.5%**, yet arithmetic-CPM mean variance fell only **1.6–2.5%** across contexts. H1 bulk log-profile sampling error fell about **1.1%**, missing the frozen 5% requirement. These measurements suggest that variation along other expression directions dominates much of the remaining profile error; balancing depth alone does not balance those directions. This is an inference from P0's results (local artifact omitted: `../reports/sampling_precision_results.json`), not proof that another sampler will succeed.

The next sampling candidate would adapt the survey-sampling [cube method](https://doc.rero.ch/record/296198/files/914893.pdf) to select representative **actual cells**, with known equal inclusion probabilities. Freeze one design, rather than a feature-count search:

- Use the existing 4,000-cell control populations and select 400 donors. Set every cell's inclusion probability to **0.1**, with fixed sample size and unit weights on emitted cells.
- Define auxiliary features from controls alone: an intercept, scaled library depth, and the first four principal-component scores of standardized arithmetic CPM on the control-mean-CPM>5 gene set. Freeze feature preprocessing and numerical rank handling before drawing samples. No response labels, source-effect outcomes or P2 residuals select these features.
- Compare this one multivariate design with uniform sampling and the frozen depth-only design, using the same predeclared draw seeds. Preserve each selected raw cell; do not repair empirical gene totals or average donors.

Equal inclusion must be proved for the chosen implementation and checked on small exhaustive examples and repeated draws. Exact auxiliary balance may be infeasible for a finite population: report the residual and the algorithm's behavior, rather than silently changing inclusion probabilities. Correct marginal inclusion also does not make selected cells independent. Repeat the **actual policy's** disjoint-control null test and the sampling-correction diagnostics.

Use P0's original precision and null requirements as the proposed advancement standard, fixed again before execution. Report arithmetic-profile and bulk-ratio errors separately: unbiased linear means do not imply an unbiased finite-sample pooled ratio or log profile. Failure stops this sampling candidate; it does not justify selecting extra PCs after viewing results. Sampling precision cannot repair an incorrect biological response direction.

## 2. Inject known signals using independent truth

Kepler's [injection–recovery study](https://www.ipac.caltech.edu/publication/2015ApJ...810...95C) supplies the useful principle: insert a known signal into a realistic background and test the complete recovery pipeline. Our adaptation would diagnose emission and detection separately from the source-response model.

The critical safeguard is **independent ground truth**. Do not generate the reference with the emitter being tested, its fitted outputs, or the same thinning/birth transport code. Supplying an oracle mean to an emitter is acceptable; calling that emitter's generated cells the oracle reference is circular.

Use two predeclared truth families with independent generation code:

| Family | Construction and exact truth | What it can establish |
| --- | --- | --- |
| Real-control composition changes | Learn a small set of state bins on inference controls; assign disjoint reference cells with that frozen rule. Change bin mixture weights by fixed amounts. Calculate population arithmetic means and expected raw-count totals from the finite reference population. | Recovery of composition-driven expression changes while retaining real cell backgrounds. This is not a general model of cell-intrinsic knockdown effects. |
| Independent count simulation | Draw reference cells from a separately specified latent-state count distribution, such as a Dirichlet-multinomial family with declared depth dependence and probability shifts. Derive population moments from that generative specification. | Stress tests for signed effects, low abundance and activation at donor-zero entries. Synthetic realism remains limited. |

Freeze true-null, weak, moderate and strong signal regimes, reference populations and seeds before testing. Use disjoint inference, comparison and reference cells wherever real cells are involved. Hold out some signal templates and seeds from emitter development. Define signs and magnitudes **after full-axis normalization**; genes outside an edited list are not automatically unchanged after compositional closure.

Feed the same known population profiles to each emitter, bypassing the biological response predictor. At 400 cells per group, report realized arithmetic/bulk error, sign recovery by effect-size stratum, activation recovery, and coverage of any claimed profile intervals. Run the pinned DE pipeline and report false calls on **genuine null scenarios**, recovery on signal scenarios, and which true signals the control-expression gate excludes. A nonzero population mean effect does not guarantee significance in one finite sample. The known bulk quantity is a ratio of expected counts, not automatically the expectation of a sampled pooled ratio.

A failed oracle-recovery test points to emission, finite sampling or detection limitations. Passing it only establishes performance on the declared synthetic/composition families; it supplies no evidence that a source model predicts real perturbations.

## 3. Use P2 to choose the next question

Keep the [P2 gate](MATCHED_RESPONSE_P2_PROTOCOL_20261021.md) unchanged; it failed. Use known-signal recovery to distinguish emission/detection limitations from the observed source-model limitation. If recovery succeeds, prioritize matched response labels, context-transfer information and model assumptions. If recovery fails too, isolate that additional emission/detection problem before combining more components. Another dispersion or sampling adjustment cannot supply missing response signs.

Weather-style [ensemble copula coupling](https://arxiv.org/abs/1302.7149) is a later candidate if marginal recovery succeeds but dependence remains inadequate. A control-derived rank template would still need integer-count, row-depth, donor-zero activation and exact zero-effect identity contracts. Unconstrained reranking can violate those contracts. This proposal does not implement that method or attribute it to a challenge leader.

All follow-ups require a separate prospective freeze. These diagnostics can identify a failure mechanism; they cannot promote a model, reopen the reserve, or replace independent biological validation.
