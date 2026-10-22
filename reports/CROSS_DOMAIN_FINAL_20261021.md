# Cross-domain sampling and injection–recovery study

This study implements the [cross-domain validation plan](../docs/CROSS_DOMAIN_VALIDATION_PLAN_20261021.md) at **400 cells per group**. Both branches and the immutable development/confirmation analyses are complete. Balanced sampling improves precision modestly but fails its advancement threshold. Independent injection–recovery identifies an additional distribution-recovery limitation: the preserving emitter avoids massive zero-effect DE calls, yet misses detection increases and many known effects. **No model is promoted or submitted.**

## Balanced sampling: a real but insufficient engineering gain

We implemented a randomized cube-flight/pivotal-landing sampler with fixed sample size and equal marginal inclusion probabilities. Each draw selects 400 unmodified cells from 4,000 controls. The frozen auxiliary variables are library depth and four control-derived expression PCs, plus an intercept. Exact auxiliary balance is not always possible; the implementation records its residuals and numerical tolerances. Equal inclusion does not make selected cells independent.

All 200 prespecified draws were retained in each context. The ratios below compare balanced with uniform sampling; lower is better. These are **profile sampling errors**, not official challenge metrics.

| Context | Arithmetic log-profile error ratio, 95% interval | Bulk log-profile error ratio, 95% interval | Precision gate | Control-only DE gate |
| --- | ---: | ---: | --- | --- |
| H1 | 0.9845 (0.9796–0.9893) | 0.9587 (0.9518–0.9652) | Fail | Pass |
| A | 0.9876 (0.9821–0.9934) | 0.9624 (0.9558–0.9687) | Pass | Pass |
| B | 0.9885 (0.9777–0.9994) | 0.9605 (0.9518–0.9690) | Pass | Pass |
| C | 0.9860 (0.9802–0.9915) | 0.9386 (0.9283–0.9487) | Pass | Pass |

H1 bulk-profile sampling error fell **4.13%**, short of the frozen **5%** requirement. The combined engineering gate therefore fails. We did not change the threshold or add PCs after seeing the result. Balanced sampling was excluded from the subsequent injection experiment before any recovery results were opened.

*Sampling error ratios and paired intervals (figure excluded from the source-only export).*

The [separate correction audit](balanced_sampling_correction_audit.md) shows why nominal independent-cell variance estimates need care: nominal variance/n exceeds the observed arithmetic-mean error by factors of 1.33–1.50 for balanced sampling, versus 1.10–1.12 for uniform sampling. Uniform sampling already has a finite-population correction. This diagnostic is **not** the official transformed bulk/global-capped correction and does not establish official score inflation.

Evidence: [frozen protocol](../docs/BALANCED_SAMPLING_PROTOCOL_20261021.md), results (local artifact omitted: `balanced_sampling_results.json`), eligibility decision (local artifact omitted: `injection_recovery_eligibility.json`), [implementation](../scripts/balanced_donor_sampling.py).

## Independent injection–recovery design

This branch bypasses the biological response predictor and supplies known arithmetic-mean and bulk-probability response deltas to the current pooled emitter and the control-preserving emitter with uniform donor sampling. Both receive identical requests. It tests whether the emission and detection pipeline can recover declared effects; it does not test whether the biological model can predict those effects.

Truth comes from two independent constructions:

- **Composition:** change the mixture of frozen control-derived cell-state bins, then sample actual cells from a disjoint real-control reference population. No emitter is used to generate reference cells.
- **Dirichlet-multinomial:** independently simulate cells from a frozen four-state/depth distribution with fixed concentration 20,000 and two 64-gene response programs. The simulator does not use the candidate count-transport implementation.

Each family has two templates, four regimes (NULL, weak, moderate and strong), and two seeds. Development is template 1/seed 20261026 only. Confirmation separately holds out template, seed, or both. The same frozen methods and analysis code apply to both stages.

Requests add the known probability delta to the empirical inference-control baseline, clip negative probabilities if necessary, and normalize. An exact zero delta returns the inference profiles unchanged. This preserves the meaning of a zero response while retaining finite baseline differences as nuisance. Candidates receive two moments, not complete target distributions or cell-state mixture weights. A recovery gap can therefore reflect insufficient distributional information, generator assumptions, finite sampling, baseline differences or numerical error.

The DE analysis uses the pinned cell-eval2 0.16.0 `compute_de` settings: arithmetic CPM, comparison-only mean CPM >5, epsilon 1e−9, per-group BH adjustment and q<0.05. Signed recovery requires a true absolute log2 fold change of at least 0.1. Tiny effects are not treated as null. NULL checks and unchanged conditional-law genes in the synthetic family are reported separately. This is not the official six-metric challenge score.

All 8 candidate arms and 64 alias groups completed under unchanged fit limits. All 8 control-preserving NULL groups exactly replay their selected raw donors. The audit verifies the profiles calculated from **actual integer output counts**, rather than accepting requested or fitted profiles as realized outputs.

An initial truth-file encoding problem stopped all four first attempts before cell generation. The gene-name arrays were repacked as NumPy Unicode; all numerical arrays were verified unchanged. The serialization amendment (local artifact omitted: `injection_recovery_serialization_amendment.json`) records this repair. Final truth-generation jobs and count/source/depth validation passed.

Evidence: [truth protocol](../docs/INJECTION_RECOVERY_PROTOCOL_20261021.md), [analysis protocol](../docs/INJECTION_RECOVERY_ANALYSIS_PROTOCOL_20261021.md), truth execution (local artifact omitted: `injection_recovery_truth_execution.json`), candidate execution (local artifact omitted: `injection_prediction_execution_receipt.json`), analysis freeze (local artifact omitted: `injection_recovery_analysis_freeze.json`).

## Development results

The immutable [development report](injection_recovery_development.md) opens only template 1/seed 20261026. Both opened independent-reference NULL groups have zero calls. Its family/seed validity was correctly reported as provisional at that stage; confirmation subsequently completed and passed all reference-NULL checks.

| Development case | Independent reference | Current pooled emitter | Preserving emitter |
| --- | ---: | ---: | ---: |
| Composition NULL: significant calls | 0 | 4,428 | 0 |
| Synthetic NULL: significant calls | 0 | 10,438 | 0 |
| Synthetic moderate: correctly signed significant effects / 64 | 56/64 | 57/64 | 30/64 |
| Synthetic strong: correctly signed significant effects / 64 | 64/64 | 64/64 | 48/64 |

The current emitter's apparent recovery comes with 10,373 calls on genuinely unchanged synthetic genes in each nonzero development program. Call count or signed sensitivity alone would therefore give a misleading impression of quality.

For the preserving emitter, the population moments fit well and request clipping is negligible in the synthetic programs. Yet no upregulated donor-zero entries become positive in any of the three nonzero programs. In the moderate case, the independent reference's mean detection fraction across the 32 upregulated genes rises 7.06 percentage points over controls, while the candidate changes −0.92 points. In the strong case, these changes are +12.15 versus +0.05 points. These are sample-level detection differences, not paired per-cell activation ground truth.

This supports investigating how a response is distributed across cells, including entry into previously zero genes. It does not prove that activation alone causes every recovery gap: the candidates receive only two moments, baseline samples differ, and 400-cell donor sampling contributes error. The count audit also separates donor-sampling error from count-realization error; in the synthetic development programs the latter is much smaller in global probability L1. Global error size does not by itself identify the cause of a particular gene's DE result.

## Confirmation and decision

The [confirmation report](injection_recovery_confirmation.md) used unchanged code, requests, thresholds and emitters. Its recorded development hash matches the original immutable report. All four family/seed reference-NULL panels pass, with **0–1 calls per group**. Across the complete study, current-emitter NULL calls range from **4,361–4,433** in composition and **10,416–10,438** in simulation. Preserving-emitter NULL calls are **0** throughout composition and **0–2** in simulation.

The held-out synthetic cases reproduce a substantial recovery gap:

| Confirmation case | Regime | Independent reference: correct-sign calls / 64 | Preserving emitter: correct-sign calls / 64 |
| --- | --- | ---: | ---: |
| Held-out template, seed 20261026 | Moderate | 15/64 | 0/64 |
| Held-out template, seed 20261026 | Strong | 61/64 | 29/64 |
| Held-out seed 20261027, template 1 | Moderate | 48/64 | 26/64 |
| Held-out seed 20261027, template 1 | Strong | 64/64 | 45/64 |
| Held-out template and seed 20261027 | Moderate | 3/64 | 0/64 |
| Held-out template and seed 20261027 | Strong | 63/64 | 23/64 |

The held-out template selects genes with more donor zeros. In its strong cases, reference detection across upregulated genes rises **11.91 and 12.11 percentage points** over controls; preserving outputs change **−0.118 and −0.109 points**. There are **zero donor-zero → positive edits among known-up genes in all 12 nonzero synthetic programs**, including development and confirmation. This is a repeated engineering observation, not a claim that the independent reference identifies individual cellular counterfactuals.

An explicitly **post hoc** up/down split, independently reconstructed from the frozen DE tables and donor counts, sharpens that interpretation: in the two strong held-out-template cases, the preserving emitter recovers **0/32 upregulated genes**, versus **29/32 and 31/32** in the independent references. It recovers **29/32 and 23/32 downregulated genes**, versus 32/32 in both references. This is descriptive evidence for the next mechanism study, not a new selection gate or proof that the support prior is the sole cause.

Synthetic request-delta distortion is negligible (arithmetic-probability L1 below **6.65e−16**). Donor-sampling L1 errors are **0.0286–0.0304**, while conditional integer-realization L1 errors are **0.000030–0.000240**. Thus count realization is not the largest global moment-error term in these cases. This does not eliminate finite-sampling effects on individual DE results, nor prove that insufficient activation is the only distributional problem.

*NULL calls and signed recovery across all fixed cases (figure excluded from the source-only export).*

The figure includes every defined fixed-case ratio. Lines are descriptive means across template/seed cases, not confidence intervals or independent biological replications. Composition-reference power is often low, particularly for weaker effects, so its recovery percentages must be interpreted alongside reference power. The pooled emitter's high signed sensitivity is accompanied by massive NULL and unchanged-gene calls; it is not evidence of a successful generator.

The practical next step is to quantify and redesign how each gene's requested added count mass reaches donor-zero cells, while retaining exact zero-response identity and both requested moments. The existing global support prior is a candidate mechanism to investigate, not a demonstrated sole cause. These experiments do not establish improved official direction fidelity or significance overlap, and they do not establish nonregression of the other challenge metrics. Biological response prediction still requires matched source labels and independent validation.

## Validation and scope

The implementation has **43 passing focused tests**: 10 sampler tests, 14 independent-truth tests, 10 prediction/request tests and 9 staged-analysis tests. They cover sampling marginals and fixed size, independent truth moments and serialization, exact zero-input requests/replay, incomplete-arm handling, true-null classification and separation of development from confirmation. All 12 DE jobs completed with exit code zero.

The [completion audit](cross_domain_completion_audit.md) passed **994 checks across 224 hashed files**, including independent reproduction of the diagnostic summaries and direct donor/count checks. The published model/configuration and production emitter are preserved. This work uses H1 control cells, independent simulations and existing source controls; it does not evaluate actual H1 perturbation responses or the reserved response set. No official submission is made. The user's conditional authorization for one sufficiently supported submission remains unused.

Audit execution required operational repairs, documented in the restart receipt (local artifact omitted: `cross_domain_completion_audit_restart.json`), matrix-loading amendment (local artifact omitted: `cross_domain_completion_audit_io_amendment.json`), and membership-lookup amendment (local artifact omitted: `cross_domain_completion_audit_membership_amendment.json`). The audit now loads each count matrix once per seed and uses an indexed gene-name lookup. The original and replacement lookup selected byte-identical gene masks on an actual 18,080-gene case and a reordered/missing/Unicode case; lookup time fell from 14.51 to 0.00749 seconds in the actual case. Progress logging identified the repeated membership comparisons as a major delay. Abandoned audit attempts and source versions are archived. No frozen experimental code, model, DE output or test result changed, and no model or DE job was rerun for these repairs.

The proposed next work is recorded separately in [the distribution-recovery plan](../docs/DISTRIBUTION_RECOVERY_NEXT_PLAN_20261021.md). It begins with an exact conditional activation-budget audit, then tests one gene-local allocation mechanism under the same zero-response, integer-depth and dual-moment contracts. Fresh confirmation cases and matched biological response labels are required before any promotion. This proposal has not been implemented or evaluated in this study.
