# Next iteration after the published Model4 result

Prepared 2026-10-21 after the October 20 submission. **Planning only: no experiments are launched by this document.** Every proposed improvement below is a hypothesis. Preserve the submitted model, its artifacts and its selection record as the fixed incumbent.

Execution update: the user subsequently authorized N1 at 400 cells per target.
The separate [preregistered protocol](DISPERSION_EXPERIMENT_20261021.md) governs
that run. It is complete: both transformed matrices pass engineering checks but
fail the primary metric objective. Fidelity/overlap point estimates decrease;
protected point estimates are preserved or favorable. Neither arm is promoted.
See [the full referee result](../reports/REFEREE_DISPERSION_20261021.md).
This tests dispersion at 400 cells and does not isolate a 200-to-400 cell-count effect. N0 and N2–N4 remain proposed. The
[top-10 referee audit](../research/current_top10_referee_20261021.md) also identifies
bulk-expression accuracy and reach as larger competitive gaps than overlap alone.

## Evidence and the question to answer

The official validation entry `[PRIVATE_IDENTIFIER_REMOVED]` is published at **0.1470186853371711**, rank **565 when checked**. Its scaled components are:

| PDS | MSE | LFC NMAE | Direction fidelity | Direction reach | Significance Jaccard |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.645588 | 0.037668 | 0.138194 | −0.018683 | 0.079313 | 0.00003234 |

These numbers support prioritizing DE calibration and expression accuracy. They do not identify a particular failing context, gene set or mechanism: the saved result contains aggregate scores, not hidden response labels. Fidelity is below its reference baseline and Jaccard is approximately at baseline, while perturbation discrimination carries much of the positive score. Official receipt (local artifact omitted: `../submissions/m4_multisource_dual_v1_status.json`), [final report](../reports/FINAL_REPORT.md).

The local evidence constrains which experiments are worth repeating:

- The complete-source model led the 12-arm H1 Validation comparison at **0.242756**, versus 0.224969 for Orion+CD4. These use local five-split anchors and are not estimates of the hidden challenge score.
- Hard source-agreement gating lost to its matched Orion baseline: **0.210776 versus 0.219109**, despite higher PDS. Half amplitudes and the existing single-profile arm also lost. Do not reinstate these settings merely because the official DE scores are weak.
- Pool 1 lost to pool 4 in the PIE experiment, but that comparison also changed donor selection and depth distribution. It did not isolate dispersion.
- Frozen H1 Test and K562 confirmations improved all six raw metrics against their respective neutral emitters. Those comparators were not the official generic-response baseline. Both confirmations used 200 predicted cells, whereas production used 400. The H1 Test report also records excess DE calls; this motivates a diagnostic, not a target call-count setting.

Evidence: [all candidates](../reports/h1_all_candidate_comparison.md), [H1 confirmation](../reports/h1_test_confirmation.md), [K562 confirmation](../reports/k562_confirmation.md), [emission review](EMISSION_REVIEW.md).

## Evaluation rules before another search

1. Keep the existing H1 Test and K562 confirmation results as a closed historical audit. **Do not choose amplitudes, dispersion, thresholds, weights or target-specific corrections from their outcomes.** Do not repeatedly rerun candidates against those references.
2. Register new development and confirmation cohorts before inspecting their response outcomes. Use public H1 Training and new public Orion cohorts for development, selected from metadata with fixed seeds. Reserve new target/cell cohorts first; exclude previously reserved target lists from new development selections where feasible, and document overlap otherwise. Exclude every representation of the destination response family, including raw and PIE variants, in each context holdout.
3. Use 400 predicted cells per target for the new calibration work. Prefer real groups with at least 400 distinct cells and up to 18,400 distinct controls where available. Otherwise retain the genuine available counts and report the power difference; never duplicate reference cells to manufacture a 400-cell experiment.
4. Use one pinned official `vcc2026` implementation/backend and fresh context-specific baselines/anchors. Preserve original source normalization and missingness. K562's 8,246-gene panel, H1's 18,080 genes, Orion's full-library CPM denominator, and CD4's published LFCs are not interchangeable expression scales. Do not normalize stored means as a substitute for normalizing individual cells.
5. Freeze a small comparison and its success rule before execution. Begin with one seed; repeat only surviving candidates and the incumbent with two additional fixed seeds. Do not run a Cartesian product of all knobs below.

Report all six metrics, their finite-target cohorts, coverage strata, DE precision/recall and direction errors, realized arithmetic/bulk errors, depth distributions, detection rates and nonzero counts. Raw metrics with different units must not be averaged. Aggregate-derived MSE must be recomputed using the official aggregation rather than treated as a mean of per-target CSV values.

## Ordered, finite experiments

### N0 — Diagnose the 200-to-400-cell and control-count difference

**Hypothesis:** part of the DE weakness arises from distribution/calibration differences that become more consequential with production-sized groups; it is not established that this explains the official result.

On development data, run the frozen model and its neutral emitter at 200 and 400 predicted cells. Keep the underlying effects, source families, seed policy and reference cohort fixed. Separately run zero-response control diagnostics against disjoint measured controls, preserving guide/batch strata where available. Compare the generator's DE-call distribution with genuine control-versus-control comparisons, rather than demanding zero calls in heterogeneous real controls.

**Control and success:** this stage selects no new model. It must distinguish a cell-number effect from a reference/control-cohort change and expose any realization error. Archive the resulting diagnostics before choosing N1's settings. Do not infer DE precision from the A/B/C controls alone: they provide a null diagnostic, not perturbed ground truth.

### N1 — Change dispersion while holding response means and depths fixed

**Hypothesis:** four-donor pooling suppresses marginal variability and changes significance calls beyond what the predicted mean effects justify.

Use the same donor pools, exact cell depths, source effects, requested arithmetic means and requested pseudobulk for every arm. Compare residual-spread multipliers **1.0, 1.5 and 2.0** relative to the current pooled template. Refit the same two moments and apply the same integer rounding policy. The 1.0 arm must reproduce the incumbent procedure. Record clipping and feasibility projections; a failed fit is not evidence against the biological hypothesis.

**Success:** improve the complete six-component development score with a supported fidelity/Jaccard benefit, rather than merely reducing the number of calls. Reject improvements bought through a clear loss of reach, NMAE, PDS or MSE outside the incumbent's measured emission-seed variability. If refitting materially changes the realized mean/bulk profiles, treat the experiment as confounded and correct it before interpretation.

### N2 — Test cell/depth dependence after the dispersion comparison

**Hypothesis:** a realistic relationship between expression, cell state and library depth can improve count marginals without changing the predicted response means.

At N1's retained setting, compare the existing pooled-control residual structure with **one** control-derived, depth-stratified residual template. Use identical donors, row depths, desired moments and gene-wise dispersion targets. A gene-wise residual permutation within depth strata is a diagnostic negative control only if its post-fit marginals and moments remain matched; otherwise it cannot isolate dependence.

**Success:** a reproducible gain under the same six-metric rule, with matched mean/bulk realization and no extra false-direction burden. Most scored statistics are gene-wise or pseudobulk statistics; do not claim that realistic cross-gene covariance itself earns points. The relevant question is whether dependence changes the evaluated normalized distributions correctly.

### N3 — Calibrate DE effects without repeating the failed hard gate

**Hypothesis:** selective attenuation based on independently estimated source uncertainty is more useful than a binary agreement filter or uniform halving.

Keep the bulk amplitude at 0.3 and retain the best supported emitter. First compare arithmetic-head amplitudes **0.45, 0.60 and 0.75** on development data only. Then, only if diagnostics justify it, compare the retained global amplitude against **one** smooth confidence-shrinkage rule fit inside development folds. Candidate predictors may use source sample size, measured effect uncertainty and agreement across distinct source families; unavailable pairs contribute no votes. Do not treat raw/PIE versions of the same context as independent evidence, or CD4 conditions as separate families.

**Control and success:** the original 0.60/0.30 model remains the comparator. The confidence rule must improve responding-gene recovery or direction reliability in addition to the aggregate; its performance must not be a disguised prediction of no response. Do not set a fixed number of called genes from the reserved H1 Test count or a public competitor's setting. CD4's converted bulk effect remains an approximation and must stay labeled accordingly.

### N4 — Test a small context-weighting model

**Hypothesis:** a limited, control-conditioned change in source weights transfers better than the current fixed weights.

Compare exactly three policies: the incumbent fixed weights; the previously tested HCT-heavy alternative; and a similarity-weighted policy with a **single temperature parameter, restricted to 0, 1 or 2**, selected in nested development context folds. Use measured source/destination control profiles on explicitly shared genes, maintain a source-weight floor for coverage, and exclude the destination family entirely. Keep CD4 at its fixed weight unless an audited source control profile is acquired; its DE table alone is not a basal expression profile. Freeze the common-response centering population so adding a new target panel does not silently change the method.

**Success:** positive gains across at least two development contexts, with no dependence on one target subset and no material coverage loss. The earlier HCT-weight increase gained only 0.001410 locally; this is not evidence of stable superiority. If agreement gating is reconsidered, validate it with the exact final source-family mixture: the 0.75 rule changes its required vote count as families are added.

## Promotion and stopping rules

For a surviving candidate, compare paired target-level outcomes and repeat emission seeds. Use paired uncertainty estimates conditional on each fixed benchmark panel; recompute nonlinear aggregate metrics correctly. Treat a gain whose uncertainty includes no improvement as inconclusive, not as proof. The intended DE improvement must be visible, not solely a PDS gain hiding a DE regression.

Freeze one candidate, its code, source hashes, centering population and configuration before a **new** reserved confirmation. Evaluate it once against the incumbent on the newly reserved cohorts. If confirmation is mixed, retain the incumbent and report the limitation; do not tune on that confirmation. Only a materially different, locally supported candidate warrants another official validation submission. No numerical leaderboard gain is promised.

## Final contexts: verified prerequisites

Arc's official announcement schedules final controls for **October 22, 2026** and final submissions for **November 5 at 23:59 UTC** (**November 6 at 07:59 HKT**). Recheck the live deadline and dataset availability before acting. [Arc announcement](https://arcinstitute.org/news/virtual-cell-challenge-2026).

The current official CLI guide says final submissions require **VCC 0.3.0 or newer**. Today's installed 0.1.0 client must therefore be upgraded and checked before final packaging. Download `controls-final` separately; its manifest identifies **D/E/F** and its own target panel. Preserve those labels and supplied axes. Final scores are blind, only the last final entry counts, and the daily quota remains two submissions with one in flight. [Official CLI guide](https://vcc-cli-wiki.virtualcellchallenge.org/).

Before final generation, rebuild source coverage for the new target identities from verified public data, retaining the frozen normalization and centering policy. The existing 575-target banks are not evidence of final-panel coverage. Audit missing pairs, library/detection distributions and moment feasibility using final controls only; do not infer hidden response outcomes or relabel contexts by guessed cell identity. Re-run complete count/axis/hash checks and strict current-CLI prep, verify the final partition, and preserve the scored package and receipt. Confirm acceptance/terminal status rather than expecting a visible final score. Do not replace the chosen final artifact with a diagnostic or unchanged resubmission.

## Priority after the referee review

The top-10 gap and failed uniform-dispersion comparison motivate the separate
[bulk-response referee plan](BULK_RESPONSE_REFEREE_PLAN.md). It proposes four
controlled mean-response arms and distinguishes uncorrected mean error from the
scorer’s sampling correction. The user subsequently authorized this four-arm study. It is complete and all
three variants fail the strict no-regression criterion; no reserved confirmation
is triggered. See [the result](../reports/BULK_ABLATION_FINAL_20261021.md).

## Direction after the zero-response test

The [completed null diagnostic](../reports/EMISSION_NULL_FINAL_20261021.md) now
supports prioritizing the [structural model redesign](NEXT_MODEL_REDESIGN_20261021.md)
over additional global scaling searches. This supplements the historical plan;
it does not revise past experimental verdicts or claim a new accepted model.
