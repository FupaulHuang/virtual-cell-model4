# Proposed next step: recover distribution changes as well as means

Status: **unexecuted proposal**. This is a follow-up to the completed [cross-domain study](../reports/CROSS_DOMAIN_FINAL_20261021.md), not an amendment to its frozen methods or results. Its motivation appeared in development and repeated in unchanged confirmation: all 12 nonzero synthetic programs show zero donor-zero → positive edits among known-up genes, while strong held-out cases recover only 23–45 of 64 effects versus 61–64 in independent references.

## Problem and why it matters

The current pooled emitter generates thousands of significant calls on zero-effect inputs. The preserving emitter avoids this behavior, but in the synthetic development programs it produces no donor-zero → positive edits among genuinely upregulated genes. Its population means fit accurately while its detection fractions remain near the control baseline. The independent reference increases both abundance and detection.

This exposes a limitation that an extra sampling PC or global dispersion adjustment cannot directly address: the same two population moments can arise from different distributions across cells. A mean increase concentrated in already-positive cells need not produce the same DE result as an increase spread into previously zero cells. The biological source model's response-prediction problem also remains unsolved.

## 1. Quantify the existing activation budget before changing it

Use the existing frozen requests and control donors to audit the preserving emitter's expected added count mass at donor-zero entries, separately for each gene and control state. Keep all observed outputs and old study reports immutable. Label this as a mechanism audit, not independent confirmation.

The existing transport permits an exact calculation conditional on the actual donor row and its **post-fit expected counts**. The total removed count is a sum of independent binomials. If deficit gene j has raw count x_j and removal probability d_j, and a surplus gene receives fraction r of all births, then a donor-zero entry's probability of remaining zero is:

`P(no birth) = product_j (1 - d_j * r) ** x_j`.

Here r is the gene's expected added count divided by the row's total expected surplus. A row with zero surplus has zero activation probability. Calculate the product stably in log space. Compare expected activation counts with observed counts, fitted moment residuals and donor-sampling error. Assess conditional uncertainty or use a preregistered transport-only simulation that retains dependence between genes: one low realization alone does not establish a transport bug. This can distinguish negligible allocation to zeros from a transport defect, but cannot identify per-cell biological counterfactuals.

Deliverable: a small table of expected/observed activation and count mass for the existing moderate and strong programs, with source hashes and a numerical check of the probability formula on a small enumerated example. Do not launch another parameter sweep.

## 2. Specify one gene-local response allocation

Investigate one replacement for the globally scaled support seed. Allocate a gene's requested positive change across donor cells using a control-derived expression/state model, including donor-zero entries. The available inputs remain the inference controls and the same arithmetic/bulk requests. Oracle detection fractions, reference cell-state weights and DE outcomes must not enter the candidate.

The key distinction is that added mass scales with that gene's requested change. It should not be suppressed merely because the total response occupies a small fraction of the transcriptome. Both heads must be respected: total expected count change is `sum(depth) * delta_bulk`, while the sum of depth-normalized expected count changes is `N * delta_arithmetic`. These changes can differ in sign. Specify a fixed rule for mixed-sign and infeasible cases; do not silently choose one head. The allocation model, control preprocessing and all fallback rules must be frozen before evaluation. A control-state weighting is a proposed assumption to test, not recovered perturbation biology.

Retain these contracts:

- An exactly zero request returns the selected raw donors bit for bit.
- Each output has nonnegative integer counts and exactly its donor's library size.
- Fit arithmetic and bulk moments on the full inference population before selecting 400 donors. Keep the existing projection and fit-error limits; reject infeasible cases rather than loosening limits after outcomes.
- Use actual uniformly selected donors. The failed balanced-sampling candidate remains excluded.
- Report allocation before fitting, fitted expected counts, selected-donor expectations, and actual integer outputs. Row normalization and fitting can change the proposed activation allocation; quantify that change.

The first implementation task is a feasibility preflight on small control panels: verify zero identity, conditional count expectations, both population moments and activation consistent with the available added mass. For a donor-zero entry, `P(output > 0) <= E[output]`; a tiny expected count cannot support a large detection increase. Keep the original emitter unchanged and version the new allocation separately. If the constraints erase the intended activation or require large distortions, stop this design before a full challenge run.

## 3. Validate one frozen candidate on fresh cases

The now-opened injection templates and seeds are diagnostic regression cases. They cannot serve again as unseen confirmation after this design choice.

Before generating new outcomes, freeze new signal templates, seeds, 400-cell sampling, reference-null limits, effect bins and a paired comparison with the unchanged preserving emitter. Use independently generated reference cells and controls; the reference must not use the new allocation/transport code. Include an independent oracle-law replicate to measure finite-sample power.

The evaluation must jointly report signed recovery, false calls on exact NULLs and unchanged synthetic genes, effect-size and abundance strata, detection changes, and arithmetic/bulk errors from actual counts. Set a minimum useful improvement and nonregression margins in the prospective protocol before execution; no new numerical acceptance threshold is asserted by this proposal. A larger DE list alone cannot pass.

A candidate that improves detection while damaging null behavior or the existing profile constraints stops here. Retain all failed cases. No selection of favorable seeds, effect sizes or templates is permitted.

## 4. Reconnect to biological prediction only after the diagnostic passes

Rebuild response effects and DE-confidence labels jointly from the same source cells, controls and batch definitions. The [source audit](../reports/matched_response_source_audit.md) showed that the existing raw banks and published DE tables are not matched populations. Begin with a metadata-selected non-H1 panel and keep its destination context out of training inputs and labels.

Compare the frozen source fusion and existing response adapter before emission; then combine a qualified predictor with a qualified emitter. The existing six-metric nonregression checks still apply. Synthetic activation recovery supplies no guarantee of better official direction fidelity or significance overlap and no guarantee that other metrics will be unchanged.

The user's conditional authorization for one official submission remains unused. A submission requires evidence of a sufficiently enhanced biological model under the established validation procedure; this proposal is not that evidence.
