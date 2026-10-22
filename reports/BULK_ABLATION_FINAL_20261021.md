# Completed bulk-response ablation

**All three variants fail the preregistered no-regression criterion.** Keep the
published Model4 incumbent. No recipe is promoted, no reserved target responses
are extracted, and no official submission is made from this study.

The experiment isolates raw-source bulk common-response subtraction and final
bulk amplitude while retaining the original arithmetic-CPM head, source files,
weights, donor template, depths, seeds and count generator. It uses the same
50-target H1 Training development cohort: 400 measured and 400 predicted cells
per target, 4,000 controls and 18,080 genes, with all H1 response sources excluded.
This is an already-observed development context, not independent confirmation.

## Registered contrasts and measured results

| Arm | Raw bulk common subtraction | Bulk amplitude | Purpose |
| --- | ---: | ---: | --- |
| Incumbent | 1.0 | 0.30 | Fixed comparator |
| common05 | 0.5 | 0.30 | Restore part of the shared raw-source bulk response; CD4 contribution fixed |
| bulk015 | 1.0 | 0.15 | Reduce final fused bulk response |
| bulk060 | 1.0 | 0.60 | Increase final fused bulk response |

All values below are local raw official-definition metrics, not leaderboard
scores. Lower is better for MSE and NMAE; higher is better for the other four.
No heterogeneous raw overall score or cross-panel anchor scaling is used.

| Metric | Incumbent | common05 | bulk015 | bulk060 |
| --- | ---: | ---: | ---: | ---: |
| PDS ↑ | 0.910204082 | 0.908979592 | 0.904081633 | 0.908979592 |
| Normalized MSE ↓ | 0.970348453 | 0.969871583 | 0.993564953 | 1.170655702 |
| LFC NMAE ↓ | 0.962901741 | 0.962755307 | 0.963253510 | 0.962400446 |
| Direction fidelity ↑ | 0.510809091 | 0.511297521 | 0.510522401 | 0.510937685 |
| Direction reach ↑ | 0.115049906 | 0.114257266 | 0.118680434 | 0.109785434 |
| Significant-set overlap ↑ | 0.092176112 | 0.092107842 | 0.092129562 | 0.091830550 |

- **common05:** MSE and true mean error improve slightly, but PDS, reach and
  overlap decrease. This is a tradeoff, not an improvement without side effects.
- **bulk015:** true mean error improves, but scored MSE, PDS, NMAE, fidelity and
  overlap worsen. Reach improves. The candidate does not qualify.
- **bulk060:** mean error and scored MSE worsen substantially, and PDS, reach
  and overlap also fall. The larger-amplitude hypothesis is unsupported here.

The point gate requires both lower raw normalized MSE and lower realized
uncorrected mean error, with all five other metrics nonworse within `1e-8`.
All three arms fail it. Small losses are not relabeled as wins after seeing the
results. Paired uncertainty describes the strength of the evidence; it does not
change that registered decision rule.

The smaller amplitude's PDS change is −0.00612245 with paired 95% interval
**[−0.01183673, −0.00163265]**. The larger amplitude's reach change is
−0.00526447 with interval **[−0.01221057, −0.00061073]**. Common05's small
protected-metric losses have intervals crossing zero; they are point-gate
failures, not established declines by those intervals. All metric-specific
finite target cohorts are identical across arms: 42 targets for NMAE and 50 for
PDS, fidelity, reach and Jaccard. Global MSE is computed from its complete
panel-level components.

## Mean response versus sampling correction

Uncorrected mean error is computed from `log1p(50000 × group gene-count fraction)`
on the complete gene axis, then averaged after removing each perturbation's own
gene. It measures predicted-versus-reference mean residuals without a sampling
correction.

| Arm | Realized uncorrected mean error | Relative change |
| --- | ---: | ---: |
| Incumbent | 0.002221432735 | — |
| common05 | 0.002220581416 | −0.0383% |
| bulk015 | 0.002149668070 | −3.2306% |
| bulk060 | 0.002582763318 | +16.2657% |

The paired 2,000-resample target-bootstrap intervals for mean-error improvement
are positive for common05 (**[2.83e−7, 1.44e−6]**, 33/50 targets improve) and
bulk015 (**[5.02e−5, 9.21e−5]**, 47/50 improve). Bulk060's interval is negative
(**[−0.00041977, −0.00030846]**, only 1/50 improves). These are conditional
development-panel intervals, not evidence of hidden-context generalization.

The smaller-amplitude result is particularly informative. Its uncorrected
error decreases by **0.000071764665**, but its effective prediction sampling
correction decreases by **0.000113148854**. The lost correction allowance is
larger than the genuine mean-error benefit. Consequently, the capped numerator
increases by **0.000041384189**, and normalized MSE worsens by **0.02321650**.
The fixed reference denominator is approximately 0.0017825335324.

This does not make the smaller amplitude a successful challenge model, but it
also does not mean its average expression prediction became worse. The model
must improve biology and satisfy the evaluated score together. The scorer's
across-perturbation correction budget is a fixed part of this experiment; no
correction rule, evaluation threshold or generated variance is changed to
offset the failed result.

The common05 MSE gain is almost entirely attributable to a lower mean residual:
its correction change is only about −1.28e−9. Bulk060's slight correction
increase offsets only a small part of its larger mean-error loss.

## Isolation and validation evidence

The new incumbent exactly reproduces all old counts and semantic axes, and
even has the identical H5AD SHA-256
`03fea48aca0cfad51a08fc12ca5be01b7c7a1430beb7f674a75f0b073d0ba3c6`.
Its existing metric results are reused only after this full identity check.

All four complete outputs passed count/axis/hash validation. Requested
arithmetic-CPM profiles are identical over all 904,000 target-gene pairs;
CD4 centering, donor template and row-depth hashes remain fixed. Variant bulk
gene totals change intentionally. Original projection and fitting limits are
unchanged; the maximum projection L1 is 0.000807 versus the 0.03 rejection limit.
The variants retain about 98.6–99.1% of their requested bulk intervention
magnitude. Realized arithmetic CPM does change through joint fitting and
integer rounding, and that change is recorded rather than assumed absent.

The builder passed ten meaningful tests. Six additional analysis tests cover
normalization, target exclusion, paired cohorts, aggregate-MSE handling and
undefined alignment. All three new CPU/scanpy metric processes exited 0.
Scoring copies, reference, cached real-DE content, source hashes and configuration
were verified. No frozen production code, model, source bank or previous score
artifact was modified.

## Evidence and disposition

- [Frozen experiment protocol](../docs/BULK_ABLATION_PROTOCOL_20261021.md)
- [Full paired analysis and uncertainty](bulk_ablation_results.md)
- [Independent referee decision](bulk_referee_review.md)
- [Requested/projected/realized profile audit](bulk_profile_diagnostics.md)
- [Builder evidence](bulk_implementation_evidence.md) and [wrapper/analysis review](bulk_protocol_wrapper_review.md)
- Exact incumbent identity (local artifact omitted: `bulk_baseline_identity.json`) and terminal statuses (local artifact omitted: `bulk_terminal_statuses.json`)

The authorized conditional confirmation is not triggered because no arm meets
the development gate. The 25-target reserve remains unextracted. The outcome
does not establish that no better bulk-response model exists; it rejects these
three isolated global adjustments under the stated requirement. Any later
context-specific or gene-specific response model needs a new finite protocol,
rather than choosing further amplitudes from these development outcomes.
