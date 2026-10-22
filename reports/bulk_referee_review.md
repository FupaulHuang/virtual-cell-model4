# Independent bulk-ablation referee review

## Engineering verdict: PASS for the registered four-arm experiment

Reviewed `scripts/build_bulk_ablation.py` SHA-256 `2819d37693e0b0c6a0365644ef441aeb2c855b8a18b65576ed31d01a848d5db6` and `tests/test_bulk_ablation.py` SHA-256 `f289c96cd17c823700781695dac552b14705e071152ac2803a96ef91148bb347`. Independently reran all 10 tests successfully in the established cell-eval2 environment. The [protocol](../docs/BULK_ABLATION_PROTOCOL_20261021.md) was written before full-arm generation or metric outcomes.

The isolated builder delegates centering to the original full-source-target implementation. It fixes arithmetic common subtraction/amplitude to 1.0/0.6, keeps CD4 centering and destination-control conversion fixed in the raw-source `common05` contrast, and changes only the final fused bulk scalar in `bulk015`/`bulk060`. Input verification checks source bytes/order/weights, original code dependencies, axes, controls and original incumbent hash. No frozen generator, model or existing scorer was edited.

The diagnostic profiler observes only the exact unchanged emitter function's return frame, copies gene-sized diagnostics, and restores any prior profiler in a `finally` block, including on failure. It does not change arrays or numerical steps. Tests verify bit-identical audited versus ordinary emission, source-population centering, missing-source masks, isolated CD4 behavior, amplitude-effect identity, actual projection capture, existing projection rejection and profiler restoration. The original L1 projection limit 0.03 and fitting limit 0.001 remain enforced by the original kernel.

Independent inspection of all four two-target engineering outputs confirmed 800 by 18,080 counts, identical requested arithmetic profiles/templates/depths/CD4 hashes, and original incumbent count/axis identity on both targets. Full 50-target incumbent identity is still required before old metric reuse. The full scored cohort must pass the existing strict validator; a partial generation cannot be scored. Genuine emitter feasibility failures must remain separate from input/provenance errors and cannot justify relaxing limits.

Fixed requested arithmetic means are not realized means. Engineering variants showed approximately 0.8 CPM maximum per-gene realized drift, explicitly recorded by the builder. This is not the tight-bound dispersion intervention. Gene totals intentionally may change, PDS is not invariant, and all six metrics must be measured. Requested/projected/fitted/realized bulk vectors are retained so projection or rounding cannot be mistaken for the intended intervention.

## Scientific decision rule

Require lower raw normalized capped MSE and lower uncorrected realized log-bulk error, with all five other raw metrics nonworse within 1e-8 and identical finite target cohorts. Uncorrected profiles use full-axis normalization to 50,000 before own-target exclusion. Mean-error and ordinary metric intervals use 2,000 paired target bootstrap draws; no naive per-target bootstrap is applied to the global capped MSE budget. The predeclared source-coverage strata and response/residual alignment are descriptive and cannot add arms or fit amplitudes.

The authorized reserve step requires a complete point-gate pass plus a positive mean-error 95% interval. Select deterministically by lowest raw MSE, then lower uncorrected mean error, then registered arm order; freeze one survivor before one 25-target confirmation. If none qualifies, stop with the reserve unopened. Closed H1 Test/K562 confirmations stay closed; no official submission is authorized by this experiment.

## Outcome status

Full generation and strict validation succeeded for all four arms. The independent full incumbent identity receipt (local artifact omitted: `bulk_baseline_identity.json`) checks all 20,000 rows and 18,080 genes; the new incumbent is also byte-identical to the old file, SHA-256 `03fea48aca0cfad51a08fc12ca5be01b7c7a1430beb7f674a75f0b073d0ba3c6`. A separate read-only audit through the new analysis script independently recomputed all 50 incumbent group profiles and reconstructed its metric aggregates. Six analysis guard tests pass, including full-axis normalization before exclusion, ratio-of-sums rather than mean ratios, finite-cohort changes, incomplete MSE panels and reference-denominator stability.

The [full profile diagnostics](bulk_profile_diagnostics.md) independently confirm all four arms' row depths and saved profile data. Maximum per-gene realized CPM drift is 1.01619 (`common05`), 0.92498 (`bulk015`) and 1.01548 (`bulk060`); the largest bulk projection L1 is 0.000806945, below the unchanged 0.03 limit. No arm is infeasible. These engineering checks do not establish metric nonregression.

## Final scientific verdict: all three candidates FAIL STRICT; retain incumbent

All three candidate scorers completed with exit code 0; the existing incumbent was reused only after full identity verification. The final [analysis](bulk_ablation_results.md) and machine-readable evidence (local artifact omitted: `bulk_ablation_results.json`) completed successfully after independently recomputing final integer-count profiles, checking scored treated counts, source/reference/cache hashes, the complete 50-target MSE component cohorts and the other identical finite metric cohorts (42 for NMAE, 50 for the remaining ordinary metrics). No scoring rerun, frozen-file edit, tolerance change or analysis error repair was needed. A second independent read-only review found no analysis blocker.

| Arm | Realized uncorrected mean error | Raw normalized MSE | Protected-metric failures |
| --- | ---: | ---: | --- |
| Incumbent | 0.002221432735 | 0.970348453 | reference |
| common05 | 0.002220581416 | 0.969871583 | PDS, reach, Jaccard |
| bulk015 | 0.002149668070 | 0.993564953 | PDS, NMAE, fidelity, Jaccard |
| bulk060 | 0.002582763318 | 1.170655702 | PDS, reach, Jaccard |

Restoring half the raw-source common response gives a small true mean-error improvement: gain `8.51319e-7`, paired 95% interval `[2.82821e-7, 1.43909e-6]`, with 33/50 targets improved. It also slightly improves scored MSE, but the three protected losses exceed the preregistered 1e-8 tolerance. This is a tradeoff, not an improvement without affecting other metrics.

Halving bulk amplitude improves uncorrected mean error by about 3.23%: gain `7.17647e-5`, paired interval `[5.01541e-5, 9.20686e-5]`, 47/50 targets improved. Its corrected MSE nevertheless worsens. Candidate-minus-incumbent changes are `ΔU = -7.1764665e-5` and `Δcapped numerator = +4.1384189e-5`, so the effective prediction correction falls by `1.1314885e-4`. That reduced correction is larger than the genuine mean-error reduction. The fixed reference denominator is `0.0017825335324`. Calling the mean biology worse solely from corrected MSE would be incorrect; declaring the arm successful despite MSE and protected losses would also violate the protocol.

Doubling bulk amplitude worsens true mean error by about 16.27%: gain `-0.000361331`, interval `[-0.000419770, -0.000308458]`, only 1/50 targets improved. It also substantially worsens scored MSE. The incumbent's mean response/residual dot is negative (`-0.000112018`), becomes less negative under half amplitude (`-0.0000184137`) and more negative under double amplitude (`-0.000517185`). This is consistent with overextension of the current transferred bulk direction on this development cohort; it does not identify a new optimum amplitude or establish the architecture's cause.

All three candidates fail the registered point gate even though the first two mean-error intervals support positive gains. None is eligible for the conditional reserve step. Keep the 25-target reserve unextracted, the previously closed H1 Test/K562 confirmations closed, and the frozen published model unchanged. No promotion or official submission follows. These results reject the tested scalar/centering recipes under the specified no-regression objective; they are not a general proof against selective mean calibration or better destination-specific response transfer.
