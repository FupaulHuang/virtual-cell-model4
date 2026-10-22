# Bulk protocol and scoring-wrapper referee review

Reviewed 2026-10-21. **PASS for engineering freeze and the registered local experiment, after the two preflight corrections below.** No remaining concrete code blocker was found in the reviewed versions. This is not evidence that all 50 targets are feasible or that any arm improves a metric. Full generation, full incumbent identity, cache creation, terminal statuses, and outcome analysis remain runtime gates.

Scope: read-only review of the [original proposal](../docs/BULK_RESPONSE_REFEREE_PLAN.md), [registered four-arm manifest](../configs/bulk_h1_400_experiment.json), [isolated builder](../scripts/build_bulk_ablation.py), [original builder](../scripts/build_prediction.py), [Atlas model](../external/atlasshift/model.py), [control template](../external/atlasshift/predict.py), [identity verifier](../scripts/verify_bulk_identity.py), and [new scoring preflight](../scripts/run_bulk_cpu_metrics.py). This reviewer did not change code, generate predictions, launch metrics, or submit anything. The parent owns the fixes and executions. Statistical analysis is separately specified in the [bulk protocol](../docs/BULK_ABLATION_PROTOCOL_20261021.md).

## Isolation and numerical units

The four arms implement the intended separate contrasts: incumbent common subtraction 1 / bulk amplitude 0.3; common05 0.5 / 0.3; bulk015 1 / 0.15; bulk060 1 / 0.6. Arithmetic common subtraction stays 1 and arithmetic amplitude stays 0.6. Raw-source common vectors are computed over each source's full available target population, not the 50-target output cohort. Missing target/gene measurements retain masks rather than becoming votes for zero response.

The original builder's shared common-subtraction flag would change both heads and CD4; the new builder avoids that contamination. It calls source-centered raw fusion separately for the arithmetic and bulk heads. CD4 is always centered with common subtraction 1 and `center_scope='source'`. For common05, the centered CD4 effect, availability, weight, destination control, and pre-amplitude conversion remain fixed. For the amplitude arms, the registered scalar intentionally scales the entire fused bulk effect, including CD4. After nonlinear inversion and probability normalization, a separate additive CD4 contribution is not identifiable; “CD4 fixed” describes its inputs and pre-amplitude effect, not a separately fixed component of final counts.

No factor-of-scale or logarithm-base mismatch was found. Raw bulk effects use natural-log differences of `log1p(50000 * probability)`; destination inversion uses that same scale and the bulk control. Arithmetic effects use log2 CPM ratios and the arithmetic control. Controls are computed differently by design: arithmetic averages per-cell normalized expression; bulk pools raw counts. The control template draws the same seeded donors, averages pool-normalized donor profiles, and uses the rounded **mean**, not sum, of donor library depths.

These checks do not remove transfer limitations. K562 raw counts and library weights cover its native 8,246 unique genes; H1 covers its released panel; Orion arithmetic CPM retains original full-library denominators while bulk probabilities use retained-gene count mass. CD4 is publisher LFC data converted through destination controls, not observed source pseudobulk. These are documented approximations, not evidence that the two heads are interchangeable.

The new builder checks original source paths/order/weights and hashes, axes, control hashes, original code hashes, incumbent bytes, and fixed base parameters before generation. Its source inputs match the H1-excluded incumbent. H1 controls condition generation; H1 treated reference rows are used for scoring. No held-out H1 Test or K562 confirmation response enters this experiment.

## Feasibility and what remains invariant

For normalized cell expression `x`, depth ratio `z = depth / mean(depth)`, requested arithmetic probability `p = E[x]`, and bulk probability `b = E[z*x]`, each supported gene requires `min(z) <= b/p <= max(z)`; `p=0` forces `b=0`. Equal depths force `b=p`. These per-gene bounds are necessary, but row-normalization and the simultaneous constraints can still prevent convergence.

The unchanged emitter tightens the depth bounds by 1%, projects the requested bulk profile, rejects projection L1 above 0.03, and rejects fitted moment L1 above 0.001. The new wrapper calls that exact function with its original defaults. Its profiler copies return-frame diagnostic vectors without altering live arrays or replacing the numerical procedure. Failed groups produce an unpromoted failure artifact; thresholds are not relaxed and incomplete arms cannot pass the full scoring shape requirement.

Holding requested arithmetic probabilities fixed does **not** hold realized arithmetic CPM, cell distributions, DE, PDS, or NMAE fixed. Changing bulk requires changing the association between depth and expression. Joint fitting has finite tolerances, and integer rounding plus bulk-total repair can shift realized arithmetic means. The audit records requested/projected/fitted/realized bulk, requested/fitted/realized arithmetic probabilities, raw gene totals, and drift against the original incumbent. It correctly reports rather than enforces the earlier dispersion experiment's unrelated CPM tolerance.

The engineering evidence (local artifact omitted: `bulk_implementation_evidence.json`) reports ten passing tests and strict validation for all four two-target outputs. ACVR1B and ATM incumbent counts match the original exactly. All four requested arithmetic, template, depth, centered-CD4 and CD4-availability hashes match. Maximum engineering projection L1 is 0.00027330025224569913. Maximum per-gene realized CPM drift versus incumbent is 0.840178 for common05, 0.849986 for bulk015 and 0.809282 for bulk060. These two targets do not establish full-cohort feasibility.

## Baseline reuse and scored-input provenance

The preflight reuses the original verified-baseline helper only after checking its SHA against the original frozen DE-cache receipt. It requires the original incumbent terminal exit 0, original implementation hashes, the exact current preset/CPU/Scanpy configuration, current reference and control hashes, unchanged baseline metadata and local-copy hashes, and unchanged exported real-DE bytes. The configuration comparison excludes only output location. Candidate prediction DE is recomputed; only the identical reference DE is supplied.

The new full incumbent must have exactly 20,000 by 18,080 values and the same ordered genes, context labels and target labels as the old incumbent. The identity verifier compares every count in row blocks. The scoring preflight then binds that identity report to the old incumbent hash, new incumbent bytes, its audited receipt, and frozen isolated-builder hash. Thus the two-target engineering check alone cannot authorize baseline reuse.

Candidate checks bind the registered arm, common subtraction, bulk amplitude, base configuration, incumbent identity, builder hash, source paths/hashes/weights and H1-exclusion validation. The scored copy must contain exactly 20,000 treated rows plus the same 4,000 measured controls on the same 18,080 genes. Every treated count is compared to the candidate and every appended control count to the reference. Actual candidate and scoring-copy hashes are recorded before the CPU CLI is called. No anchors or official submission are invoked.

The original cache identifies reference SHA `6c2f8f2becaa5f6a4fbddd884c62027fe604bc4d171064f6d917c76ea7adbdfd`, controls SHA `aea7f3e22fc7b0eae4feca2253d16360c6dafb8cd93da9cf4f4e8547cd044e0c`, incumbent SHA `03fea48aca0cfad51a08fc12ca5be01b7c7a1430beb7f674a75f0b073d0ba3c6`, and exported real-DE SHA `b7b99a2e83df6203758b82887aa7553cbcd2600516747d1450b613d28a98c38b`. These were read from the existing cache (local artifact omitted: `dispersion_real_de_cache.json`), not freshly rehashed by this reviewer.

## Findings resolved before freeze

1. **Blocking treated-block boundary bug — fixed.** The original new preflight sliced the last treated block through row 20,480, including 480 appended controls, while the candidate ended at row 20,000. It now caps both slices at `min(start + 1024, 20000)`. The parent's full scoring-copy guard check (local artifact omitted: `bulk_scoring_copy_guard_test.json`) records PASS on the existing 20,000-treated/4,000-control copy, including this final partial block.
2. **Receipt identity gap — fixed.** The preflight initially recorded the current new-incumbent receipt without comparing it to the receipt hashed during full identity verification. It now requires equality to `identity['new_receipt_sha256']` and requires the identity's builder SHA to match the new freeze. This closes provenance drift between full identity verification and cache creation.

The final outcome report must assess all six raw metrics on identical finite cohorts, with requested and realized profile diagnostics. A lower official normalized MSE must be separated into changed uncorrected mean error and changed capped sampling correction. The fixed-bulk correction-only interpretation from the dispersion study does not apply here. No performance conclusion or candidate selection follows from this review.

Reviewed snapshot SHA-256 values (small files only; Python AST and JSON parse checks passed):

| File | SHA-256 |
|---|---|
| `scripts/build_bulk_ablation.py` | `2819d37693e0b0c6a0365644ef441aeb2c855b8a18b65576ed31d01a848d5db6` |
| `scripts/run_bulk_cpu_metrics.py` | `92a483194b29e5a67d81a1c12dc3fc9ce4090297e1fe860596f49745104e52f6` |
| `scripts/verify_bulk_identity.py` | `7238ec5d32a629f6f19d08d3b6ce015f73a210e1a5f8dc3e56720802f4990af2` |
| `configs/bulk_h1_400_experiment.json` | `972476b7d36d5d90cc069b7b2862c8878b65905892d1a1258fc808550d3276fa` |

## Bounded final-summary review

Independent read-only review of `scripts/summarize_bulk_ablation.py` and its six guard tests: **PASS**, with no concrete blocker. Analysis SHA-256 reviewed: `3e95ef6edfbc6e5bef0a8bbe1d9504c9c2d6cadc21c8fac4211e863a94494ce4`. The mean-error endpoint normalizes whole-axis pooled counts before own-target exclusion; its per-target denominator matches the official scored-gene count. The correction decomposition is `delta effective prediction correction = delta uncorrected error - delta capped numerator` under the verified fixed reference. Complete finite MSE component cohorts and the ratio of component sums are checked, with no naive per-target bootstrap of global MSE. Protected directions/tolerance, paired mean-error uncertainty, and deterministic MSE-then-mean-error-then-arm-order reserve selection match the frozen protocol. No criteria, code, or predictions were changed by this review, and no metrics were run by this reviewer. Parent subsequently confirmed all three metric jobs exited zero; independent final analysis owns outcome interpretation.
