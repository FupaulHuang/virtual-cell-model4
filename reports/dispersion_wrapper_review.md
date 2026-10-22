**PASS for the fixed experiment when candidate metrics use the verified entrypoint.** Reviewed at `2026-10-21T02:46:54.338124+00:00`. No metric process was launched or restarted, and no generation/scoring/model file was modified by this review.

Reviewed code: `scripts/apply_dispersion_experiment.py`, the frozen `scripts/run_dispersion_cpu_metrics.sh`, and the new `scripts/run_verified_dispersion_metrics.py`. The module algorithm/tests were reviewed separately; this review focused on integration and provenance.

The two original shell-wrapper gaps are addressed by the new preflight:

- Existing `local_scoring_input.h5ad` copies cannot be silently reused for different prediction/reference inputs: the preflight compares receipt paths, current input hashes, expected dimensions, gene order, treated labels and the appended 4,000 control labels, then records the local-copy hash.
- Reused real DE requires a successful incumbent terminal exit, completed metric artifacts, freshly computed real DE, the exact reference/control hashes and CPU/scanpy preset/configuration. A separate immutable cache receipt pins the DE table and baseline identities; candidates must match it before invoking the frozen shell wrapper.

The bare shell script still contains its original existence-only shortcuts. Candidate runs must therefore enter through `run_verified_dispersion_metrics.py sd15` or `sd20`, after `freeze-real` succeeds. Calling the shell directly bypasses the approved provenance guards.

The final preflight also checks each candidate's recorded incumbent hash, frozen module/wrapper hashes, arm-specific requested multiplier, seed, CPM tolerances, maximum passes and pair window. The actual sd15/sd20 receipts record 1.5/2.0 respectively, seed 20261021, relative tolerance 1e-4, absolute tolerance 1e-8, eight passes and window four. Their recorded generation-code hashes match the engineering freeze. All files listed in that freeze remained unchanged during this review.

`EvalConfig` comparison is correct for the installed scorer: `to_dict()` returns `dataclasses.asdict()`. Reconstructing the intended vcc2026 CPU/scanpy configuration and comparing it with an existing run from the same CLI/preset, with only thread count adjusted from two to four and output directory removed, produced no differences. The new incumbent's `run_params.yaml` was not yet present, so its exact runtime comparison remains a required preflight gate after terminal success.

Observed data and receipt checks passed:

- Reference: 24,000 × 18,080, comprising 50 target groups with exactly 400 genuine cells each and 4,000 controls.
- Incumbent prediction: 20,000 × 18,080, with 400 cells per target; independent validation passed with H1 response-source exclusion.
- Reference, controls and prediction gene axes match. Control observation IDs match the reference control block; context-control and current scoring-copy receipts match the audited manifest and incumbent hashes.
- Current incumbent target order is sorted and observation IDs match Writer's ordering, so the wrapper's sorted group traversal preserves the actual incumbent cell/gene indexing.
- Development targets do not overlap previous H1 Validation, H1 Test, K562 confirmation, or the 25 reserved targets. The cohort remains development data in an already studied H1 context, as its manifest states.
- The derivative wrapper verifies integer/nonnegative output, exact cell depths, exact gene totals and bounded arithmetic-CPM drift for each group. It retains source/configuration provenance and records the explicit derivation and final output hash.
- Derived moment keys correctly describe drift **to the incumbent**, rather than pretending to remain errors against the original requested profiles. The experiment-specific summarizer must use these keys.

Real-DE reuse intentionally changes `de_real_fingerprint` from null in the computing incumbent to a supplied-table fingerprint in candidate runs. This difference is expected only after the frozen real-DE identity checks; supplied prediction-side DE is not part of the approved path.

Limits: this is an engineering/provenance approval, not evidence of a favorable dispersion or metric result. It does not rerun module tests or independently rescan every candidate count matrix. First-time local-copy provenance relies on the trusted copy construction and its source receipt; subsequent identity is hash-pinned. The independent result summarizer is responsible for final realized-margin/CPM/LFC comparisons and finite-cohort checks. No other concrete wrapper blocker was found for the inspected fixed inputs.

Reviewed preflight SHA-256: `ef2f3b94d37e458dfd56a81e4d75ce475004d226ef6a2f7da4493e3f3d9936d4`.
Frozen generation wrapper SHA-256: `acd77c64e9364e1f5e200cf923011d7b1f68767cf44354bbe0ee453f425d3fe4`.
Frozen metric shell SHA-256: `f0eb41a002e8d69c32d885703ba612019dfb887faa490899362000a6afffde9d`.
