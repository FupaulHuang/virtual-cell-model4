# Candidate generation for independent injection–recovery

This builder consumes the independently frozen truth manifest (local artifact omitted: `../evaluations/injection_recovery_400/manifest.json`) defined by the [truth protocol](INJECTION_RECOVERY_PROTOCOL_20261021.md). It does not generate reference cells, fit a biological response model, inspect actual H1 perturbed responses, access the reserve, run official six-metric scores, or submit anything.

## Fixed methods, inputs and stages

For each composition/Dirichlet–multinomial family and seed 20261026/20261027, use that seed's 4,000 inference controls and the family's eight frozen aliases. These are two seeds per family, not a cross-product of truth and emission seeds. Current and preserving/uniform are mandatory. Preserving/balanced is included only if `reports/balanced_sampling_results.json` records both the combined gate and injection eligibility as true. Freeze that eligibility decision before opening candidate outcomes; otherwise the balanced arm is explicitly excluded.

Stages are `freeze`, `prepare`, `generate --family --seed`, and `finalize`. Freeze builder/tests/this protocol, unchanged numerical dependencies, truth inputs and the independent balanced-gate receipt. Prepare all four family/seed request banks and any balancing features before generating any candidate. Generation checks their hashes. Source role/gene/alias axes and input hashes must match the truth manifest; inference files must contain only non-targeting cells. Comparison/reference files are bound by provenance but their expression outcomes do not determine requests or seeds.

## Identical anchored requests

Compute one canonical inference arithmetic and bulk baseline with the frozen preserving emitter's `population_profiles(raw)` on the complete ordered gene axis. Retain this baseline separately from the oracle's own unperturbed baseline. For each alias/head, let `delta = oracle - oracle_baseline` and form `before_clip = inference_baseline + delta`, then `after_clip = max(before_clip,0)` and `request = after_clip / sum(after_clip)`.

If both oracle deltas are exactly zero, return **copies of the inference baselines without normalization**. A NULL flag with nonzero delta is an input error, not permission to erase a signal. Verify the declared delta arrays agree with the oracle/baseline arrays. Every eligible method receives exactly the same saved requests. Current-template reduction baselines are independently recorded and must agree with the canonical baseline within 1e-12; they do not replace requests by a method-specific target.

For every gene preserve absolute oracle, oracle baseline, inference baseline, oracle delta, before/after clipping, removed negative mass, normalized request, requested delta, closure/clipping delta error, and request-minus-oracle offsets. Clipping and closure can change signs; this is an audited request transformation, not exact LFC transfer.

## Sampling and emission seeds

Use the frozen Atlas seed function with the family and alias strings:

- Template seed: `seed_for('injection:template:' + family, seed)`. One frozen pool-four template per family/seed is shared across aliases.
- Donor seed: `seed_for('injection:donor:' + family + ':' + alias, seed)`. Uniform and balanced use the same integer seed, with their respective fixed sampling algorithms; they need not select the same cells.
- Transport/rounding seed: `seed_for('injection:transport:' + family + ':' + alias, seed)`, identical across methods. This fixes the RNG input, not equivalence between different numerical algorithms.

Current calls unchanged `control_template(..., cells=400, pool_k=4)` and `emit_audited`, which invokes the frozen dual-moment count kernel. Record its actual projected bulk and fitted 400-cell-template moments; full-population fit is **not applicable**.

The current method's selection log records all 1,600 source-cell indices/IDs, grouped four per template row, for every attempted alias. Its donor-baseline columns describe the unfitted pooled template, not 400 actual raw cells. Preserving logs describe individual selected raw donors. Audits mark which seed domains were actually used, so the alias-specific donor seed is not misattributed to the current method's context-shared template.

Preserving methods call unchanged `fit_population` once per alias on all 4,000 raw inference cells, then sample 400 actual donors and call unchanged `emit_counts`. Uniform samples without replacement. Balanced uses the independently approved cube/pivotal sampler. It fits the same fixed four-PC control-only features with seed20261026 for each distinct inference count population; identical populations may reuse a byte-verified feature matrix. Freeze features and reconstruction state before candidate generation. No donor averaging, empirical gene-total repair, selected-seed optimization or parameter change is introduced.

Keep original projection L1<=0.03 and fit L1<=0.001 limits and all other numerical defaults. The preserving NULL fit must be identity, and every output count must equal the selected donor count. This is a zero-delta identity contract; no identity claim is made for absolute oracle profiles from a different finite baseline.

## Artifacts, failures and downstream contract

Per-arm output: `outputs/injection_recovery_400/{family}/{method}_s{seed}/`. Save `manifest.json`, `prediction.receipt.json`, `alias_audit.json`, and `gene_profiles.parquet`. A complete arm has `prediction.h5ad`, 400 cells for every one of eight aliases, independently validated counts/axes/hash, and `complete=true`. Per-gene stages include `oracle_p/b`, `oracle_baseline_p/b`, `inference_p/b`, `oracle_delta_p/b`, `anchored_before_clip_p/b`, `clipped_negative_p/b`, `post_clip_p/b`, `request_p/b`, `request_delta_p/b`, `request_delta_error_p/b`, `request_minus_oracle_p/b`, `population_fitted_p/b`, `sample_expected_p/b`, `bulk_projected_b`, `realized_p/b`, `realized_request_drift_p/b`, and `realized_count_sums`.

A read-only return-frame observer records the preserving kernel's projected bulk; it neither replaces nor alters fitting. Conditional sample expectations use the selected donors' actual depths. Raw donor baselines and their counts/indices are recorded. Full-population fit columns are NaN for current; unavailable realized stages are NaN for failed aliases, with explicit status.

Attempt all prespecified aliases and retain each failure. Recognized original feasibility errors are distinguished from engineering failures. No failure relaxes a limit, changes a seed, adds a method, or drops an alias from an apparent complete result. Any missing alias leaves a closed `prediction.partial.h5ad`, `complete=false`, and no scoreable prediction path. Incomplete methods are unscored, and no favorable subset-performance claim is permitted. Final `outputs/injection_recovery_400/manifest.json` lists every planned arm and completeness. The independent downstream runner owns DE and recovery interpretation, including reference-NULL and held-out-template/seed gates.
