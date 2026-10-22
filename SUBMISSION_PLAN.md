# Model 4 generation and submission plan

Updated: 2026-10-21. The user authorized the October 20 submission; no additional
confirmation was required for this challenge submission.

Completed: entry `[PRIVATE_IDENTIFIER_REMOVED]` is **published**, with official score
**0.1470186853** and rank **565** at approximately 23:54 HKT. All gates below
were satisfied. The unchanged artifact must not be submitted again. See
[the final report](reports/FINAL_REPORT.md) and terminal receipt (local artifact omitted: `submissions/m4_multisource_dual_v1_status.json`).

## Target and quota

Submit the selected **2026 validation** candidate on October 20, Asia/Hong_Kong,
subject to source acquisition, validation and service availability. Final test
contexts are released October 22 and are not available for the October 20 artifact.

The challenge permits two scoring submissions per team per UTC day and one
submission in flight. Inspect the current production state before starting.
Do not consume a second slot on an unchanged artifact. Use a second slot only
if the first result and local evidence justify a materially different candidate.

## Gate S1 — freeze the candidate

- [x] Complete `reports/MODEL_SELECTION.md` using the experiment plan.
- [x] Save a machine-readable configuration and hashes of every source bank.
- [x] Confirm source coverage across all 300 target identities and document any gaps.
- [x] Record approximated quantities and correlated-source mixing explicitly.
- [x] Verify the reserved confirmation result is consistent with the stated limitations.

Why this matters: without a frozen configuration and provenance, the scored
submission cannot be attributed to a reproducible model.

Scheduling note: full generation and packaging may proceed on an explicitly
provisional, frozen candidate configuration alongside evaluation. Its artifacts
remain unselected until S1 is satisfied. Final submission never precedes model
selection, confirmation review or the required validation gates.

## Gate S2 — generate the full prediction

Run `scripts/build_prediction.py` with the selected banks and exact official
gene/target CSVs, all three context control files, 400 cells per target, the
selected pooling/amplitudes and a recorded seed. Use bounded worker concurrency
and one BLAS thread per worker.

The output must contain:

| Property | Required value |
| --- | --- |
| Shape | 360,000 cells × 18,533 genes |
| Contexts | A, B, C; 120,000 cells each |
| Targets | Exact official 300 in every context |
| Groups | 900 context-target pairs |
| Cells per group | Exactly 400 |
| Gene order | Exact `data/gene_names.csv` order |
| Matrix | Finite, nonnegative, integral raw counts |
| Control cells | None |
| Per-cell counts | At most 1,000,000 |
| Stored sparse entries | At most 4,750,000,000 |

- [x] Generation finishes without missing blocks or moment-fitting errors.
- [x] Final receipt records dimensions, group diagnostics, realized moments,
      run arguments, elapsed time, file size and SHA-256.
- [x] Inspect realized expression errors and density across all groups, not
      only a smoke subset.

Why this matters: the earlier successful smoke test proves code execution on
3 targets; it does not prove that the complete challenge panel is valid.

## Gate S3 — independent validation and official packaging

- [x] Run `scripts/validate_prediction_receipt.py` against the complete H5AD.
- [x] Verify every group and count block; reject unexpected labels, explicit
      sparse zeros, fractional counts, nonfinite values or cap violations.
- [x] Run production `vcc prep` with official gene and perturbation lists.
- [x] Preserve package/prep logs and hash the resulting `.vcc`.

Recorded strict packaging arguments are in
`submissions/m4_fullsource_candidate.manifest.json`. For a future run, replace
`SELECTED_RUN` below and retain all flags:

```bash
.venv/bin/python scripts/validate_prediction_receipt.py \
  outputs/SELECTED_RUN/prediction.h5ad

vcc prep outputs/SELECTED_RUN/prediction.h5ad \
  -g data/gene_names.csv --perts data/pert_counts.csv \
  -o submissions/SELECTED_RUN.vcc \
  --contexts A,B,C --cells-per-pert 400 --encoding 32 \
  --require-counts --reject-controls --verify-targets --check-cell-counts \
  --expected-gene-dim 18533 --max-cell-dim 400000 \
  --max-nnz 4750000000 --max-counts-per-cell 1000000 --json
```

Keep all official target, count and raw-count checks enabled. The package
validator's success is format evidence; it does not establish biological quality.

The expected roughly 2.45 billion stored entries require 64-bit CSR row pointers,
which both our writer and the installed VCC codec support. Default encoding 32
controls expression-value precision, not pointer width. Set `TMPDIR` to a
workspace directory with at least 60 GiB free for the uncompressed temporary
H5AD and compressed package; the shared `/tmp` volume has less spare capacity.

## Gate S4 — submit and follow the existing entry

- [x] Confirm production authentication, quota and absence of another active entry.
- [x] Submit with a descriptive model name and concise source/method description.
- [x] Save the returned entry ID and upload verification receipt immediately.
- [x] Upload completed in one verified transfer; no interrupted entry or duplicate remains.
- [x] Poll that same entry until published or failed; an observation timeout
      alone does not mean the scoring job failed.
- [x] Save the complete terminal status and all six raw/scaled metrics.

```bash
vcc submit submissions/SELECTED_RUN.vcc -m "model4-selected-name"
vcc status RETURNED_ENTRY_ID --json
```

The installed CLI already has a production login. Credentials stay in its
existing credential store and must not be printed or copied into this workspace.

## Gate S5 — report the measured outcome

- [x] Compare the official result with Model3's **−0.2946503069** using the
      same validation panel/anchor version.
- [x] State the score difference, rank at retrieval time and component changes.
- [x] Distinguish local holdout scores from official validation scores.
- [x] Link the summary, model selection, configuration, artifact hashes and entry receipt.
- [x] Update `README.md`, `PLAN.md` and the continuation state with completed work.

If a result is worse, report that result directly and diagnose which metrics
changed. Do not label it an improvement based on a favorable subset of metrics.
If the service has changed its panel or anchor version, disclose that difference
and avoid a same-scale improvement claim until comparable scores are available.
If scoring fails, retain the failure response and correct the specific defect
before another authorized attempt.

## Completion evidence

The requested work is complete only when the prior-work summary, reproducible
selected model, full valid prediction and actual submitted scoring outcome have
all been checked. If the service or data remains unavailable, preserve a concrete
continuation record and keep the goal active rather than claiming submission.

## Completed gate evidence

| Gate | Evidence |
| --- | --- |
| S1 | [Frozen decision](reports/MODEL_SELECTION.md), current release (local artifact omitted: `models/m4_multisource_dual/release.json`), [H1 Test](reports/h1_test_confirmation.md), [K562](reports/k562_confirmation.md). |
| S2 | Prediction receipt (local artifact omitted: `outputs/m4_fullsource_candidate/prediction.receipt.json`): 360,000 × 18,533, 900 groups, 2,421,598,713 stored entries. |
| S3 | Full validation (local artifact omitted: `outputs/m4_fullsource_candidate/validation.json`), 64-bit pointer audit (local artifact omitted: `outputs/m4_fullsource_candidate/csr_pointer_audit.json`), strict prep (local artifact omitted: `reports/m4_fullsource_candidate_prep.json`), package manifest (local artifact omitted: `submissions/m4_fullsource_candidate.manifest.json`). |
| S4 | Create (local artifact omitted: `submissions/m4_multisource_dual_v1_entry.json`), verified upload (local artifact omitted: `submissions/m4_multisource_dual_v1_upload.json`), launch (local artifact omitted: `submissions/m4_multisource_dual_v1_launch.json`), published status (local artifact omitted: `submissions/m4_multisource_dual_v1_status.json`). |
| S5 | [Final results](reports/FINAL_REPORT.md), [official comparison](reports/official_result_review.md), updated README, plans and continuation state. |

The native submission path could not use its default lock directory, and its
writable-runtime invocation encountered DNS failure. The accepted upload used
permitted HTTPS requests to the official API and returned storage session,
retaining authentication, quota, lock, size/checksum and launch checks. No
endpoint override, quota bypass or escalation was used. Temporary private
transport files were removed. Transport record (local artifact omitted: `reports/submission_transport.json`).

These are validation A/B/C instructions. A later final D/E/F submission needs
its own current CLI, official manifest, controls and target panel; see
[the next iteration plan](docs/NEXT_ITERATION_PLAN.md).
