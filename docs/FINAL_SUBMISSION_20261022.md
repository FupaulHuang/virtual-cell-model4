# October 22 final submission

The user explicitly requested one final submission of the latest version today. This authorizes a new final D/E/F entry independently of the earlier conditional authorization for a model improvement. Use the latest validated release, `m4_multisource_dual_v1`, with unchanged source banks, source weights, numerical code and hyperparameters. The later experimental emitter has not qualified for promotion.

The final controls, gene list, perturbation list and manifest come from the authenticated production `controls-final` download. Preserve labels D/E/F exactly. The final panel is `vcc2026-test-1`, with 300 targets, 400 cells per target/context, and 18,533 genes. Generate a fresh 360,000-cell prediction; the old validation artifact is not reused.

The frozen source banks cover 141 of the 300 final targets. The remaining 159 use the existing zero-response fallback. This is a material model limitation, disclosed before submission; no new response sources or unvalidated predictor are added to hide it. The final manifest will preserve the exact missing-target list.

Required checks, in order:

1. Verify the official download checksum, control labels/counts/axes, and unchanged frozen model/source hashes.
2. Generate all 900 groups with the existing builder and final inputs. Save actual-count diagnostics and output hash.
3. Run the independent complete prediction validator, then official released VCC CLI 0.3.0 prep with all gene/target/count/density/control checks enabled.
4. Verify package metadata identifies D/E/F and the final panel; preserve package size, SHA-256 and upload checksum.
5. Refresh account eligibility, season, deadline, daily allowance and in-flight state. Create at most one new final entry. If transfer is interrupted, resume that entry instead of creating another.
6. Verify the uploaded object size and checksum, launch scoring, and follow that same entry to a terminal state. A successfully scored final submission has hidden scores; no metric-improvement claim follows.

The global CLI is 0.1.0 and its installation is outside writable roots. The official 0.3.0 wheel is installed under `.tools/vcc-cli-0.3.0`, verified against PyPI SHA-256 `1ebfd26c4fb3e0a58f61b9a15a271b95f662f651291425286c9c2b8db5dc6afb`, and runs with the existing CLI Python environment. Its refreshed skill is stored in `.tools/vcc-skill-0.3.0/vcc`.

Python DNS is unavailable in this environment. If native submission cannot reach production, use the working HTTPS transport to the same official API and its returned Google Cloud Storage upload session, preserving the released CLI's checks and request sequence. No endpoint override or quota bypass is permitted. Credentials and signed session URLs stay in a mode-0700 temporary directory, never in public receipts or logs, and are removed after completion.

The server reports final submissions open, an exclusive deadline of `2026-11-06T00:00:00Z` (last accepted minute November 5, 23:59 UTC), and two submissions per UTC day. Only the last final entry counts for prizes. Today’s work submits once; it does not promise an unseen score or another automatic submission.

References: [official CLI guide](https://vcc-cli-wiki.virtualcellchallenge.org/), final preflight (local artifact omitted: `../reports/final_20261022_preflight.json`), execution state (local artifact omitted: `../reports/final_20261022_state.json`).
