# Model4: completed implementation and official result

**m4_multisource_dual_v1 scored 0.1470186853 in the official 2026 validation
challenge**, improving on Model3's −0.2946503069 by **0.4416689922**. Both entries
use `vcc2026-val-1` and the identical `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`
anchor version. The Model4 status is terminal `published`, with no error.

Entry `[PRIVATE_IDENTIFIER_REMOVED]` was submitted on **October 20, 2026 at 23:30:54 HKT**.
Rank was **565** when checked at approximately 23:54 HKT. The saved leaderboard
snapshot from 23:13 HKT has a leader score of 0.4733651969; the new model improves
substantially on our earlier work but remains below leading performance.
Rank and leaderboard comparisons are timestamped observations.

Evidence: official result (local artifact omitted: `../submissions/m4_multisource_dual_v1_status.json`),
[independent result comparison](official_result_review.md),
leaderboard snapshot summary (local artifact omitted: `leaderboard_late_summary.json`).

## What changed and why

Earlier models had sparse same-target coverage and converted differences of
mean-log expression into raw-count multipliers. Those quantities differ from
arithmetic expression fold changes, and multiplying existing counts cannot
activate genes at donor-zero entries. Local vector-space proxies also failed
to predict actual challenge scores. [Prior-work analysis](../docs/RESEARCH_SUMMARY.md)
and [count-emission review](../docs/EMISSION_REVIEW.md) explain the evidence.

The selected model combines measured responses from public K562 GWPS (weight 2),
H1 2025 (2), X-Atlas Orion HCT116 (1), HEK293T (1), and CD4 DE summaries (0.5).
All 300 official targets have at least 20 measured source cells in one eligible
source. This means target coverage, not complete measurement of every gene in
every source. Missing target-gene measurements remain masked.

Separate arithmetic-CPM and pseudobulk response heads use amplitudes 0.6 and
0.3, prior strength 100,000 counts, source-common subtraction 1 and clipping 3.
Four destination-control donors are pooled per generated cell. Moment fitting
and seeded integer emission preserve cell library totals. No agreement filter
or promoter prior is used. CD4's bulk head is an explicit approximation because
the published DE summaries do not provide raw library-weighted pseudobulk.

AtlasShift's pinned MIT source provides the numerical kernels. Model4 adds
audited source preparation, experiment orchestration, count/axis validation and
provenance. [Upstream sources](upstream_sources.json) and
environment record (local artifact omitted: `../models/m4_multisource_dual/environment.json`) preserve attribution.

## Selection and confirmation

Twelve candidates were compared on the same public H1 Validation reference,
CPU/scanpy metric backend and fresh five-split reference anchors. All H1 response
sources were excluded. The selected model's **local score was 0.2427562224**;
the runner-up, Orion plus CD4, scored 0.2249690328. The initial filtered-PIE
dual-profile model scored 0.105982 and the neutral emitter −0.023070.
These local scores use different data and anchors from the official score.

The original one-split selection was frozen before reserved confirmation
outcomes. The independent five-split comparison confirmed the same winner.
H1 is added as an external response source for production A/B/C; it remains
excluded from both H1 local evaluations.

| Confirmation | Response family excluded | Outcome versus fixed neutral emitter | Important scope limit |
| --- | --- | --- | --- |
| H1 Test | All H1 responses | All six aggregate raw metrics improve | 100 targets; NMAE uses 72 and reach 97; same context as selection, new target outcomes. |
| K562 | All K562/GWPS responses | All six aggregate raw metrics improve | 50 targets; 8,246 genes; NMAE uses 12 and reach 34; pooled controls without batch-adjusted DE. |

Both confirmations used 200 predicted cells per target; production uses 400.
No settings were retuned from confirmation results. The neutral emitter is a
control-moment comparator, not the official generic-response baseline. No
cross-panel averaging or reuse of H1 Validation anchors was performed.

Evidence: [selection](MODEL_SELECTION.md),
[all candidates](h1_all_candidate_comparison.md),
[H1 Test confirmation](h1_test_confirmation.md),
[K562 confirmation](k562_confirmation.md).

## Official components and remaining limits

All reference-scaled components below are higher-is-better. Raw MSE and NMAE
are lower-is-better; the other four raw metrics are higher-is-better.

| Component | Model4 raw | Model3 scaled | Model4 scaled |
| --- | ---: | ---: | ---: |
| Perturbation discrimination | 0.791568 | 0.018404 | 0.645588 |
| Expression MSE | 0.952085 | 0.000000 | 0.037668 |
| LFC NMAE | 0.917743 | −0.000361 | 0.138194 |
| Direction fidelity | 0.507158 | −1.697095 | −0.018683 |
| Direction reach | 0.149042 | −0.007310 | 0.079313 |
| Significant-set Jaccard | 0.030320 | −0.081540 | 0.000032 |
| **Overall: mean of six scaled components** | — | **−0.294650** | **0.147019** |

All six raw and scaled components improved over Model3. Most of the aggregate
gain comes from recovering its severe fidelity failure and improving
discrimination. The remaining gene-level accuracy is weak: fidelity is below
the reference baseline, Jaccard is essentially at it, and expression MSE gains
are small. No official per-context breakdown was returned.

H1 Test also showed 206,159 predicted significant target-gene pairs versus
115,411 in the reference. This motivates testing effect and variance calibration;
it does not establish the cause of hidden-context errors. Source normalization
and gene panels differ, and broad coverage alone does not guarantee reliable
gene-level transfer. Five anchor splits calibrate the score, not uncertainty
across independent model runs. Multiple changes were bundled, so the official
improvement cannot be assigned causally to one change.

## Artifact and submission verification

The prediction contains **360,000 cells × 18,533 ordered genes**, 900 groups of
400 cells, and 2,421,598,713 stored entries. Every count is finite, nonnegative
and integral; no control cells are included. Maximum library size is 51,786.
Independent checks covered the entire matrix and 64-bit sparse row pointers.
Strict production `vcc prep` accepted the artifact without warnings, drops or
gene reordering and preserved counts.

The package is **3,967,170,560 bytes**. Uploaded size and storage MD5 were checked
against the local package before launching scoring. Both MD5 values were
`XET//2QKJCPhd8oY2oWvfQ==`. Source/count/axis tests and prior exhaustive artifact
checks passed; the final audit checks receipt consistency without rerunning
expensive generation or scoring.

| Artifact | SHA-256 |
| --- | --- |
| Frozen model definition | `44d70efadcd12f7ba190f2d6935bceee3553719702cbe7e13673d2a2fc250993` |
| Production configuration | `207175a3586028582b5782a31c61118528daa8cd11748ff3cdc053881d7efb59` |
| Prediction H5AD | `57e5a06aa99067ef5a2cdce2211d371a466753694877968e8cb5dba2fbb87d42` |
| Submitted VCC package | `d7c4d014060cb54a176c22adb76ea28b70af5462e97147745d28de1e262dfaae` |

Evidence: full validation (local artifact omitted: `../outputs/m4_fullsource_candidate/validation.json`),
CSR audit (local artifact omitted: `../outputs/m4_fullsource_candidate/csr_pointer_audit.json`),
strict prep (local artifact omitted: `m4_fullsource_candidate_prep.json`),
package manifest (local artifact omitted: `../submissions/m4_fullsource_candidate.manifest.json`),
verified upload (local artifact omitted: `../submissions/m4_multisource_dual_v1_upload.json`),
launch (local artifact omitted: `../submissions/m4_multisource_dual_v1_launch.json`),
transport details (local artifact omitted: `submission_transport.json`).

## Reproducing and extending the work

The original [production config](../configs/abc_fullsource_candidate.json) and
frozen model (local artifact omitted: `../models/m4_multisource_dual/model.json`) retain their historical
pre-confirmation labels because their exact bytes were hashed before confirmation.
release.json (local artifact omitted: `../models/m4_multisource_dual/release.json`) is the current published
state; those old labels do not mean the model remains unselected.

From this workspace, inspect the exact inference command without rewriting the
submitted artifact:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_prediction_config.py \
  configs/abc_fullsource_candidate.json --dry-run
```

For a new reproduction, copy the config to a new run name and change only its
`name` and `output` fields before executing without `--dry-run`. The original
output snapshot prevents an accidental overwrite. Record the new config hash,
then independently validate the completed H5AD and package using all strict
flags in [SUBMISSION_PLAN.md](../SUBMISSION_PLAN.md). Reproduction is expensive;
the existing artifact is already verified and published.

Source preparation is documented by the individual bank provenance, the
GWPS acquisition note (local artifact omitted: `../data/raw/gwps_subset/ACQUISITION.md`),
[GWPS subset builder](../scripts/gwps_subset_bank.py),
[CD4 ranged-file builder](../scripts/build_cd4_ranged_bank.py),
[Orion raw builder](../scripts/build_orion_raw_bank.py), and
[H1 builder](../scripts/build_h1_bank.py). GWPS uses verified byte ranges for
the selected target/control subset; it is not a claimed full 65.8-GB download.
The complete alternate CD4 file was verified against its published hash.

For selection reproduction, use the per-arm configurations and receipt chains
in the comparison JSON (local artifact omitted: `h1_all_candidate_comparison.json`) and the frozen
`evaluations/h1_validation/` reference. Source exclusions are part of the
protocol, not optional flags. The H1 Test and K562 results are now observed;
they cannot be reused as fresh confirmation for a later model.

The requested plans are [implementation](../PLAN.md),
[experiments](../EXPERIMENT_PLAN.md), and [submission](../SUBMISSION_PLAN.md).
The [next iteration plan](../docs/NEXT_ITERATION_PLAN.md) specifies new work,
why it matters and the evidence needed before another candidate is promoted.
No additional submission is required to complete this October 20 run.
