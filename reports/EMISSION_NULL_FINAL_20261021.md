# A concrete next direction after the unsuccessful ablations

**The zero-response generator fails the null sanity check on this H1 panel.**
It reports thousands of significant genes even when no perturbation response is
supplied. This gives a concrete reason to redesign the cell emitter before more
global amplitude/dispersion searches. It is a diagnostic advance, not a new
leaderboard gain.

## Completed controlled comparison

The reference uses only public H1 non-targeting cells: ten groups of 400 genuine
controls compared with a separate 4,000-cell control set. Both sets are disjoint
from the 4,000 inference controls. Batch and guide totals match exactly. Exact
joint matching was infeasible in one rare stratum, so a metadata-only amendment
minimized joint reallocation: 14/4,000 cells (0.35%). No new expression rows had
been read for that failed matching attempt. Real null and comparison cells have
the same recorded joint allocation.

Ten existing gene names are technical aliases only; the measured null cells have
no gene knockdown. The official vcc2026 DE settings, CPU/scanpy, real comparison
controls, per-target BH and strict adjusted p < 0.05 are held fixed. Only DE-count
diagnostics are run; no challenge performance score is computed on null labels.

| Arm | Median significant genes per group | Range |
| --- | ---: | --- |
| Independent measured controls | **0** | 0–0 |
| Raw donor replay | **0** | 0–0 |
| Four-donor pooling, without moment fitting | **4,672** | 4,610–4,700 |
| Unchanged dual emitter, zero response | **4,323.5** | 4,244–4,378 |
| Existing response probe | 5,603 | 5,260–5,914 |

The neutral generator exceeds donor replay in all ten aliases and calls roughly
40% of the tested genes. It meets the preregistered engineering flag. Since the
empirical null medians are zero, the registered ratio uses a one-call floor; it
must not be interpreted as an ordinary biological fold change relative to zero.

Pooling reduces the variance-trace ratio to approximately 0.285 and increases
detection-distribution error from about 0.012 for measured/replayed controls to
about 0.082. Matching mean profiles cannot undo this distributional mismatch.
The rank-based DE test correctly detects the synthetic distribution change;
interpreting it as a predicted perturbation effect is the modeling problem.

## Relation to the full response model

Under the same comparison controls, 29,933 of 55,696 existing response-probe DE
calls (**53.7%**) are also called by the generator with its response set to zero.
Pooling without fitting shares 55.2%. These are matched diagnostic call-set
overlaps. They do not prove that those individual genes are false against real
perturbation outcomes, and the probe is not evaluated against the null groups
as if they were perturbed-cell truth.

This clarifies why scalar adjustments were insufficient: a substantial call
background can originate in emission itself. However, the earlier dispersion
experiment also showed that suppressing calls can lose genuine responses.
Fixing null behavior alone is therefore not sufficient for a better challenge
model; gene-specific response transfer still needs improvement.

## What should be done next

1. Build an emitter whose zero-effect output preserves selected measured control
   cells exactly, without donor averaging or forced gene totals at that identity
   point. Verify its null behavior on H1 and the supplied A/B/C controls.
2. Fit a small context-aware response-confidence/magnitude adapter using public
   source contexts and strict destination-response exclusion.
3. Test old/new response models and old/new emitters in a frozen four-arm ablation
   on actual perturbation data. Keep all six metrics and mean error visible.

The detailed [next-model plan](../docs/NEXT_MODEL_REDESIGN_20261021.md) defines
implementation contracts, isolation tests, data boundaries and selection rules.
It has not yet been implemented or scored. The existing published model remains
unchanged, and the reserved perturbation responses remain unextracted.

## Verification and limits

Seven engineering tests and a small official-DE pipeline smoke test passed.
All five complete prediction matrices passed count/axis/hash validation. A separate
source check compared every one of the 8,000 selected reference rows with the
original counts and verified population disjointness and batch/guide margins.
All five diagnostic runs exited 0; their configurations, tested gene sets and
reference-DE tables match. Scored copies were checked against actual prediction
and comparison-control counts before execution.

Review agents were unavailable because of usage limits. Root performed code
inspection, meaningful tests and complete data/provenance checks; this is not
claimed to be an external peer review. The ten aliases share controls and donor
pools and are correlated technical diagnostics. Donor replay retains original
donor depths, while pooled methods use pooled depths, so the comparison concerns
complete emission policies. The magnitude of the finding is established here
for H1, not for every possible context.

Evidence: [protocol](../docs/EMISSION_NULL_PROTOCOL_20261021.md),
[full results](emission_null_results.md), JSON diagnostics (local artifact omitted: `emission_null_results.json`),
source audit (local artifact omitted: `emission_null_source_audit.json`),
terminal statuses (local artifact omitted: `emission_null_terminal_statuses.json`).
