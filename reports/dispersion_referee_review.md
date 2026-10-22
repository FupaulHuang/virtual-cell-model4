# Independent referee review: constrained dispersion experiment

Reviewed 2026-10-21 before dispersion-arm metric outcomes.

**Engineering verdict: PASS for the preregistered development experiment.**
This approves freezing and evaluating the proposed implementation at residual
SD multipliers 1.5 and 2.0. It does not approve a new predictive model, establish
biological realism, certify unchanged metrics, or authorize a submission.

The fixed decision rules are in
`docs/DISPERSION_EXPERIMENT_20261021.md`. The performance verdict is pending the
complete paired six-metric comparison. Small MSE, NMAE or reach tradeoffs must
not be excused after seeing results.

## Reviewed implementation and evidence

- `scripts/constrained_dispersion.py` SHA-256:
  `6d7eebc5c0e4e0ed0418c204aa9e391f680020ab835fa717fb45c647fa160981`.
- `tests/test_constrained_dispersion.py` SHA-256:
  `571d89e0d3210b7a8f3c96cb548746bde1e9eb62e662876616c3594c4343b518`.
- New experiment wrapper inspected at SHA-256
  `acd77c64e9364e1f5e200cf923011d7b1f68767cf44354bbe0ee453f425d3fe4`.
- Independently reran the 12 implementation tests: all passed.
- Ran 24 additional synthetic cases using exact rational arithmetic for every
  gene's CPM mean, independently of the implementation's floating drift cache.
  All 24 changed the matrix and passed exact row/column totals and declared CPM
  bounds; maximum exact drift was 0.841060 of its allowed bound. Aggregate CPM
  variance did not decrease.

The source validates counts before integer conversion, uses int64 arithmetic,
and requires positive bounded row depths. Each pair move has integer bounds
from nonnegativity and the remaining cumulative CPM budget. The zero-sum box
projection and balanced rounding preserve each row depth; opposite updates in
the paired row preserve every raw gene column total. Both depth orderings are
handled. Independent recomputation after passes and at final output checks the
actual count matrix, rather than trusting accumulated increments.

Identity returns unchanged values, dtype and bytes. Zero-total genes remain
zero. Constant inputs, a single cell and infeasible movement are explicitly
reported as unattainable; requested dispersion is not silently substituted for
achieved dispersion. There is no expression jitter: seeds determine pair order
and rounding-tie resolution. Duplicate-cell and density changes are reported.

No blocking defect was found for the validated incumbent and the two specified
multipliers. Two nonblocking general-purpose hardening observations remain:
extremely large finite multipliers can overflow derived squared values, and the
wrapper casts the input block before calling the module. Neither affects the
present hash-locked, independently validated integer incumbent and multipliers
1.5/2.0. They are not grounds to alter running frozen experiment code.

## Independent cohort audit

Reconstructed the selection from the source's observation metadata only, without
reading the expression matrix for this audit. Eligibility, seeded permutation,
all 50 development targets, all 25 reserve targets and every sampled row matched
the saved records exactly. There are 114 eligible targets, 400 rows per
development target and 4,000 unique controls: 24,000 unique source cells and
18,080 genes. No reserved or previously excluded target entered extraction.

- Cohort freeze SHA-256:
  `77a93b85cdd455598ea06894a1597db11539d40169499f3ab231aa968be4666c`.
- Source-row table SHA-256:
  `a16b766712f7cb9ec04de014c80899ecd10d0cddccf947d504cb04e80489e1bb`.
- Incumbent configuration excludes H1 response sources and uses the fixed
  K562/HCT116/HEK293T/CD4 source ensemble.

The baseline scorer initially rejected unused categorical levels inherited from
the full Training file. The repair is documented in
`evaluations/dispersion_h1_training_400/reference_category_fix.json`: only unused
observation categories were removed, actual row identities/labels/gene order
were preserved, and every repaired count was checked against the immutable
selected Training rows. The freeze and row-table hashes above did not change.
The repaired reference SHA-256 is
`6c2f8f2becaa5f6a4fbddd884c62027fe604bc4d171064f6d917c76ea7adbdfd`.
The failed pre-metric run and old file hash are retained. This is an input
metadata repair, not an additional modeling arm or outcome-driven resampling.

## Mathematical and scientific limits

Exact raw margins establish pseudobulk/PDS invariance for the same panel, but
weighted arithmetic CPM is only bounded to 0.01% relative tolerance. The
resulting LFC bound is approximately 1.443e-4 log2 units; actual NMAE still has
to meet the much stricter preregistered metric nonregression test.

The scored expression MSE includes a capped sampling correction, and reach uses
the prediction's confidence ordering. Neither is guaranteed invariant by count
margins or CPM bounds. Wilcoxon is sensitive to distributional ranks and ties,
so changing dispersion need not reduce significance monotonically.

The requested multiplier concerns the square root of the sum of per-gene CPM
sample variances. It does not imply that each gene's SD increases by the same
factor. High-expression genes can dominate the aggregate. Per-gene summaries,
decreased-variance genes, density and duplicate-cell flags must accompany the
global achieved multiplier. Low-count genes may already be maximally dispersed
at their fixed totals.

Even strict development dominance would be evidence on one development cohort.
H1 Training is not previously untouched biological data; the reserve remains
unread and closed H1 Test/K562 confirmation panels remain closed. Any future
promotion needs the preregistered paired uncertainty and separate confirmation.

## Final outcome assessment

**Final verdict: FAIL STRICT for both dispersion arms; promote neither.**
All three metric processes completed with exit 0. Independent checks verified
the complete final matrices, the counts actually passed to the scorer, exact
cell/gene totals, bounded CPM drift, identical scoring/reference identities and
the exact reused real-DE table. Finite target cohorts were identical across all
arms: 50 for PDS/fidelity/reach/Jaccard and 42 for NMAE. Expression MSE is an
aggregate ratio and has no per-target normalized cohort.

The engineering premise was achieved. PDS stayed identical, MSE changed by at
most 1.11e-16, and all four protected metrics passed their preregistered point
nonregression rules. The intended DE improvements were not achieved:

| Raw metric | Incumbent | SD1.5 | SD2.0 |
| --- | ---: | ---: | ---: |
| PDS, higher better | 0.910204082 | 0.910204082 | 0.910204082 |
| MSE, lower better | 0.970348453 | 0.970348453 | 0.970348453 |
| NMAE, lower better | 0.962901741 | 0.962897978 | 0.962899227 |
| Fidelity, higher better | 0.510809091 | 0.500658067 | 0.340200986 |
| Reach, higher better | 0.115049906 | 0.129836861 | 0.152301126 |
| Jaccard, higher better | 0.092176112 | 0.070284180 | 0.030675343 |

Paired 2,000-resample target-bootstrap intervals, using the fixed seed:

| Primary gain versus incumbent | Point gain | 95% interval |
| --- | ---: | --- |
| SD1.5 fidelity | -0.010151 | [-0.022504, +0.000760] |
| SD1.5 Jaccard | -0.021892 | [-0.031890, -0.012517] |
| SD2.0 fidelity | -0.170608 | [-0.226366, -0.116077] |
| SD2.0 Jaccard | -0.061501 | [-0.089242, -0.034811] |

The SD1.5 fidelity decrease is not established by its interval, but its point
estimate fails the positive-gain rule and its Jaccard decrease is clear. Both
primary decreases for SD2.0 have wholly negative intervals. Reach improves in
both arms with positive intervals. NMAE's small point improvements have
intervals crossing zero; they are not evidence of a reliable accuracy gain.

### Why fewer calls did not help overlap

The new development cohort itself—not the earlier H1 Test cohort—has 59,997
reference-significant target-gene pairs under strict adjusted p < 0.05 with the
own target gene removed.

| Arm | Predicted calls | TP | FP | FN |
| --- | ---: | ---: | ---: | ---: |
| Incumbent | 282,674 | 33,512 | 249,162 | 26,485 |
| SD1.5 | 145,309 | 18,786 | 126,523 | 41,211 |
| SD2.0 | 23,817 | 3,922 | 19,895 | 56,075 |

Here TP/FP/FN are relative to this finite experimental call set, not definitive
biological truth. Overcalling was present, but the transformation was not
selective enough: reducing false calls also lost many reference-supported calls.
Pooled call precision increases from 11.86% to 12.93%/16.47%, while pooled call
recall falls from 55.86% to 31.31%/6.54%. These pooled diagnostics are not the
official macro-averaged Jaccard or directional fidelity. They explain why lower
false-positive totals alone cannot establish better overlap; fidelity additionally
uses its direction-match numerator and reference-budget denominator.

Closer marginal variability, if observed in the separate realism audit, does
not ensure correct DE calls when the predicted mean responses remain fixed and
imperfect. This test increases an aggregate dispersion target across all
perturbations; it is not a gene-specific uncertainty model. The result rejects
these two tested policies, not all possible selective or context-dependent
calibration approaches. No additional arm, tolerance, source or target was
selected from these outcomes.

Achieved aggregate SD ranges were 1.48505–1.48835 and 1.98015–1.98923, with no
new duplicate-cell flags. Every group met the requested multiplier within the
declared 1% spread tolerance. Maximum realized LFC drift was approximately
1.443e-4 log2 units. Global MSE was unchanged; in any fixed-bulk experiment a
material MSE change would need to be attributed to its capped sampling
correction rather than improved mean-response biology.

Full paired intervals, win/loss/tie counts, finite target IDs, provenance and
independent matrix audits are preserved in
`reports/dispersion_experiment_results.json` and its Markdown companion. The
25-target reserve remains unextracted, closed confirmation panels were not
reopened, and the published model remains unchanged.

### Analysis provenance correction

The initial analysis guard expected `resolved_de_backend=scanpy` on every run.
Pinned cell-eval2 instead records this field as null when a real-DE table is
supplied, even though predicted DE is computed with explicit scanpy configuration
(`baseline.py:1065–1067`). The analysis failed closed and was corrected only for
that documented provenance representation, after checking configuration, fresh
predicted DE, exact real-table values/fingerprint and preflight/cache hashes.
No scorer was rerun and no metric or acceptance tolerance changed. Original
analysis SHA-256 `5c663db3242fa98228fb30936dc3e99cc16057d9802e0eaa59f7938b7a95cd57`,
its source copy and failure log are retained; the narrow correction is recorded
in `reports/dispersion_analysis_corrections.json`.

## Reviewed freeze synchronization addendum

Engineering approval also applies to the final module SHA-256
`945d102322e9e2688f64bbfae98b1f815622b354086d3f08f6e3805bc79000ea`
and final test SHA-256
`c65cb64f288b9fdfbbadcc0c963927528db1ef9c0fd8ef19048322b5c8060df5`.
I independently reran all 13 tests successfully. Removing only the new
two-line feasibility docstring and two diagnostic fields from the final source
reconstructs the exact previously reviewed module SHA-256 above. Thus the
transformation algorithm did not change. The additional test checks that a zero
CPM drift budget is not relaxed for unequal depths.

Root's hash gate detected the synchronization mismatch. Two prematurely launched
generation attempts were stopped with exit 130 before any metric run; their
partial files and logs are retained under `_aborted_freeze_sync` names. They are
not experimental arms or results. Clean generation must pin the final hash and
must not reuse those partial outputs. No metric outcome informed this addendum.
