# Zero-response emission diagnostic — 2026-10-21

The user authorized proceeding after two unsuccessful ablations. This diagnostic
asks whether the generator itself creates apparent DE without a perturbation
response. It is not a candidate selection benchmark and produces no challenge
score. Existing models, metric rules and past acceptance decisions remain intact.

## Fixed inputs and split

Use only the public H1 Training non-targeting cells for new reference data.
Retain the exact 4,000 inference controls and gene axis used in the completed
400-cell development study. Select 4,000 separate comparison controls and ten
disjoint groups of 400 additional genuine control cells, without replacement.
Comparison controls and the pooled ten null groups each match the inference
controls' batch totals and guide totals exactly. Minimize joint batch-by-guide
count reallocation subject to independent-cell availability, using a fixed
integer transportation problem. Report the complete joint deviation. Allocate
each stratum across ten groups as evenly as possible while keeping each group
at 400. Freeze source rows,
identities, stratum allocation and seed before reading selected expression rows.
Require enough unused cells for the recorded joint allocation; do not reuse cells.

Metadata-only amendment before expression extraction: exact joint matching
proved infeasible for a rare stratum (11 inference cells but only 12 unused,
where 22 were required). The initial attempt produced no cohort freeze or
expression output. It is retained as a failed metadata-selection attempt.
The revised rule preserves both marginal distributions exactly and finds the
smallest feasible joint change. It is not an outcome-driven adjustment.

Use the first ten frozen development target names as *technical aliases only*,
so the same feature-resolution convention can be used. These measured groups
contain no gene knockdown. Preserve original non-targeting labels in metadata
and clearly mark every artifact as a null diagnostic. No reserved or previously
closed perturbation outcomes are read. Selection seed is 20261023; the generator
uses the existing control-template seed 20261021.

## Arms

1. **Measured null:** the ten independent real-control groups, compared with
   the separate 4,000 comparison controls. This measures the empirical null.
2. **Donor replay:** choose one actual donor from each of the original generator's
   400 four-cell donor pools, using a fixed per-alias seed. Replay raw counts.
3. **Pool4 without moment fitting:** use the original pooled probability template
   and pooled depths, applying the original systematic row rounding and bulk
   repair but omitting moment fitting.
4. **Neutral dual emitter:** the unchanged two-moment emitter, requesting the
   inference-control arithmetic and bulk profiles, with no response-source effects.
5. **Existing response probe:** the original Model4 counts for the same ten
   development labels. This only measures overlap between its DE calls and
   zero-response emission calls under the same comparison controls. It is not
   judged against the null groups as if they were genuine perturbation outcomes.

The two pooled arms use exactly the same donor template and row depths. Donor
replay retains each real donor's actual depth, so it is a comparison of complete
emission policies, not a dispersion-only intervention. Real null cells are
independent of the inference/comparison cells. Generator replicates share a
template, inference pool and comparison controls; they are correlated technical
diagnostics, not ten independent biological experiments.

## Evaluation and interpretation

Use the pinned vcc2026 normalization, control-mean-CPM gate, Wilcoxon/scanpy,
per-target BH and strict adjusted p < 0.05. Run only the diagnostic significant
gene-count metrics and export DE tables. Do not compute normalized MSE, fidelity
or an overall score on fake perturbations. The measured-null run supplies a
verified reference-DE table for subsequent arms. Verify all data/config/cache
identities before reuse and all terminal statuses before reporting results.

Report per-alias significant calls, common tested-gene counts, strict q calls
excluding the alias gene, effect magnitudes, detection/variance diagnostics, and
overlap of neutral-emitter calls with the existing response probe. An overlap
does not establish that a particular full-model call is biologically false.

Preregistered diagnostic flag: neutral dual median calls exceed both measured
null and donor replay by at least tenfold (denominators floored at one), exceed
1% of the shared tested-gene universe, and exceed the corresponding donor replay
in at least eight of ten aliases. This is an engineering triage threshold, not
a calibrated biological hypothesis-test significance level. A high measured
null background prevents a clean attribution and must be reported.

If neutral emission strongly exceeds the empirical null, prioritize an emitter
that preserves the control distribution at zero effect before further response
parameter searches. If it does not, prioritize learning response identity and
direction across source contexts. Neither branch licenses a new model promotion
or submission without actual perturbation validation.

## Reassessment of progress

The original Model4 improved the official score, but the two later ablation
studies produced no accepted replacement. Negative experiments narrow the
explanation; they are not new leaderboard gains. The all-metric 1e-8 point gate
is stricter than a practical multi-objective selection rule. Future work should
report effect sizes, uncertainty and the tradeoff frontier explicitly, while
preserving the user's requested no-regression objective and all historical
verdicts. This null diagnostic changes no promotion criterion.
