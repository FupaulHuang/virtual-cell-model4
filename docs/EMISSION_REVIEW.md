# Count emission and H1 validation review

Reviewed on 2026-10-20. AtlasShift is pinned at
`d24ce4fdae0cd7cba3ba29546cd8737094eae98a`. External code was not edited.
Its model.py SHA-256 is
`8ad68112dbfa9c5c8c8b21a717dbe475fe23613f84c25e844fbfbf8bb47cff3b`.

## Main finding

The emitter resolves an important error in Model3's count generation. Model3
learned differences in **mean(log1p(CPM))**, used `exp(delta)` as raw-count fold
changes, and diagnosed **log1p(mean(CPM))**. Those three quantities are not
interchangeable. It also modified only already observed nonzeros, so it could
not increase detection for genes absent from a sampled donor cell. Model3's
submitted fidelity was 0.008386 and Jaccard was 0.000394; the generator mismatch
is a plausible contributor, not proof of the sole cause.

AtlasShift separately estimates an arithmetic mean CPM distribution and a
library-size-weighted pseudobulk distribution. Its count emitter fits both
distributions, adds support for genes absent from the sampled donor pool, then
rounds while preserving each cell's integer library depth. A final repair
enforces nearest integer pseudobulk gene totals. This matches the distinct
arithmetic-CPM DE and pseudobulk expression estimands in the local official
scorer.

Twenty-two independent tests passed using the existing cell-eval2 Python 3.12
environment. They cover distinct feasible moments, exact depth conservation,
zero-gene activation, requested zeros, integral nonnegative counts, seed
reproducibility, impossible moment rejection, H5AD roundtrips, the root Writer
with all 18,533 official genes, gene order, group completeness, count limits,
control-context provenance, source exclusion and receipt/hash corruption.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python \
  -m unittest discover -s tests -p test_count_emission.py -v
```

## CPU benchmark and operational configuration

Each benchmark emitted only one 400-cell neutral group using the complete
official 18,533-gene axis and four pooled donors per synthetic cell. No complete
challenge prediction was created by this review. One BLAS/OMP thread was used.

| Context | Emission seconds | Nonzeros | Mean-probability L1 error | Bulk-probability L1 error |
| --- | ---: | ---: | ---: | ---: |
| A | 1.158 | 2,805,325 | 0.006598 | 0.000434 |
| B | 1.171 | 2,599,365 | 0.012753 | 0.000460 |
| C | 1.155 | 2,745,156 | 0.006417 | 0.000453 |

All three preserved every requested cell depth exactly. Loading and preparing
each context took approximately 8.4–8.8 seconds. Peak process RSS including
controls was about 2.9 GiB. A synthetic 500-gene response with mean log2FC
amplitude 0.6 and bulk amplitude 0.3 also passed, taking 1.17 seconds with
2,792,557 nonzeros and similar errors. Full records are in
`reports/count_emission_benchmark.json`.

Neutral extrapolation gives 2.445 billion nonzeros, 51.5% of the official
4.75-billion-entry cap, and about 17.4 minutes of single-worker emission for
900 groups, excluding writing, packaging and hashing. This is an estimate,
not a full-run measurement. Difficult moment constraints can take longer or
fail; real perturbation density must still be checked.

The root wrapper's `--cells 400 --pool 4 --workers 4` is feasible on the
observed CPU environment if approximately 8–12 GiB of RAM is available to the
parent plus workers. Keep BLAS/OMP/MKL at one thread per worker and bound queued
dense count blocks. Stream CSR output, keep int64 CSR pointers, and validate the
final file through `scripts/validate_prediction_receipt.py` and official prep.
The tested numerical defaults are 100 iterations, fitting tolerance 0.0002,
mean effect amplitude 0.6, bulk effect amplitude 0.3 and effect clip 3; only the
first two are numerical settings. The effect amplitudes need source-context
validation and are not established as optimal by these tests.

## Limits that must remain visible

- The fitting tolerance applies before integer rounding. In the real-axis
  benchmark the final mean-probability error is 0.0064–0.0128, much larger than
  0.0002. The lower-depth cells in B amplify integer error. Record final moments;
  do not claim the pre-rounding tolerance was preserved.
- Simultaneously requesting mean and bulk profiles imposes depth-dependent
  feasibility bounds. The implementation projects bulk toward those bounds
  and rejects L1 projection changes over 0.03. Validate every candidate's
  realized moments before full submission.
- Pooling four donors reduces biological variance; this can alter Wilcoxon
  significance even when means are correct. A pooling ablation must accompany
  claims about improved DE accuracy. Cell-count invariants alone cannot prove
  realistic dispersion or biological predictive accuracy.
- Activation is possible when a gene has positive requested mean but is absent
  from sampled donor cells. The default multiplicative mean-profile model does
  not activate genes that are zero in *all* context controls. A bulk-only
  positive value cannot override a zero requested arithmetic mean.
- `desired_mean(..., space='log2fc')` applies a multiplicative effect derived
  using a CPM pseudocount; it does not invert that pseudocount exactly. This is
  deliberate shrinkage for low-expression genes, and differs from the Model3
  log-mean mismatch, but still needs validation.

## Frozen H1 comparison protocol

`scripts/prepare_h1_validation.py` materializes
`evaluations/h1_validation/` from Model3's existing frozen reference. It retains
50 public H1 Validation targets, 18,080 genes, 2,000 controls and 200 real cells
per target except MAX, which has 161. The copied reference SHA-256 is
`130bf27286ebbcc76424d6beb55609f2aa2ceb67af0e06d9032cc5988d61ac7a`.
Generate 200 predicted cells per target. The 2,000 controls support either
pooling one or four donors. H1 Test remains untouched.

For a **leave-context-out** comparison, exclude the entire H1 response bank.
Use only H1 controls, the known target names and external-source responses.
The independent validator accepts `--exclude-source-token H1` for this gate.
Known targets in other source contexts are allowed: this matches the 2026
task. Audit direct source coverage by target and report covered/uncovered
strata. The source filenames alone are a useful gate, not a complete data
leakage audit; source manifests must support their stated provenance.

Compare these arms with identical source effects, control draws and seeds:

1. Neutral control-moment emission (`scripts/build_control_prediction.py`).
2. An old-style donor multiplier using correctly expressed fold changes,
   with pooling one, to isolate count-generation effects.
3. Dual-moment emission with pooling one.
4. Dual-moment emission with pooling four, as a separate variance ablation.

The saved historical Model3 prediction is useful as an additional comparator,
but it used H1 Training responses and therefore is not a leave-context-out
model. It must not be described as having the same source exclusion.

Recompute every arm, the mean-response baseline and five replicate anchors
using the same CPU `scanpy` backend and the exact same reference. Existing
Model3 metrics used CUDA/gpudge, so do not directly compare those old numbers
to new CPU runs. Use the official `vcc2026` profile; report all six raw/scaled
components, significant-set sizes and final moment errors. Once a candidate
is frozen on the 50-target validation panel, evaluate the 100-target H1 Test
panel once; do not tune repeatedly on Test.

Actual CLI testing revealed that cell-eval2 0.16.0's initial `validate_pair`
requires the control label in both inputs, even with `control_source=real`.
The metric run wrapper now uses `scripts/prepare_metric_input.py` to append
the exact reference controls to an isolated local scoring copy. Original
prediction files and their hashes remain unchanged. The copied controls pass
the input gate; the scorer subsequently uses its real-side controls as
specified. These scoring copies must never be submitted to the challenge.
