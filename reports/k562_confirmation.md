The frozen selected model improves **6 of six raw metrics** over the zero-response control comparator on the reserved K562 cohort (0 worse, 0 tied, 0 undefined). No anchor-scaled overall score is computed.

The model was frozen at `2026-10-20T14:46:16.606761+00:00`; confirmation started at `2026-10-20T14:52:42.181326+00:00`. Both independent validations excluded K562 and GWPS response sources, and both metric processes exited 0.

| Metric | Better | Selected | Neutral | Selected − neutral | Result |
| --- | --- | ---: | ---: | ---: | --- |
| PDS | higher | 0.753061 | 0.500000 | 0.253061 | improved |
| MSE | lower | 1.206125 | 1.307675 | -0.101550 | improved |
| NMAE | lower | 0.876944 | 1.002672 | -0.125728 | improved |
| Fidelity | higher | 0.498891 | 0.443680 | 0.055211 | improved |
| Reach | higher | 0.458733 | 0.289018 | 0.169715 | improved |
| Jaccard | higher | 0.027044 | 0.010656 | 0.016387 | improved |

Selected source coverage is 50/50; neutral coverage is 0/50. The selected predictor uses H1 × 2, HCT116 × 1, HEK293T × 1 and CD4 DE × 0.5. The comparator keeps an HCT116 input at weight 0; it contributes no perturbation response.

The frozen cohort contains 8,827 genuine treated cells, 2,000 measured controls, 50 targets and 8,246 unique genes. Target selection used only the fixed seed and the eligibility rule of at least 100 cells. Source row identities and GEM batches are preserved.

Limits:

- This is one public K562 context on 50 targets selected with a fixed seed among those with at least 100 cached cells; it does not represent every target or the hidden 2026 A/B/C contexts.
- The native panel has 8,246 unique measured genes, versus 18,080 in H1 and 18,533 in the challenge. HSPA14/TBCE duplicate source columns were summed; library normalization uses this native gene axis.
- Real groups contain 103–200 cells; predictions contain 200 cells per target. Exactly 2,000 measured controls are shared between inference and the reference.
- NMAE is defined for only 12/50 targets and reach for 34/50; their improvements apply to those eligible subsets. Finite target identities are recorded for both arms. Normalized expression MSE is aggregate-derived and has no per-target output rows.
- The reference retains 267 GEM batch labels, but this vcc2026 metric run compares against pooled real controls and does not adjust its DE model for batch effects.
- K562 raw responses were a public source for the earlier H1 validation comparison. All K562/GWPS response representations are excluded here; these metrics were unseen when parameters were frozen. This evaluates transfer across contexts using existing public data.
- The neutral comparator has zero response-source weights and preserves control moments. It is not the official generic-response baseline.
- No anchors were built or reused. There is no reference-scaled overall score, and heterogeneous raw metrics are not averaged.
- This confirmation is evaluated once after model freeze; no parameter retuning, model selection, or submission was performed by this workflow.

The JSON contains all six values, finite-target counts, coverage, frozen parameters, source/prediction/reference hashes, source-row provenance, and authoritative process exit codes. No H1 anchors or K562-based retuning were used.
