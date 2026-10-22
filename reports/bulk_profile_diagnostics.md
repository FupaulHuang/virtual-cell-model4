# Frozen bulk-arm profile diagnostics

All four complete 50-target outputs passed this independent engineering audit (20,000 cells × 18,080 genes each). The new incumbent is identical to the original incumbent in every integer count and observation/gene axis. Actual H5AD means, raw gene totals and row depths were recomputed; saved realized/old-incumbent profiles and CPM changes agree exactly. Input artifact and frozen dependency hashes were verified. No reference-response values, metric results or reserved targets were read.

Requested arithmetic profiles are bit-identical across all arms for all 904,000 gene-target pairs. CD4, control-template and depth hashes also match. Every actual row depth matches the original. Variant realized arithmetic means differ after fitting and rounding; this experiment does not impose the separate dispersion study’s CPM bound.

| Arm | Projection L1 median / max | CPM drift median across gene-target pairs / max | Requested → projected → realized bulk intervention L1 (target median) |
| --- | ---: | ---: | ---: |
| incumbent | 0.0001154767 / 0.0002876678 | 0 / 0 | 0 → 0 → 0 |
| common05 | 0.0001134962 / 0.0002850087 | 0.0753597 / 1.01619 | 0.0007840394 → 0.0007808384 → 0.0007766747 |
| bulk015 | 4.192659e-05 / 0.0001063298 | 0.0793849 / 0.924975 | 0.006408952 → 0.006361223 → 0.00630419 |
| bulk060 | 0.0003206967 / 0.0008069445 | 0.0783569 / 1.01548 | 0.01274143 → 0.01261612 → 0.01254972 |

Bulk intervention magnitudes compare each arm with the new incumbent at the corresponding stage. They use the full gene axis, including the on-target gene, and are not biological error measures. Exact per-target probability L1/L2 and log1p(50,000×probability) RMS changes are in the JSON report (local artifact omitted: `bulk_profile_diagnostics.json`).

The unchanged projection limit is L1 ≤ 0.03 and each fitting limit is ≤ 0.001. An independent projection calculation reproduced the saved projection within floating-point precision. Projection effects are nonuniform across genes and targets:

- **incumbent:** requested genes outside bounds per target range 112–1723 (median 1376); the ten largest gene changes account for a median 7.14% of projection L1. Median L1 distortion of the intervention attributable to projection is 0. Absolute pre-round-to-integer CPM change has pooled median 0.060167, maximum 0.697574.
- **common05:** requested genes outside bounds per target range 112–1718 (median 1368); the ten largest gene changes account for a median 7.21% of projection L1. Median L1 distortion of the intervention attributable to projection is 4.170152e-06. Absolute pre-round-to-integer CPM change has pooled median 0.060275, maximum 0.738883.
- **bulk015:** requested genes outside bounds per target range 85–1199 (median 945); the ten largest gene changes account for a median 8.64% of projection L1. Median L1 distortion of the intervention attributable to projection is 7.434348e-05. Absolute pre-round-to-integer CPM change has pooled median 0.0609442, maximum 0.679311.
- **bulk060:** requested genes outside bounds per target range 160–2362 (median 1908); the ten largest gene changes account for a median 6.08% of projection L1. Median L1 distortion of the intervention attributable to projection is 0.0002015659. Absolute pre-round-to-integer CPM change has pooled median 0.0583543, maximum 0.688485.

The complete result supports mechanical isolation of the requested bulk intervention, not unchanged downstream DE scores. Projection, finite fitting and integer rounding can alter its realized magnitude and arithmetic CPM. Six-metric and real-reference residual conclusions belong to the separate analysis; no arm is selected by this report.
