# VCC 2026 public methods and actionable data, checked 2026-10-20

The immediate opportunity is to combine measured, same-target responses from several public cell contexts with correct arithmetic-CPM and raw-pseudobulk summaries, then generate integer cells that preserve those summaries. The public methods inspected do not justify claiming that a newly trained deep network will outperform this approach today.

## Current competition and leaderboard

The official leaderboard snapshot in `reports/leaderboard_20261020.json` contains 1,591 entries. Illumina AI / iCell leads with 0.473365; rank 20 is 0.309431 and rank 100 is 0.262157. The leader's scaled metrics are PDS 0.873147, MSE 0.792176, NMAE 0.439717, fidelity 0.185043, reach 0.470785 and Jaccard 0.079322. The public API supplied no method description for that entry. I found no released iCell implementation or checkpoint in the targeted search. These current scores are substantially above the inspected open implementations below. [Official leaderboard](https://virtualcellchallenge.org/leaderboard), [public API](https://virtualcellchallenge.org/api/leaderboard?get_final=false).

This year's task holds out entire cell contexts, provides only non-targeting controls, and accepts predictions rather than requiring published weights. Final contexts are scheduled for October 22; the final deadline is November 5. The public validation leaderboard therefore does not establish final-context generalization. [Official challenge](https://virtualcellchallenge.org/), [Arc announcement](https://arcinstitute.org/news/virtual-cell-challenge-2026).

The authoritative implementation is `ArcInstitute/cell-eval2` with the complete `vcc2026` preset. Six scaled metrics cover pseudobulk separation/error and differential-expression magnitude, directions and called sets. Correct normalization and integer cell emission affect the metrics independently of response-direction quality. Local scoring must keep the exact software, preset and reference-anchor provenance; earlier `cell-eval` metrics do not reproduce 2026 ranking. [Official scorer](https://github.com/ArcInstitute/cell-eval2).

## Inspected open methods

| Method | Reported submitted score | Useful implementation | Limitation |
|---|---:|---|---|
| AtlasShift | 0.154562, September 10 | Source-centered response ensemble; separate mean-CPM and pseudobulk heads; constrained integer emission | Historical author-reported score; five source banks required |
| Sidechain SER-19 | 0.1684, October 16 | Shrinkage, depth-weighted bulk anchor, real-control shape with selected DE-callable genes | Many successful settings were selected through repeated official submissions |
| rpinho v5b | 0.1528, through October 18 | AtlasShift anchors; spread calibration; source-sign gating; control-axis removal | v6 generalization has no recorded official score |
| zero-shot-perturbation | 0.137, October 2 | Broad source coverage; norm-restored transfer; explicit held-out-context studies | Local metric improvements do not always predict leaderboard changes |

AtlasShift uses source weights K562/HCT116/HEK293T/H1 = 2/1/1/2, with CD4 weight 0.5. It shrinks source counts toward matched controls using 100,000 prior counts, subtracts each source's common response, and uses amplitudes 0.6 for log2 arithmetic-CPM fold change and 0.3 for log1p pseudobulk change. Templates pool four depth-sorted control donors per predicted cell. The count solver fits both moments and repairs integer column totals. These are starting parameters, not validated optima for our source banks. [AtlasShift source](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026).

Sidechain's latest reported entry narrows variation for 2,700 selected genes per perturbation while preserving a control-derived cell shape elsewhere. Its previous real-cell-shaped submission had a strongly negative fidelity term, despite unchanged mean-based terms; the result demonstrates why matching means alone is insufficient. Its bulk anchor explicitly distinguishes the depth-weighted control pseudobulk from the equal-cell mean CPM. [Sidechain models](https://saberhq.com/sidechain/models/), [MIT implementation](https://github.com/saberhq/sidechain).

rpinho's workbench records improvements from spread calibration, cross-source direction agreement and projection away from control variation. It explicitly separates synthetic verification, local proxy estimates and officially submitted scores. Its public release offers useful bounded-memory scoring and source-bank acquisition code. [Repository and score history](https://github.com/rpinho/vcc2026).

The zero-shot-perturbation study reports that the essential-gene panels add no same-target coverage for the 300 validation targets. The Replogle genome-wide source covers most targets; Orion adds complete target coverage and many missing output genes. Its cross-context tests also show that correcting weak knockdowns by simply dividing effects by efficacy magnifies noise, and that tuning source weights yields limited gains once the transferable response has been recovered. [Repository and empirical results](https://github.com/MMatinGerami/zero-shot-perturbation).

Stage1-lab provides mathematical analysis of how DE significance can change at fixed means and why decoding errors can dominate scores. No code license was found in the inspected repository, and none of its implementation was copied. Its numerical claims are author-reported, not our independent verification. [Stage1-lab](https://github.com/danielchen26/vcc2026-stage1-lab).

## Newly acquired source data

The official Arc PIE release provides two Orion DE tables totaling 5.22 GB, with arithmetic CP10k target/control means, cell counts, fold changes and p-values. The pipeline filters for substantial knockdown before computing these tables. The resulting union-target banks cover 146 official targets in HCT116 and 11 in HEK293T; the lower HEK293T coverage is a measured limitation, not missing values to fill with zeros. Full-screen common responses are computed before selecting the 575-target union. [Official Arc PIE Orion data](https://huggingface.co/datasets/arcinstitute/PIE_x_atlas_orion), [processing configuration](https://github.com/ArcInstitute/pie/blob/main/src/pie/configs/process.yaml).

The source table gives arithmetic CPM but does not give library-weighted pseudobulk. Accordingly, `PIE_*_DE_statistics.npz` explicitly marks its bulk effect as an approximation. Its `measured` matrix is per target/gene; absent records are excluded from fusion. Both source files were verified against Hugging Face SHA256 metadata. Reports are in `reports/PIE_*_DE_statistics.json`.

The original Orion release contains 109 HCT116 and 223 HEK293T Parquet files, totaling 126.26 GB. `build_orion_raw_bank.py` builds true raw-count and arithmetic-CPM statistics, matches controls by experimental batch, and retains the union of all 300 H1 targets and all 300 official targets (575 unique genes). Each downloaded batch is checksum verified; source libraries are checked against count sums; duplicate symbol loci are summed explicitly. The dataset license is CC BY-NC-SA 4.0. [Original Xaira release](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion).

An independent scalar/list-based check of the first batch in each line reproduced raw sums exactly and CPM sums within 1e-12 of the streaming sparse implementation. Evidence: `reports/orion_raw_builder_independent_check.json`. Completed batches are resumable in `data/processed/orion_batches/`; final source reports establish actual coverage and cell counts.

All original-count acquisition and processing finished on October 20. Both final banks passed count/probability/mass/mask checks and every source-file SHA256 receipt, and both loaded successfully through the AtlasShift response interface. Their reports contain the artifact SHA256 and an explicit completed audit.

| Completed bank | Treated cells retained | Controls | Official targets present / ≥20 cells | H1 targets present / ≥20 cells | Median cells per target |
|---|---:|---:|---:|---:|---:|
| `data/processed/HCT116_full_statistics.npz` | 89,839 | 165,777 | 300 / 290 | 299 / 289 | 133 |
| `data/processed/HEK293T_full_statistics.npz` | 120,987 | 218,838 | 300 / 294 | 299 / 292 | 185 |

Both banks use the same 18,536-gene union axis, with 18,106 measured source symbols. Targets below the cell threshold remain explicitly counted and can be excluded by the model; their measurements are never silently recoded as absent biological effects. Final receipts: `reports/HCT116_full_statistics.json` and `reports/HEK293T_full_statistics.json`. The observed control counts independently agree with those in Arc's published Orion dataset card.

The optional CD4 bank also finished at approximately 22:08 HKT on October 20, before its 22:20 cutoff. Public S3 range requests retrieved 251 parts with eight concurrent transfers; assembly verified HTTP 206, exact byte ranges, lengths and a consistent ETag for every part. The complete 16,786,240,107-byte file matched the published SHA256 `c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62`. It resides at `data/raw/GWCD4i.DE_stats.ranged.h5ad`; the pre-existing full-file download was not modified.

`data/CD4_DE_statistics.npz` retains 575 target slots and 10,282 publisher genes across Rest, Stim8hr and Stim48hr. After the publisher's guide/knockdown/off-target/expression quality gates and a 20-cell minimum, it covers 251 official targets and 223 H1 targets. The finite-effect, availability and quality-flag checks passed, as did the actual AtlasShift CD4 adapter. Acquisition and bank evidence are in `reports/cd4_range_acquisition.json`, `reports/CD4_DE_statistics.json`, `reports/cd4_quality_flag_schema.json`, and `reports/cd4_atlasshift_integration.json`. The raw acquisition URL and immutable S3 version are retained in those receipts. Its contribution remains a model-selection question; source-bank completion is not evidence of a higher challenge score.

## Concrete model recommendation

1. Use all public H1 data for the final training bank, but exclude it completely when evaluating transfer to held-out H1 contexts.
2. Combine H1, Replogle GWPS and both raw Orion banks. Evaluate the strongly filtered PIE Orion bank as an additional, separately weighted source instead of silently replacing weak raw targets.
3. Use arithmetic CPM changes for DE and raw-pseudobulk changes for PDS/MSE; keep missing gene/target measurements masked.
4. Calibrate effect amplitudes and emission dispersion on genuine held-out H1 responses using the official scorer. Do not infer superiority from a reconstruction test or same-context training score.
5. Generate all 360,000 count profiles, verify the realized means and distribution-derived DE behavior, package with the current official CLI and submit an artifact selected by this validation.

These steps directly address the measured weaknesses of the earlier models. They do not establish that the result will match the current 0.473365 leader; only an official submission can measure that comparison.

## Source revisions retained locally

| Directory | Revision | Code license |
|---|---|---|
| `external/atlasshift` | `d24ce4fdae0cd7cba3ba29546cd8737094eae98a` | MIT |
| `external/sidechain` | `6a84a26980af45f72c83fe51c4848679c56ea07f` | MIT |
| `external/zero-shot-perturbation` | `2a198d5fa165179633ed0c0becc0d2b050b27187` | MIT |
| `external/rpinho-vcc2026` | `e00c9b5575b89c0e150e0480918afc6218861739` | MIT for included original/permissive code |
| Arc PIE Orion data | `ecb4491ea49fcf28df463fac1fc9ef00c0c3ac17` | Original Orion noncommercial/share-alike terms apply |

Cloning source did not execute it. Our two source-bank scripts are new local implementations; they read public data and preserve its provenance.
