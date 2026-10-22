As of **2026-10-21T02:18:22.675372+00:00** (10:18 HKT on October 21), the current top 10 mean is **0.341648**, versus Model4 **0.147019**. The clearest measured explanation is stronger expression-profile accuracy and direction reach, together with better perturbation discrimination. Exact current-entry methods remain largely undisclosed; score superiority is more strongly supported than any architectural explanation.

**Ranking audit**

The fresh [official API](https://virtualcellchallenge.org/api/leaderboard?get_final=false) returned 1,608 matching published rows and 1,608 unique team IDs, with **zero duplicate team rows**. The top 10 returned entries and top 10 unique teams therefore coincide in this snapshot. This feed does not expose all historical submissions, so it cannot establish an all-history top 10-entry ranking. The leader now displays a new iCell entry at 0.471613; the saved October 20 snapshot displayed 0.473365. We do not silently substitute a historical best.

All compared rows use `vcc2026-val-1` and `vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4`. Each top 10 overall score matches the mean of the six scaled components to floating-point precision. Model4 is rank 552 in this snapshot; its stored publication rank 565 is an earlier observation, not a new model result.

**Top 10 unique teams: all scaled components are higher-is-better**

| Rank | Team / current model | Overall | PDS | MSE | NMAE | Fidelity | Reach | Jaccard | API submissions |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | Illumina AI / iCell | 0.471613 | 0.877375 | 0.792614 | 0.432061 | 0.180695 | 0.466971 | 0.079961 | 93 |
| 2 | X-team / X-97 | 0.340832 | 0.773266 | 0.507952 | 0.294379 | -0.007735 | 0.378019 | 0.099114 | 53 |
| 3 | cqawesome / m1020c | 0.332096 | 0.831500 | 0.461532 | 0.288508 | 0.088982 | 0.296349 | 0.025704 | 120 |
| 4 | FlowState / flowstate | 0.331676 | 0.817751 | 0.513521 | 0.257460 | 0.090255 | 0.307276 | 0.003789 | 37 |
| 5 | realcell / test | 0.329543 | 0.838449 | 0.512292 | 0.239948 | 0.072221 | 0.310692 | 0.003655 | 41 |
| 6 | TAO / v72 | 0.326720 | 0.830083 | 0.499470 | 0.241512 | 0.057995 | 0.314609 | 0.016652 | 60 |
| 7 | Aginglab.com / VirtualOmics | 0.324455 | 0.813830 | 0.524436 | 0.253895 | 0.044103 | 0.285906 | 0.024560 | 88 |
| 8 | Chenhao Tan / j94new4 | 0.323554 | 0.830969 | 0.511142 | 0.256510 | 0.054017 | 0.274058 | 0.014631 | 8 |
| 9 | Hopkin Jiang / bbbbbbbbulk | 0.318462 | 0.831974 | 0.454417 | 0.251918 | 0.067855 | 0.286106 | 0.018499 | 15 |
| 10 | C-team / The_Way53 | 0.317528 | 0.832494 | 0.469382 | 0.214200 | 0.092907 | 0.276186 | 0.020000 | 93 |
| — | Model4 | 0.147019 | 0.645588 | 0.037668 | 0.138194 | -0.018683 | 0.079313 | 0.000032 | — |

**Top 10 returned entries, independently ranked**

| Entry rank | Entry ID | Team | Overall |
| ---: | --- | --- | ---: |
| 1 | `[PRIVATE_IDENTIFIER_REMOVED]` | Illumina AI | 0.471613 |
| 2 | `[PRIVATE_IDENTIFIER_REMOVED]` | X-team | 0.340832 |
| 3 | `[PRIVATE_IDENTIFIER_REMOVED]` | cqawesome | 0.332096 |
| 4 | `[PRIVATE_IDENTIFIER_REMOVED]` | FlowState | 0.331676 |
| 5 | `[PRIVATE_IDENTIFIER_REMOVED]` | realcell | 0.329543 |
| 6 | `[PRIVATE_IDENTIFIER_REMOVED]` | TAO | 0.326720 |
| 7 | `[PRIVATE_IDENTIFIER_REMOVED]` | Aginglab.com | 0.324455 |
| 8 | `[PRIVATE_IDENTIFIER_REMOVED]` | Chenhao Tan | 0.323554 |
| 9 | `[PRIVATE_IDENTIFIER_REMOVED]` | Hopkin Jiang | 0.318462 |
| 10 | `[PRIVATE_IDENTIFIER_REMOVED]` | C-team | 0.317528 |

**Where the score advantage comes from**

The gap to the top 10 mean is **0.194629**. For each component, contribution = (top 10 mean component − Model4 component)/6. These contributions sum exactly to the overall gap; they describe scores, not causes. The cutoff for this top 10 is 0.317528, which is 0.170509 above Model4. [Primary score source](https://virtualcellchallenge.org/api/leaderboard?get_final=false).

| Component | Model4 | Top 10 mean | Contribution to overall gap | Share of gap |
| --- | ---: | ---: | ---: | ---: |
| MSE | 0.037668 | 0.524676 | 0.081168 | 41.7% |
| Reach | 0.079313 | 0.319617 | 0.040051 | 20.6% |
| PDS | 0.645588 | 0.827769 | 0.030364 | 15.6% |
| NMAE | 0.138194 | 0.273039 | 0.022474 | 11.5% |
| Fidelity | -0.018683 | 0.074129 | 0.015469 | 7.9% |
| Jaccard | 0.000032 | 0.030657 | 0.005104 | 2.6% |

MSE plus reach explain **62.3%** of the mean gap. Every top 10 MSE component is at least 0.4544, while Model4 is 0.0377. Jaccard contributes only 2.6% of the mean gap, despite remaining a serious biological limitation. Many top teams also have small scaled Jaccard values; a high overall score is not evidence that all affected genes are accurately identified.

The official scorer uses perturbation-profile ranks, a sampling-corrected expression error, and DE-based magnitude/direction/set metrics. Its values are baseline-to-replicate scaled, not percentages of correct genes. PDS excludes the panel target genes, so simply reproducing the direct knockdown is not a sufficient explanation of high discrimination. Expression error also includes a bounded sampling correction; an exact score does not by itself distinguish improved mean prediction from better distribution calibration. [Official definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md).

*Contribution of each metric to the gap versus Model4 (figure excluded from the source-only export).*

**What primary sources actually substantiate**

All 10 current API descriptions are null. No exact, reproducible recipe tied to the current top 10 entry IDs was verified in this bounded search. The following related evidence must not be promoted into unsupported entry attribution:

| Team | Primary evidence | What remains unverified |
| --- | --- | --- |
| Illumina AI | A listed member describes genomic deep-learning research. No exact iCell architecture, weights, training list, or submission reconstruction was found. [team-member profile](https://anikethjr.github.io/) | Research experience does not identify the submitted algorithm. |
| X-team | The API lists Xaira. Its X-Cell repository describes a diffusion language model trained on 25.6 million perturbed cells with biological priors and iterative refinement. The README says model weights/inference code are coming soon. [X-Cell repository](https://github.com/Xaira-Therapeutics/X-Cell); [X-Cell model card](https://huggingface.co/Xaira-Therapeutics/X-Cell/blob/main/README.md?code=true) | No source ties entry X-97 to a specific X-Cell checkpoint, dataset version, ensemble or emitter. Do not equate the names. |
| cqawesome | An institutional page naming Chenxi Qiu describes VCC participation using specialized agents and scientific evaluation harnesses. [institutional research description](https://research.bidmc.org/student-experience/opportunities-research-students-constantly-updated) | The name matches the API-listed member, but the page does not link m1020c or supply its prediction recipe; this supports a workflow claim only. |
| FlowState | No exact current-entry primary method disclosure found; API description is null. | The name FlowState does not prove flow matching or any particular model family. |
| realcell | No exact current-entry primary method disclosure found; API description is null. | Architecture, data, calibration and ensembling remain unknown. |
| TAO | No exact current-entry primary method disclosure found; API description is null. | Affiliation and version labels are not architecture evidence. |
| Aginglab.com | The company describes harmonizing nine public sources, 10.72 million measured-cell records and a generative expression platform. [AgingLab platform](https://aginglab.com/platform.html); [AgingLab validation](https://aginglab.com/validation.html) | These are company-level disclosures, not a reconstructible CRISPRi VirtualOmics submission. Its internal PCC claim is not a VCC score. |
| Chenhao Tan | No exact current-entry primary method disclosure found; API description is null. | Do not infer identity or architecture from name similarity. |
| Hopkin Jiang | No exact current-entry primary method disclosure found; API description is null. | The model name bbbbbbbbulk is not evidence of a particular bulk-response algorithm. |
| C-team | API-listed members Dewei Hu and Marc Pielies Avelli authored SLIM, a compact STRING-embedding bilinear ridge model. [SLIM repository](https://github.com/RasmussenLab/SLIM); [SLIM preprint](https://www.biorxiv.org/content/10.64898/2026.08.07.743481v1) | SLIM evaluates held-out perturbations in its paper; no source identifies The_Way53 as SLIM or establishes its 2026 context-transfer recipe. |

The Xaira evidence makes large perturbational corpora, biological priors and iterative conditional prediction plausible mechanisms to investigate. The C-team-related SLIM work demonstrates that useful perturbation predictors can also be compact linear models. Neither finding tells us which mechanism produced X-97 or The_Way53. Company-scale data claims and internal PCC results are not challenge ablations or independent score reproductions.

**Referee assessment of robustness and validation feedback**

The top 10 API submission_count values span 8–120 (median 56.5). Their precise counter semantics and full search histories are not exposed here. Repeated validation feedback can favor settings adapted to A/B/C; it does not demonstrate that any team overfit, nor that any team violated rules. Model version names are not a substitute for search history. A referee should request a precommitted final recipe and independent-context evidence rather than infer misconduct from success or iteration.

The organizer allows public and team-owned training data. The final three contexts are scheduled for October 22 and differ from the three validation contexts, so current validation rank does not establish final-context generalization. [Official task and timeline](https://arcinstitute.org/news/virtual-cell-challenge-2026).

Evidence needed for stronger scientific claims:

- Exact frozen checkpoint, source-data/version inventory and inference recipe linked to each scored entry.
- Raw and scaled metrics per context and perturbation, with uncertainty and sensitivity to seeds/effect-strength/source-coverage strata.
- Ablations separating predicted mean/bulk-profile changes from cell-shape/dispersion changes at controlled means, library sizes and count constraints.
- Independent final-context or external-context results from a precommitted recipe; do not reuse public Validation adaptively and call it independent confirmation.
- Complete submission/selection history if making a claim of robustness to leaderboard feedback; counts alone cannot quantify selection bias.

**Implications for our next experiment**

- Prioritize held-out calibration of the bulk profile and effect magnitude: expression MSE contributes the largest measured gap, while preserving PDS and direction reach.
- Evaluate context-dependent response transfer rather than assuming a universal source fold change; this is a model hypothesis, not an attribution to top teams.
- Keep DE precision/fidelity as guardrails; matching the top 10 mean Jaccard alone would close only 2.6% of the observed overall gap.
- Broader related public models and compact biological representations are both plausible candidates; no evidence supports assuming that every leader is a giant foundation model.

This research used only public read-only sources. No model, frozen configuration, submission, authentication, or official-score receipt was changed. The machine-readable report (local artifact omitted: `../reports/top10_referee_20261021.json`) contains exact per-team raw/scaled gaps, source-evidence levels, historical comparisons and snapshot SHA256; the fresh public snapshot (local artifact omitted: `../reports/leaderboard_20261021.json`) preserves the observed data.
