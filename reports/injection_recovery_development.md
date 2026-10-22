# Moment-oracle development results

This stage evaluates known response moments against independently generated reference cells. It is not a full-distribution oracle, a biological-transfer test or an official score. All candidate methods and requests were frozen before DE outcomes; no adaptation is permitted between stages.

| Family | Seed | Alias | Method | Reference calls | Candidate calls | Reference signed recovery | Candidate signed recovery | NULL/marginal-null calls |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| composition | 20261026 | COMP_T1_NULL | current | 0 | 4428 | undefined | undefined | 4428 |
| composition | 20261026 | COMP_T1_NULL | preserving_uniform | 0 | 0 | undefined | undefined | 0 |
| composition | 20261026 | COMP_T1_WEAK | current | 1 | 4669 | undefined | undefined | 0 |
| composition | 20261026 | COMP_T1_WEAK | preserving_uniform | 1 | 0 | undefined | undefined | 0 |
| composition | 20261026 | COMP_T1_MODERATE | current | 1 | 5206 | 0.020 | 0.857 | 0 |
| composition | 20261026 | COMP_T1_MODERATE | preserving_uniform | 1 | 2 | 0.020 | 0.020 | 0 |
| composition | 20261026 | COMP_T1_STRONG | current | 18 | 5719 | 0.047 | 0.787 | 0 |
| composition | 20261026 | COMP_T1_STRONG | preserving_uniform | 18 | 28 | 0.047 | 0.054 | 0 |
| dirichlet_multinomial | 20261026 | DM_T1_NULL | current | 0 | 10438 | undefined | undefined | 10438 |
| dirichlet_multinomial | 20261026 | DM_T1_NULL | preserving_uniform | 0 | 0 | undefined | undefined | 0 |
| dirichlet_multinomial | 20261026 | DM_T1_WEAK | current | 4 | 10430 | 0.062 | 0.891 | 10373 |
| dirichlet_multinomial | 20261026 | DM_T1_WEAK | preserving_uniform | 4 | 2 | 0.062 | 0.031 | 0 |
| dirichlet_multinomial | 20261026 | DM_T1_MODERATE | current | 57 | 10430 | 0.875 | 0.891 | 10373 |
| dirichlet_multinomial | 20261026 | DM_T1_MODERATE | preserving_uniform | 57 | 34 | 0.875 | 0.469 | 4 |
| dirichlet_multinomial | 20261026 | DM_T1_STRONG | current | 65 | 10437 | 1.000 | 1.000 | 10373 |
| dirichlet_multinomial | 20261026 | DM_T1_STRONG | preserving_uniform | 65 | 51 | 1.000 | 0.750 | 3 |

Signed recovery uses oracle |log2FC|≥0.1 within the realized comparison-CPM>5 gate. Below-threshold calls are reported without calling all of them false. NULL/marginal-null calls use the independent true-law mask, not merely a small aggregate mean. Undefined denominators remain undefined.

Reference NULL validity: composition/20261026: provisional_opened_null_pass (maximum0calls); dirichlet_multinomial/20261026: provisional_opened_null_pass (maximum0calls).

The detailed JSON/CSV retain effect strata, out-of-gate true signals, reference call power, overlap, finite-baseline nuisance, per-gene request/clipping provenance and population/sample/realized profile errors. Current-emitter population-fit statistics are unavailable rather than fabricated. Candidate profile values were checked against actual integer counts.

A weak result can reflect anchoring distortion, finite-control uncertainty, missing higher-order distribution information or generator-family limitations. It is not automatically an implementation bug. No scalar success metric or biological promotion decision is inferred.

This immutable report opens only template1/seed20261026. Full family/seed reference-null validity remains pending confirmation; held-out templates/seeds are not used for a development decision.
