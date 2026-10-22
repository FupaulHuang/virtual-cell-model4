# Provisional one-split H1 comparison

**Diagnostic only. The independent five-split calibration is still running.**

| Arm | Provisional six-component mean | pds_cosine | expr_mse_unbiased_capped_norm | de_wilcoxon_lfc_nmae | de_wilcoxon_direction_fidelity_yield_raw | de_wilcoxon_direction_reach_raw | de_wilcoxon_sig_jaccard |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| neutral_cpu | -0.021361 | -0.011916 | 0.000000 | -0.120518 | 0.094889 | -0.095732 | 0.005109 |
| pie_dual_cpu | 0.107191 | 0.368469 | 0.000000 | -0.077034 | 0.259883 | 0.014014 | 0.077813 |
| pie_single_cpu | 0.070968 | 0.366636 | 0.000000 | -0.100503 | 0.212882 | -0.096941 | 0.043734 |
| pie_dual_pool1_cpu | 0.046211 | 0.257562 | 0.000000 | -0.076974 | 0.050845 | 0.009253 | 0.036582 |

All values use official per-component replica scaling against the existing CPU baseline and exactly one official split. The reported mean is explicitly provisional; it is not a production challenge score or the completed five-split result.
