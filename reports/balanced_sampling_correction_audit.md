# Post hoc arithmetic sampling-correction audit

For each policy, compare the actual error over all200 frozen draws with the nominal independent-cell sample-variance/n term on population-CPM>5 genes. This is **not the official log-bulk/global-capped MSE correction** and changes no gate.

| Context | Policy | Actual mean-error sum (CPM²) | Nominal iid term (CPM²) | Nominal / actual |
| --- | --- | ---: | ---: | ---: |
| H1 | uniform | 94122.425 | 103631.857 | 1.1010 |
| H1 | depth_stratified | 89528.501 | 103424.196 | 1.1552 |
| H1 | balanced | 77357.616 | 103723.095 | 1.3408 |
| A | uniform | 380713.010 | 422341.313 | 1.1093 |
| A | depth_stratified | 370641.043 | 421699.908 | 1.1378 |
| A | balanced | 281334.437 | 421539.072 | 1.4984 |
| B | uniform | 433433.558 | 475820.525 | 1.0978 |
| B | depth_stratified | 416136.663 | 477999.598 | 1.1487 |
| B | balanced | 359024.638 | 477445.717 | 1.3298 |
| C | uniform | 357199.142 | 399876.729 | 1.1195 |
| C | depth_stratified | 356058.770 | 396093.144 | 1.1124 |
| C | balanced | 299923.929 | 399932.978 | 1.3334 |

Without-replacement and auxiliary-balanced draws are dependent designs. Individual-cell sample variance/n is an iid estimator, not the exact variance of these sample means. Uniform sampling already needs the0.9 finite-population correction; balancing changes additional covariance terms. A similar cell-variance trace does not imply the same group-mean sampling error.

The nominal term is recovered without rereading counts: the stored variance-trace ratio multiplies the known population variance. Specifically, `nominal = mean(trace_ratio) × sum(uniform_exact_variance) /0.9`. Actual error is the sum of per-gene Monte Carlo mean squared errors, with matching axes and gates verified.

An iid correction can exceed actual balanced-sample mean error even when cell-level variances look unchanged. This does not demonstrate inflation of the official score: its transformed bulk jackknife and global cap are different calculations. The finite-draw mean bias in the JSON is simulation noise, not a finding of systematic selection bias.

Full formulas, values and input hashes (local artifact omitted: `balanced_sampling_correction_audit.json`). No held-out oracle performance was read.
