# Posthoc nonzero-response distribution audit

All eight frozen predictions were checked against their receipts and the actual 400-cell-per-target measured reference. All count values, target groups and ordered gene axes passed. This is descriptive analysis after the experiment outcomes; it adds no promotion criterion and does not change any previous failure or selection decision. No metric was rerun and no reserved responses were accessed.

The shared gene gate is real-control arithmetic mean CPM >5, excluding each target’s own gene. CPM variances use sample ddof1 and full-axis per-cell normalization. Variance-trace ratios compare the sum of predicted gene variances with the sum of measured gene variances on that gate. Means may differ, so these are not pure dispersion effects.

| Response | Emitter | Seed | CPM variance trace ratio, mean (median) | Gene detection MAE | Detected genes/cell, prediction / real | Depth q05 / q50 / q95, prediction |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| adapter | current | 20261021 | 0.2870 (0.2912) | 0.08191 | 9647.7 / 8729.2 | 29458 / 54964 / 94326 |
| adapter | current | 20261022 | 0.2876 (0.2918) | 0.08017 | 9625.7 / 8729.2 | 28558 / 54521 / 94173 |
| adapter | preserving | 20261021 | 0.9906 (0.9991) | 0.02015 | 8724.1 / 8729.2 | 30280 / 54493 / 94792 |
| adapter | preserving | 20261022 | 0.9849 (0.9942) | 0.02029 | 8723.8 / 8729.2 | 30303 / 54504 / 95191 |
| current | current | 20261021 | 0.4378 (0.4077) | 0.08003 | 9621.2 / 8729.2 | 29458 / 54964 / 94326 |
| current | current | 20261022 | 0.4381 (0.4069) | 0.07839 | 9599.8 / 8729.2 | 28558 / 54521 / 94173 |
| current | preserving | 20261021 | 1.1874 (1.1601) | 0.02357 | 8621.4 / 8729.2 | 30280 / 54493 / 94792 |
| current | preserving | 20261022 | 1.1731 (1.1565) | 0.02365 | 8621.5 / 8729.2 | 30303 / 54504 / 95191 |

Measured depth quantiles averaged across targets are 29711 / 54050 / 93328. Depth table values average per-target quantiles; they are not pooled-cell quantiles. Detection MAE is a fraction, so0.01 is one percentage point.

The CSV (local artifact omitted: `redesign_distribution_audit.csv`) contains all400arm/seed/target rows, full depth quantiles, gene detection, CPM variance traces and supplementary raw-count variance traces. Raw-count variance is explicitly in squared-count units and includes depth effects. The JSON (local artifact omitted: `redesign_distribution_audit.json`) preserves artifact/receipt/freeze hashes, gate and normalization definitions, per-arm summaries and paired emitter differences.

This audit describes nonzero predictions, separately from the zero-effect identity/null audit. A more realistic distribution does not by itself establish improved DE directions, response specificity, expression error or the six-metric promotion gate.
