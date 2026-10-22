# One response-balanced hurdle objective: P2

**Magnitude proxies only; no DE-significance or deployment-probability claim.** The fixed prevalence operating point expresses class-balanced costs, not probability >0.5.

| Context | Head | Signed precision | Signed recovery | Call fraction | Unweighted MSE | Proxy Brier |
|---|---|---:|---:|---:|---:|---:|
| K562 | log2fc | 0.203106 | 0.096022 | 0.131342 | 0.106474 | 0.221956 |
| K562 | bulk_delta | 0.127662 | 0.231863 | 0.441011 | 0.018539 | 0.204956 |
| HCT116 | log2fc | 0.154536 | 0.493519 | 0.791464 | 0.282540 | 0.295058 |
| HCT116 | bulk_delta | 0.066412 | 0.340139 | 0.746369 | 0.024526 | 0.162251 |
| HEK293T | log2fc | 0.133914 | 0.493701 | 0.723819 | 0.238478 | 0.250880 |
| HEK293T | bulk_delta | 0.053727 | 0.400065 | 0.773815 | 0.021064 | 0.134030 |

Both-head advancement against matched fusion: **False**.

The JSON records fitting/recalibration weights and target IDs, original prevalences, monotone Platt parameters, proxy reliability bins, conditional sign/magnitude models, fixed decisions, paired uncertainty and failures. Recalibration uses a target-label subset within training contexts with full-source common preprocessing; it is not independent unseen-target calibration. The outer held-out context is absent from every training/calibration input and label.

The exact [P1 fusion gate](../docs/MATCHED_RESPONSE_BENCHMARK_PROTOCOL_20261021.md) is unchanged; stable false-call rate and matched-control CPM>5/<=5 strata remain diagnostic. [P2 protocol](../docs/MATCHED_RESPONSE_P2_PROTOCOL_20261021.md) was fixed before fitting. No method or threshold rescue is performed.

**Stop the next H1 scoring cycle.** This one objective did not establish sufficient source-only response evidence. Do not promote, extract reserved responses, or submit this candidate.
