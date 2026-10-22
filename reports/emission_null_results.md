# Zero-response emission diagnostic

Prioritize a control-distribution-preserving zero-response emitter before further global response tuning.

| Arm | Median DE calls | Range | Detection MAE | Variance-trace ratio |
| --- | ---: | --- | ---: | ---: |
| measured_null | 0.0 | 0–0 | 0.01206 | 0.9939 |
| donor_replay | 0.0 | 0–0 | 0.01197 | 1.0075 |
| pool4_no_fit | 4672.0 | 4610–4700 | 0.08258 | 0.2848 |
| neutral_dual | 4323.5 | 4244–4378 | 0.08208 | 0.2844 |
| response_probe | 5603.0 | 5260–5914 | 0.07847 | 0.4707 |

Registered emitter-excess flag: **True**; ratio versus the larger empirical-null median: **4323.5×**; neutral exceeds donor replay in **10/10** aliases. The ratio uses the preregistered one-call floor because both empirical-null medians are zero; it is not an ordinary biological fold change.

Joint reallocation: 14/4000 cells; batch and guide margins remain exact.

## Overlap with existing response calls

- donor_replay: 0/55,696 response-probe calls shared (0.0%). This is diagnostic overlap, not attribution of biological false positives.
- pool4_no_fit: 30,733/55,696 response-probe calls shared (55.2%). This is diagnostic overlap, not attribution of biological false positives.
- neutral_dual: 29,933/55,696 response-probe calls shared (53.7%). This is diagnostic overlap, not attribution of biological false positives.

## Limits

- Ten aliases are correlated technical replicates with shared inference/comparison controls, not independent biological experiments.
- Real controls are independent of inference cells and comparison cells; batches/guides match marginally, with minimal recorded joint reallocation.
- Donor replay retains original donor depths; pooled arms use pooled depths. This compares full policies, not variance alone.
- Response-probe calls are not classified as false using null pseudo-outcomes; overlap with zero-response calls is only diagnostic.
- No real perturbation or reserved response outcomes were used to build this null reference.
