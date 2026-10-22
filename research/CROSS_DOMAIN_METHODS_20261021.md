# Relevant methods outside biology

The user's question about geography and astronomy fits three distinct parts of
the observed problem. These are hypotheses for adaptation, not evidence that an
external method will improve the challenge score.

| Domain | Primary method | Mapping to this project | Boundary |
|---|---|---|---|
| Geographic/survey sampling | Deville & Tillé, balanced sampling/cube method | Select a representative subset of real control cells while retaining known inclusion probabilities and raw cell variation. | Current P0 tests depth stratification, not the full cube method. No guarantee every gene improves. |
| Weather prediction | Schefzik et al., ensemble copula coupling | Calibrate marginal predictive distributions while preserving a multivariate rank/dependence template. | Integer counts, row-depth constraints, zero-effect identity and biological transfer require separate treatment. |
| Astronomy | Christiansen et al., Kepler injection–recovery | Inject known signals into real control backgrounds; measure signal recovery and false calls across signal strengths. | Synthetic recovery tests emission/detection, not whether the response predictor learns real biology. |

Primary sources checked 21 October 2026:

- [Efficient balanced sampling: the cube method](https://doc.rero.ch/record/296198/files/914893.pdf), Biometrika 91(4), 893–912 (2004).
- [Uncertainty quantification in complex simulation models using ensemble copula coupling](https://arxiv.org/abs/1302.7149), Schefzik, Thorarinsdottir & Gneiting (2013).
- [Measuring transit signal recovery in the Kepler pipeline II](https://www.ipac.caltech.edu/publication/2015ApJ...810...95C), Christiansen et al. (2015). Artificial transit signals were injected into real pixel observations and recovered through the operational pipeline.

Priority follows our measured failures: sampling precision first, then a known-signal
recovery diagnostic to separate emitter sensitivity from response-model weakness.
Dependence-preserving ensemble calibration is a later candidate once response
information improves. Existing experiment freezes are unchanged by these notes.

## Named benchmarks and practical priority

[WeatherBench 2](https://research.google/pubs/weatherbench-2-a-benchmark-for-the-next-generation-of-data-driven-global-weather-models/)
is a particularly close benchmark analogy in our assessment. Its official
[evaluation implementation](https://github.com/google-research/weatherbench2/blob/main/scripts/evaluate.py)
includes ensemble mean MSE, debiased ensemble mean MSE, ensemble variance and CRPS
components. This separates mean accuracy from distributional quality, directly
relevant to the tradeoff seen in our generator. These are similarities in the
evaluation problem, not proof that a weather model transfers to gene perturbations.

[PLAsTiCC](https://plasticc.org/) is an actual astronomical machine-learning
challenge, with classification as its goal. It supplies a challenge-design
analogy; Kepler's injection–recovery procedure is the more immediate diagnostic
idea for our count generator.

The completed P0 depth-stratification test reduces H1 bulk-profile sampling error
by only about1.1%, below its frozen5% requirement. Exact arithmetic-mean sampling
variance falls only1.6–2.5% across the four populations, even though depth-mean
variance falls by at least99.5%. Depth is therefore a weak balancing feature for
these normalized gene profiles. A future cube-method experiment could balance
several control-derived expression features as well as depth, while preserving
known inclusion probabilities and actual donor cells. It would require a new
freeze, mathematical tests and exact-policy null checks; it has not been run or
selected here. See [P0 results](../reports/sampling_precision_results.md).

There is also a direct model-design analogy to the P2 hurdle prototype:
[two-stage precipitation forecasting](https://arxiv.org/abs/0901.3484) models
occurrence separately from the amount conditional on occurrence. Our P2 design
likewise separates magnitude-proxy responsivity from conditional sign/magnitude.
It is a project-specific small model, not an implementation of that weather
paper's spatial Gaussian processes or precipitation distribution. The analogy
does not establish calibration or predictive improvement; the frozen source
benchmark remains the test.
